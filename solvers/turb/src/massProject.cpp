/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "mass.hpp"
#include "timer.hpp"
#include <limits>

void mass_t::BoundarySetup(){

  //check all the bounaries for a Dirichlet
  allNeumann = 0;
  allNeumannPenalty = 0;

  //translate the mesh's element-to-boundaryflag mapping
  EToB.malloc(mesh.Nelements*mesh.Nfaces, 0);
  for (dlong e=0;e<mesh.Nelements;e++) {
    for (int f=0;f<mesh.Nfaces;f++) {
      int bc = mesh.EToB[e*mesh.Nfaces+f];
      if (bc>0) {
        int BC = BCType[bc];         //translate mesh's boundary flag
        EToB[e*mesh.Nfaces+f] = BC;  //record it
        if (BC!=2) allNeumann = 0;   //check if its a Dirchlet
      }
    }
  }
  o_EToB = platform.malloc<int>(EToB);

  //collect the allNeumann flags from other ranks
  mesh.comm.Allreduce(allNeumann, comm_t::Min);

  //translate the mesh's node-wise bc flag
  Nmasked = 0;
  mapB.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np, 0);
  for (int n=0;n<mesh.Nelements*mesh.Np;n++) {
    int bc = mesh.mapB[n];
    if (bc>0) {
      int BC = BCType[bc];     //translate mesh's boundary flag
      mapB[n] = BC;  //record it

      if (mapB[n] == 1) Nmasked++;   //Dirichlet boundary
    }
  }
  o_mapB = platform.malloc<int>(mapB);

  maskIds.malloc(Nmasked);
  Nmasked =0; //reset
  for (dlong n=0;n<mesh.Nelements*mesh.Np;n++) {
    if (mapB[n] == 1) maskIds[Nmasked++] = n;
  }
  o_maskIds = platform.malloc<int>(maskIds);

  //make a masked version of the global id numbering
  maskedGlobalIds.malloc(mesh.Nelements*mesh.Np);
  maskedGlobalIds.copyFrom(mesh.globalIds);
  for (dlong n=0;n<Nmasked;n++) {
    maskedGlobalIds[maskIds[n]] = 0;
  }

  //use the masked ids to make another gs handle (signed so the gather is defined)
  bool verbose = false;
  bool unique = true; //flag a unique node in every gather node
  ogsMasked.Setup(mesh.Nelements*mesh.Np, maskedGlobalIds,
                  mesh.comm, ogs::Signed, ogs::Auto,
                  unique, verbose, platform);

  // TW - need Nfields  ?
  //setup normalization constant
  //note that we can use the mesh ogs, since there are no masked nodes
  allNeumannScale = 0./sqrt((dfloat)ogsMasked.NgatherGlobal);

  /* use the masked gs handle to define a global ordering */
  dlong Ntotal  = mesh.Np*mesh.Nelements; // number of degrees of freedom on this rank (before gathering)
  hlong Ngather = ogsMasked.Ngather;     // number of degrees of freedom on this rank (after gathering)

  // build inverse degree vectors
  // used for the weight in linear solvers (used in C0)
  weight.malloc(Ntotal, (pfloat)1.0);

  weightG.malloc(Ngather);
  ogsMasked.Gather(weightG, weight, 1, ogs::Add, ogs::Trans);

  for(dlong n=0;n<ogsMasked.Ngather;++n) {
    if (weightG[n]>0.0) weightG[n] = 1./weightG[n];
  }

  ogsMasked.Scatter(weight, weightG, 1, ogs::NoTrans);

  // TW
  o_weight  = platform.malloc<pfloat>(weight);
  o_weightG = platform.malloc<pfloat>(weightG);

  // create a global numbering system
  memory<hlong> globalIds(Ngather);

  // every gathered degree of freedom has its own global id
  hlong globalOffset=static_cast<hlong>(Ngather);
  comm.Scan(Ngather, globalOffset);
  globalOffset = globalOffset-Ngather;

  //use the offsets to set a consecutive global numbering
  for (dlong n =0;n<ogsMasked.Ngather;n++) {
    globalIds[n] = n + globalOffset;
  }

  //scatter this numbering to the original nodes
  maskedGlobalNumbering.malloc(Ntotal, -1);
  ogsMasked.Scatter(maskedGlobalNumbering, globalIds, 1, ogs::NoTrans);

  /* Build halo exchange for gathered ordering */
  gHalo.SetupFromGather(ogsMasked);

  GlobalToLocal.malloc(mesh.Nelements*mesh.Np);
  ogsMasked.SetupGlobalToLocalMapping(GlobalToLocal);

  o_GlobalToLocal = platform.malloc<dlong>(GlobalToLocal);
}

// TW need to replace with mass and non-constant viscosity

void mass_t::BuildOperatorDiagonal(memory<dfloat>& diagA){

  if(comm_t::world().rank()==0) {printf("Building diagonal...");fflush(stdout);}

  // assume C0
  memory<dfloat> diagAL(mesh.Np*mesh.Nelements*Nfields);
  
  switch(mesh.elementType){
  case Mesh::TRIANGLES:
    BuildOperatorDiagonalContinuousTri2D(diagAL);
    break;
  case Mesh::QUADRILATERALS:
    BuildOperatorDiagonalContinuousQuad2D(diagAL);
    break;
  case Mesh::TETRAHEDRA:
    BuildOperatorDiagonalContinuousTet3D(diagAL);
    break;
  case Mesh::HEXAHEDRA:
    BuildOperatorDiagonalContinuousHex3D(diagAL);
    break;
  }
  
  //gather the diagonal to assemble it
  ogsMasked.Gather(diagA, diagAL, Nfields, ogs::Add, ogs::Trans);

  if(comm_t::world().rank()==0) printf("done.\n");
}

void mass_t::BuildOperatorDiagonal(deviceMemory<pfloat> &o_invDiagA ){

  deviceMemory<dfloat> o_diagAL =
    platform.reserve<dfloat>(mesh.Nelements*Nfields*mesh.Np);
  deviceMemory<dfloat> o_diagA =
    platform.reserve<dfloat>(Ndofs);
  
  buildOperatorDiagonalKernel(mesh.Nelements,
			      mesh.o_wJ,
			      mesh.o_MM,
			      o_diagAL);
  
  ogsMasked.Gather(o_diagA, o_diagAL, Nfields, ogs::Add, ogs::Trans);

  reciprocalKernel(Ndofs, o_diagA, o_invDiagA);

  if(0){
    memory<pfloat> tmp(Ndofs);
    o_invDiagA.copyTo(tmp);
    memory<dfloat> MM(mesh.Np*mesh.Np);
    mesh.o_MM.copyTo(MM);
    for(int n=0;n<mesh.Np;++n){
      for(int m=0;m<mesh.Np;++m){
	printf("%g, ", MM[n*mesh.Np+m]);
      }
      printf("\n");
    }
    
  }
  
}



void mass_t::BuildOperatorDiagonalContinuousTri2D(memory<dfloat>& A) {
  exit(-1);
  
  for(dlong eM=0;eM<mesh.Nelements;++eM){
    dfloat J   = mesh.wJ[eM];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      A[Nfields*(eM*mesh.Np+n)+0] = J*mesh.MM[n+n*mesh.Np];
      A[Nfields*(eM*mesh.Np+n)+1] = J*mesh.MM[n+n*mesh.Np];
      A[Nfields*(eM*mesh.Np+n)+2] = J*mesh.MM[n+n*mesh.Np];
      A[Nfields*(eM*mesh.Np+n)+3] = J*mesh.MM[n+n*mesh.Np];
    }
  }

  
}

void mass_t::BuildOperatorDiagonalContinuousQuad2D(memory<dfloat>& A) {

  /*
    \sum_j d_j ( nu (d_j u_i + d_i u_j ))
    => -\sum_j (d_j phi, nu (d_j u_i + d_i u_j ))
    
    [2.*Dx'*nu*Dx*u + Dy'*nu*Dy*u] + Dy'*nu*Dx*v
    Dx'*nu*Dy*u + [Dx'*nu*Dx*v + 2.*Dy'*nu*Dy*v]
    
  */
  
  for(dlong e=0;e<mesh.Nelements;++e){
    for (int m=0;m<mesh.Nq;++m) {
      for (int n=0;n<mesh.Nq;++n) {
        dlong iid = n+m*mesh.Nq;
	dlong lid = iid + e*mesh.Np;
	dlong uid = Nfields*lid + 0;
	dlong vid = Nfields*lid + 1;
	dlong kid = Nfields*lid + 2;
	dlong tauid = Nfields*lid + 3;

	dlong vbase = e*mesh.Np*mesh.Nvgeo;
	dfloat JW = mesh.vgeo[vbase + n + m*mesh.Nq + mesh.JWID*mesh.Np];
	A[uid] = JW;
	A[vid] = JW;
	A[kid] = JW;
	A[tauid] = JW;
      }
    }
  }
}

void mass_t::BuildOperatorDiagonalContinuousTet3D(memory<dfloat>& A) {
  
  for(dlong eM=0;eM<mesh.Nelements;++eM){
    dfloat J   = mesh.wJ[eM];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      dlong lid = n + eM*mesh.Np;
      dlong uid = Nfields*lid + 0;
      dlong vid = Nfields*lid + 1;
      dlong wid = Nfields*lid + 2;
      
      A[Nfields*(eM*mesh.Np+n)+0] = J*mesh.MM[n+n*mesh.Np];
      A[Nfields*(eM*mesh.Np+n)+1] = J*mesh.MM[n+n*mesh.Np];
      A[Nfields*(eM*mesh.Np+n)+2] = J*mesh.MM[n+n*mesh.Np]; // fix for k-tau later      
    }
  }
}

void mass_t::BuildOperatorDiagonalContinuousHex3D(memory<dfloat>& A) {

  std::cout << "BuildOperatorDiagonalContinuousHex3D not implemented" << std::endl;
  exit(-1);
}

void mass_t::Operator(deviceMemory<double> &o_q, deviceMemory<double> &o_Aq){

  deviceMemory<double> o_MM, o_JW, o_vgeo;

  if constexpr (std::is_same_v<dfloat,double>) {
    o_MM   = mesh.o_MM;
    o_JW   = mesh.o_wJ;
  } else if (std::is_same_v<pfloat,double>) {
    o_MM   = mesh.o_pfloat_MM;
    o_JW   = mesh.o_pfloat_wJ;
  } else {
    LIBP_FORCE_ABORT("mass_t::Operator called on type double, but double not set in types.h");
  }

  // assume C0
  //buffer for local Ax
  deviceMemory<double> o_AqL = platform.reserve<double>(Nfields*mesh.Np*mesh.Nelements);

  gHalo.ExchangeStart(o_q, 1); // Nfields);
  
  if(mesh.NlocalGatherElements/2){
    massPartialAxKernel(mesh.NlocalGatherElements/2,
			mesh.o_localGatherElementList,
			o_GlobalToLocal,
			o_JW,
			o_MM,
			o_q,
			o_AqL);
  }
  
  // finalize halo exchange
  gHalo.ExchangeFinish(o_q, 1); // Nfields);
  
  if(mesh.NglobalGatherElements) {
    
    massPartialAxKernel(mesh.NglobalGatherElements,
		    mesh.o_globalGatherElementList,
		    o_GlobalToLocal,
		    o_JW,
		    o_MM,
		    o_q,
		    o_AqL);
  }

  //gather result to Aq
  ogsMasked.GatherStart(o_Aq, o_AqL, Nfields, ogs::Add, ogs::Trans);

  if((mesh.NlocalGatherElements+1)/2){
    massPartialAxKernel((mesh.NlocalGatherElements+1)/2,
			mesh.o_localGatherElementList+(mesh.NlocalGatherElements/2),
			o_GlobalToLocal,
			o_JW,
			o_MM,
			o_q,
			o_AqL);
  }

  ogsMasked.GatherFinish(o_Aq, o_AqL, Nfields, ogs::Add, ogs::Trans);
  
}


void mass_t::Operator(deviceMemory<float> &o_q, deviceMemory<float> &o_Aq){

  deviceMemory<float> o_MM, o_JW;
  
  if constexpr (std::is_same_v<dfloat,float>) {
    o_MM   = mesh.o_MM;
    o_JW   = mesh.o_wJ;
  } else if (std::is_same_v<pfloat,float>) {
    o_MM   = mesh.o_pfloat_MM;
    o_JW   = mesh.o_pfloat_wJ;
  } else {
    LIBP_FORCE_ABORT("mass_t::Operator called on type float, but float not set in types.h");
  }

  // assume C0
  //buffer for local Ax
  deviceMemory<float> o_AqL = platform.reserve<float>(Nfields*mesh.Np*mesh.Nelements);
  
  gHalo.ExchangeStart(o_q, 1); // Nfields);
  
  if(mesh.NlocalGatherElements/2){
    floatMassPartialAxKernel(mesh.NlocalGatherElements/2,
			     mesh.o_localGatherElementList,
			     o_GlobalToLocal,
			     o_JW,
			     o_MM,
			     o_q,
			     o_AqL);
  }
  
  // finalize halo exchange
  gHalo.ExchangeFinish(o_q, 1); // Nfields);
  
  if(mesh.NglobalGatherElements) {
    floatMassPartialAxKernel(mesh.NglobalGatherElements,
			     mesh.o_globalGatherElementList,
			     o_GlobalToLocal,
			     o_JW,
			     o_MM,
			     o_q,
			     o_AqL);
  }
  
  //gather result to Aq
  ogsMasked.GatherStart(o_Aq, o_AqL, Nfields, ogs::Add, ogs::Trans);
  
  if((mesh.NlocalGatherElements+1)/2){
    floatMassPartialAxKernel((mesh.NlocalGatherElements+1)/2,
			     mesh.o_localGatherElementList+(mesh.NlocalGatherElements/2),
			     o_GlobalToLocal,
			     o_JW,
			     o_MM,
			     o_q,
			     o_AqL);
  }
  
  ogsMasked.GatherFinish(o_Aq, o_AqL, Nfields, ogs::Add, ogs::Trans);
  
}


void mass_t::BlockInverseOperator(deviceMemory<double> &o_q, deviceMemory<double> &o_Aq){

  deviceMemory<double> o_tmp_invMM, o_tmp_invJW;

  if constexpr (std::is_same_v<dfloat,double>) {
    o_tmp_invMM   = o_invMM;
    o_tmp_invJW   = o_invJW;
  } else if (std::is_same_v<pfloat,double>) {
    o_tmp_invMM   = o_pfloat_invMM;
    o_tmp_invJW   = o_pfloat_invJW;
  } else {
    LIBP_FORCE_ABORT("mass_t::BlockInverseOperator called on type double, but double not set in types.h");
  }

  // assume C0
  //buffer for local Ax
  deviceMemory<double> o_AqL = platform.reserve<double>(Nfields*mesh.Np*mesh.Nelements);

  gHalo.ExchangeStart(o_q, 1); // Nfields);
  
  if(mesh.NlocalGatherElements/2){
    massPartialAxKernel(mesh.NlocalGatherElements/2,
			mesh.o_localGatherElementList,
			o_GlobalToLocal,
			o_tmp_invJW,
			o_tmp_invMM,
			o_q,
			o_AqL);
  }
  
  // finalize halo exchange
  gHalo.ExchangeFinish(o_q, 1); // Nfields);
  
  if(mesh.NglobalGatherElements) {
    
    massPartialAxKernel(mesh.NglobalGatherElements,
		    mesh.o_globalGatherElementList,
		    o_GlobalToLocal,
		    o_tmp_invJW,
		    o_tmp_invMM,
		    o_q,
		    o_AqL);
  }

  //gather result to Aq
  ogsMasked.GatherStart(o_Aq, o_AqL, Nfields, ogs::Add, ogs::Trans);

  if((mesh.NlocalGatherElements+1)/2){
    massPartialAxKernel((mesh.NlocalGatherElements+1)/2,
			mesh.o_localGatherElementList+(mesh.NlocalGatherElements/2),
			o_GlobalToLocal,
			o_tmp_invJW,
			o_tmp_invMM,
			o_q,
			o_AqL);
  }

  ogsMasked.GatherFinish(o_Aq, o_AqL, Nfields, ogs::Add, ogs::Trans);
  
}


void mass_t::BlockInverseOperator(deviceMemory<float> &o_q, deviceMemory<float> &o_Aq){

  deviceMemory<float> o_tmp_invMM, o_tmp_invJW;

  if constexpr (std::is_same_v<dfloat,float>) {
    o_tmp_invMM   = o_invMM;
    o_tmp_invJW   = o_invJW;
  } else if (std::is_same_v<pfloat,float>) {
    o_tmp_invMM   = o_pfloat_invMM;
    o_tmp_invJW   = o_pfloat_invJW;
  } else {
    LIBP_FORCE_ABORT("mass_t::BlockInverseOperator called on type float, but float not set in types.h");
  }

  // assume C0
  //buffer for local Ax
  deviceMemory<float> o_AqL = platform.reserve<float>(Nfields*mesh.Np*mesh.Nelements);
  
  gHalo.ExchangeStart(o_q, 1); // Nfields);
  
  if(mesh.NlocalGatherElements/2){
    floatMassPartialAxKernel(mesh.NlocalGatherElements/2,
			 mesh.o_localGatherElementList,
			 o_GlobalToLocal,
			 o_tmp_invJW,
			 o_tmp_invMM,
			 o_q,
			 o_AqL);
  }
  
  // finalize halo exchange
  gHalo.ExchangeFinish(o_q, 1); // Nfields);
  
  if(mesh.NglobalGatherElements) {
    floatMassPartialAxKernel(mesh.NglobalGatherElements,
			 mesh.o_globalGatherElementList,
			 o_GlobalToLocal,
			 o_tmp_invJW,
			 o_tmp_invMM,
			 o_q,
			 o_AqL);
  }
  
  //gather result to Aq
  ogsMasked.GatherStart(o_Aq, o_AqL, Nfields, ogs::Add, ogs::Trans);
  
  if((mesh.NlocalGatherElements+1)/2){
    floatMassPartialAxKernel((mesh.NlocalGatherElements+1)/2,
			 mesh.o_localGatherElementList+(mesh.NlocalGatherElements/2),
			 o_GlobalToLocal,
			 o_tmp_invJW,
			 o_tmp_invMM,
			 o_q,
			 o_AqL);
  }
  
  ogsMasked.GatherFinish(o_Aq, o_AqL, Nfields, ogs::Add, ogs::Trans);
  
}



void MassJacobiPrecon::Update(){
  mass.BuildOperatorDiagonal(o_invDiagA);
}


int mass_t::Solve(linearSolver_t<dfloat>& linearSolver,
		  deviceMemory<dfloat> &o_x,
		  deviceMemory<dfloat> &o_r,
		  const dfloat tol, const int MAXIT, const int verbose){

  int Niter = linearSolver.Solve(*this, precon, o_x, o_r, tol, MAXIT, verbose);

  return Niter;
}



// Jacobi preconditioner
MassJacobiPrecon::MassJacobiPrecon(mass_t& _mass):
  mass(_mass) {

  o_invDiagA = mass.platform.malloc<pfloat>(mass.Ndofs);
  
  mass.BuildOperatorDiagonal(o_invDiagA);
  
}

void MassJacobiPrecon::Operator(deviceMemory<pfloat>& o_r, deviceMemory<pfloat>& o_Mr) {

  linAlg_t& linAlg = mass.platform.linAlg();

  // Mr = invDiag.*r
  linAlg.amxpy(mass.Ndofs, (pfloat)1.0, o_invDiagA, o_r, (pfloat)0.0, o_Mr);

}

void mass_t::Run(){

  //setup linear algebra module
  platform.linAlg().InitKernels({"set"});

  //setup linear solver
  hlong NglobalDofs;
  NglobalDofs = ogsMasked.NgatherGlobal*Nfields;

  std::cout << "MASS LINEARSOLVER" << std::endl;
  
  linearSolver_t<dfloat> linearSolver;
  linearSolver.Setup<LinearSolver::pcg<dfloat> >(Ndofs, Nhalo, platform, settings, comm);

  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add standard boundary functions
  int Nmax = std::max(mesh.Np, mesh.Nfaces*mesh.Nfp);
  kernelInfo["defines/" "p_Nmax"]= Nmax;
  kernelInfo["defines/" "p_Nfields"]= Nfields;
  kernelInfo["defines/" "p_NVfields"]= Nfields;

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();
  
  std::string oklFilePrefix = DINS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  fileName   = oklFilePrefix + "massRhs" + suffix + oklFileSuffix;
  kernelName = "massRhs" + suffix;
  kernel_t forcingKernel = platform.buildKernel(fileName, kernelName,
                                                    kernelInfo);

  kernel_t rhsBCKernel, addBCKernel;
  fileName   = oklFilePrefix + "massRhsBC" + suffix + oklFileSuffix;
  kernelName = "massRhsBC" + suffix;
  
  rhsBCKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
  
  fileName   = oklFilePrefix + "massAddBC" + suffix + oklFileSuffix;
  kernelName = "massAddBC" + suffix;
  
  addBCKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  //create occa buffers
  dlong Nall = Nfields*mesh.Np*(mesh.Nelements+mesh.totalHaloPairs);
  memory<dfloat> rL(Nall);
  memory<dfloat> xL(Nall);
  deviceMemory<dfloat> o_rL = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_xL = platform.malloc<dfloat>(Nall);

  deviceMemory<dfloat> o_r, o_x;
  dlong Ng = ogsMasked.Ngather;
  dlong Nghalo = gHalo.Nhalo;
  dlong Ngall = Nfields*(Ng + Nghalo);
  o_r = platform.malloc<dfloat>(Ngall);
  o_x = platform.malloc<dfloat>(Ngall);

  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  //populate rhs forcing
  forcingKernel(mesh.Nelements,
                mesh.o_wJ,
                mesh.o_MM,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                o_rL);

  //Set x to zero
  platform.linAlg().set(mesh.Nelements*mesh.Np*Nfields, (dfloat)0.0, o_xL);

  rhsBCKernel(mesh.Nelements,
	      mesh.o_wJ,
	      mesh.o_x,
	      mesh.o_y,
	      mesh.o_z,
	      o_rL);

  // gather rhs to globalDofs if c0
  ogsMasked.Gather(o_r, o_rL, Nfields, ogs::Add, ogs::Trans);
  ogsMasked.Gather(o_x, o_xL, Nfields, ogs::Add, ogs::NoTrans);

  int maxIter = 5000;

  timePoint_t start = GlobalPlatformTime(platform);
  
  //call the solver
  dfloat tol = (sizeof(dfloat)==sizeof(double)) ? 1.0e-8 : 1.0e-5;

  bool verbose=false;
  int iter = Solve(linearSolver, o_x, o_r, tol, maxIter, verbose);

  //add the boundary data to the masked nodes
  // scatter x to LocalDofs if c0
  ogsMasked.Scatter(o_xL, o_x, Nfields, ogs::NoTrans);

  //fill masked nodes with BC data
  addBCKernel(mesh.Nelements,
	      mesh.o_x,
	      mesh.o_y,
	      mesh.o_z,
	      o_mapB,
	      o_xL);
  
  timePoint_t end = GlobalPlatformTime(platform);
  double elapsedTime = ElapsedTime(start, end);

  if ((mesh.rank==0) && verbose){
    printf("%d, " hlongFormat ", %g, %d, %g, %g; global: N, dofs, elapsed, iterations, time per node, nodes*iterations/time %s\n",
           mesh.N,
           NglobalDofs,
           elapsedTime,
           iter,
           elapsedTime/(NglobalDofs),
           NglobalDofs*((dfloat)iter/elapsedTime), "JACOBI");
  }


  // output norm of final solution
  { // NEED TO FIX MASS MATRIX
    //compute q.M*q
    dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
    deviceMemory<dfloat> o_MxL = platform.reserve<dfloat>(Nentries);
    mesh.MassMatrixApply(o_xL, o_MxL);

    dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_xL, o_MxL, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);
  }
}

void mass_t::Setup(platform_t& _platform, mesh_t& _mesh, settings_t& _settings,
		   const int _NBCTypes, const memory<int> _BCType){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  //  settings.report();
  
  Nfields = mesh.dim + 2; // velocity + k-tau

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);
  ogs::InitializeKernels(platform, ogs::Pfloat, ogs::Add);

  //setup linear algebra module
  platform.linAlg().InitKernels({"add", "sum", "scale",
        "axpy", "zaxpy",
        "amx", "amxpy", "zamxpy",
        "adx", "adxpy", "zadxpy",
        "innerProd", "norm2", "d2p", "p2d"});

  /*setup trace halo exchange */
  traceHalo = mesh.HaloTraceSetup(Nfields);

  // Boundary Type translation. Just defaults.
  NBCTypes = _NBCTypes;
  BCType.malloc(NBCTypes);
  BCType.copyFrom(_BCType);

  std::cout << "BOUNDARY SETUP" << std::endl;
  
  //setup boundary flags and make mask and masked ogs
  BoundarySetup();

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();

  std::string oklFilePrefix = DINS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  //add standard boundary functions

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 1024;

  kernelInfo["defines/" "p_Nfields"]= Nfields;
  kernelInfo["defines/" "p_NVfields"]= Nfields;
  
  int NblockV = std::max(1,blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  properties_t kernelInfoDouble = kernelInfo;
  kernelInfoDouble["defines/dfloat"] = "double";
  kernelInfoDouble["defines/dfloat4"] = "double4";

  properties_t kernelInfoFloat = kernelInfo;
  kernelInfoFloat["defines/dfloat"] = "float";
  kernelInfoFloat["defines/dfloat4"] = "float4";

  // Ax kernel (assume C0)
  fileName   = oklFilePrefix + "massKernels" + suffix + oklFileSuffix;
  kernelName = "massAx" + suffix;
  
  massAxKernel = platform.buildKernel(fileName, kernelName, kernelInfoDouble);
  floatMassAxKernel = platform.buildKernel(fileName, kernelName, kernelInfoFloat);

  // Ax kernel (assume C0)
  kernelName = "massPartialAx" + suffix;
  
  massPartialAxKernel = platform.buildKernel(fileName, kernelName, kernelInfoDouble);
  floatMassPartialAxKernel = platform.buildKernel(fileName, kernelName, kernelInfoFloat);

  kernelName = "massScatter" + suffix;
  massScatterKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "massWeight" + suffix;
  weightKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
  
  
  /* Preconditioner Setup */
  Ndofs = ogsMasked.Ngather*Nfields;
  Nhalo = gHalo.Nhalo*Nfields;

  kernelName = "massBuildOperatorDiagonal" + suffix;

  buildOperatorDiagonalKernel = platform.buildKernel(fileName, kernelName,
						     kernelInfo);

  // diagonal inverse (dfloat=>(pfloat)(1/float))
  kernelName = "massReciprocal";

  reciprocalKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

#if 0
  // TET
  memory<dfloat> invMM(mesh.MM);
  memory<pfloat> pfloat_invMM(mesh.Np*mesh.Np);
  memory<pfloat> pfloat_invJW(mesh.Nelements);

  linAlg_t::matrixInverse(mesh.Np, invMM);
  
  for(int n=0;n<mesh.Np;++n)
    for(int m=0;m<mesh.Np;++m)
      pfloat_invMM[n*mesh.Np+m] = invMM[n*mesh.Np+m];

  for(dlong e=0;e<mesh.Nelements;++e){
    pfloat_invJW[e] = 1./mesh.wJ[e];
  }
  
  o_pfloat_invJW = platform.malloc<pfloat>(mesh.Nelements, pfloat_invJW);
  o_pfloat_invMM = platform.malloc<pfloat>(mesh.Np*mesh.Np, pfloat_invMM);
#endif

  
  // assume Jacobi

  precon.Setup<MassJacobiPrecon>(*this);
  if(0)
    precon.Setup<MassInversePrecon>(*this);
  
  
}


// block mass
MassInversePrecon::MassInversePrecon(mass_t& _mass):
  mass(_mass) {
  
}

void MassInversePrecon::Operator(deviceMemory<pfloat>& o_r, deviceMemory<pfloat>& o_Mr) {

  linAlg_t& linAlg = mass.platform.linAlg();
  
  mass.BlockInverseOperator(o_r, o_Mr);

  dlong Ngather = mass.ogsMasked.Ngather;     // number of degrees of freedom on this rank (after gathering)

  mass.weightKernel(Ngather, mass.o_weightG, o_Mr);
}
