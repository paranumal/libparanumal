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

#include "stress.hpp"
#include "timer.hpp"
#include <limits>

void stress_t::BoundarySetup(){

  //check all the bounaries for a Dirichlet
  allNeumann = (lambda==0) ? 1 : 0; //if lambda>0 we don't care about all Neumann problem
  allNeumannPenalty = 1.;

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
  bool verbose = settings.compareSetting("VERBOSE", "TRUE") ? true : false;
  bool unique = true; //flag a unique node in every gather node
  ogsMasked.Setup(mesh.Nelements*mesh.Np, maskedGlobalIds,
                  mesh.comm, ogs::Signed, ogs::Auto,
                  unique, verbose, platform);

  // TW - need Nfields ?
  //setup normalization constant
  if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    allNeumannScale = 1./sqrt((dfloat)mesh.Np*mesh.NelementsGlobal);
  } else {
    //note that we can use the mesh ogs, since there are no masked nodes
    allNeumannScale = 1./sqrt((dfloat)ogsMasked.NgatherGlobal);
  }

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


void stress_t::Operator(deviceMemory<double> &o_q, deviceMemory<double> &o_Aq){

  deviceMemory<double> o_MM, o_D, o_S;
  deviceMemory<double> o_wJ, o_vgeo;

  if constexpr (std::is_same_v<dfloat,double>) {
    o_MM   = mesh.o_MM;
    o_D    = mesh.o_D;
    o_S    = mesh.o_S;

    o_wJ   = mesh.o_wJ;
    o_vgeo = mesh.o_vgeo;
    
  } else if (std::is_same_v<pfloat,double>) {
    o_MM   = mesh.o_pfloat_MM;
    o_D    = mesh.o_pfloat_D;
    o_S    = mesh.o_pfloat_S;

    o_wJ   = mesh.o_pfloat_wJ;
    o_vgeo = mesh.o_pfloat_vgeo;
  } else {
    LIBP_FORCE_ABORT("stress_t::Operator called on type double, but double not set in types.h");
  }

  // assume C0
  //buffer for local Ax
  deviceMemory<double> o_AqL = platform.reserve<double>(Nfields*mesh.Np*mesh.Nelements);

  gHalo.ExchangeStart(o_q, Nfields);
  
  if(mesh.NlocalGatherElements/2){
    partialAxKernel(mesh.NlocalGatherElements/2,
		    mesh.o_localGatherElementList,
		    o_GlobalToLocal,
		    o_nut,
		    o_wJ,
		    o_vgeo,
		    o_D,
		    o_S,
		    o_MM,
		    static_cast<double>(lambda),
		    o_q,
		    o_AqL);
  }
  
  // finalize halo exchange
  gHalo.ExchangeFinish(o_q, Nfields);
  
  if(mesh.NglobalGatherElements) {
    
    partialAxKernel(mesh.NglobalGatherElements,
		    mesh.o_globalGatherElementList,
		    o_GlobalToLocal,
		    o_nut,
		    o_wJ,
		    o_vgeo,
		    o_D,
		    o_S,
		    o_MM,
		    static_cast<double>(lambda),
		    o_q,
		    o_AqL);
  }

  //gather result to Aq
  ogsMasked.GatherStart(o_Aq, o_AqL, Nfields, ogs::Add, ogs::Trans);

  if((mesh.NlocalGatherElements+1)/2){
    partialAxKernel((mesh.NlocalGatherElements+1)/2,
		    mesh.o_localGatherElementList+(mesh.NlocalGatherElements/2),
		    o_GlobalToLocal,
		    o_nut,
		    o_wJ,
		    o_vgeo,
		    o_D,
		    o_S,
		    o_MM,
		    static_cast<double>(lambda),
		    o_q,
		    o_AqL);
  }

  ogsMasked.GatherFinish(o_Aq, o_AqL, Nfields, ogs::Add, ogs::Trans);

}


void stress_t::Operator(deviceMemory<float> &o_q, deviceMemory<float> &o_Aq){
  // verifies that we do not need pfloat o_nut
  printf("stress_t::Operator float\n");
  exit(-1);
  
  deviceMemory<float> o_MM, o_D, o_S;
  deviceMemory<float> o_wJ, o_vgeo;

  if constexpr (std::is_same_v<dfloat,float>) {
    o_MM   = mesh.o_MM;
    o_D    = mesh.o_D;
    o_S    = mesh.o_S;

    o_wJ   = mesh.o_wJ;
    o_vgeo = mesh.o_vgeo;
  } else if (std::is_same_v<pfloat,float>) {
    o_MM   = mesh.o_pfloat_MM;
    o_D    = mesh.o_pfloat_D;
    o_S    = mesh.o_pfloat_S;

    o_wJ   = mesh.o_pfloat_wJ;
    o_vgeo = mesh.o_pfloat_vgeo;
  } else {
    LIBP_FORCE_ABORT("stress_t::Operator called on type float, but float not set in types.h");
  }

  // assume C0
  //buffer for local Ax
  deviceMemory<float> o_AqL = platform.reserve<float>(Nfields*mesh.Np*mesh.Nelements);
  
  gHalo.ExchangeStart(o_q, Nfields);
  
  if(mesh.NlocalGatherElements/2){
    floatPartialAxKernel(mesh.NlocalGatherElements/2,
			 mesh.o_localGatherElementList,
			 o_GlobalToLocal,
			 o_nut,
			 o_wJ,
			 o_vgeo,
			 o_D,
			 o_S,
			 o_MM,
			 static_cast<float>(lambda),
			 o_q,
			 o_AqL);
  }
  
  // finalize halo exchange
  gHalo.ExchangeFinish(o_q, Nfields);
  
  if(mesh.NglobalGatherElements) {
    floatPartialAxKernel(mesh.NglobalGatherElements,
			 mesh.o_globalGatherElementList,
			 o_GlobalToLocal,
			 o_nut,
			 o_wJ,
			 o_vgeo,
			 o_D,
			 o_S,
			 o_MM,
			 static_cast<float>(lambda),
			 o_q,
			 o_AqL);
  }
  
  //gather result to Aq
  ogsMasked.GatherStart(o_Aq, o_AqL, Nfields, ogs::Add, ogs::Trans);
  
  if((mesh.NlocalGatherElements+1)/2){
    floatPartialAxKernel((mesh.NlocalGatherElements+1)/2,
			 mesh.o_localGatherElementList+(mesh.NlocalGatherElements/2),
			 o_GlobalToLocal,
			 o_nut,
			 o_wJ,
			 o_vgeo,
			 o_D,
			 o_S,
			 o_MM,
			 static_cast<float>(lambda),
			 o_q,
			 o_AqL);
  }
  
  ogsMasked.GatherFinish(o_Aq, o_AqL, Nfields, ogs::Add, ogs::Trans);
  
}

int stress_t::Solve(linearSolver_t<dfloat>& linearSolver,
                      deviceMemory<dfloat> &o_x,
		      deviceMemory<dfloat> &o_r,
                      const dfloat tol, const int MAXIT, const int verbose){

  // if there is a nullspace, remove the constant vector from r
  if(allNeumann) ZeroMean(o_r);
  
  int Niter = linearSolver.Solve(*this, precon, o_x, o_r, tol, MAXIT, verbose);

  return Niter;
}


void stress_t::BuildOperatorDiagonal(memory<dfloat>& diagA){
  printf("Wrong stress_t::BuildOperatorDiagonal\n");
  exit(-1);
  
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

void stress_t::BuildOperatorDiagonal(deviceMemory<pfloat> &o_invDiagA ){

  deviceMemory<dfloat> o_diagAL =
    platform.reserve<dfloat>(mesh.Nelements*Nfields*mesh.Np);
  deviceMemory<dfloat> o_diagA =
    platform.reserve<dfloat>(Ndofs);
  
  dfloat neumannBoost = (allNeumann) ? allNeumannPenalty*
    allNeumannScale*allNeumannScale: 0;

  buildOperatorDiagonalKernel(mesh.Nelements,
				     o_nut,
				     mesh.o_mapB,
				     neumannBoost,
				     mesh.o_wJ,
				     mesh.o_vgeo,
				     mesh.o_D,
				     mesh.o_S,
				     mesh.o_MM,
				     lambda,
				     o_diagAL);
  
  ogsMasked.Gather(o_diagA, o_diagAL, Nfields, ogs::Add, ogs::Trans);

  reciprocalKernel(Ndofs, o_diagA, o_invDiagA);
}



void stress_t::BuildOperatorDiagonalContinuousTri2D(memory<dfloat>& A) {
  
  std::cout << "BuildOperatorDiagonalContinuousTri2D not implemented" << std::endl;
  exit(-1);
}

void stress_t::BuildOperatorDiagonalContinuousQuad2D(memory<dfloat>& A) {

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
	dlong uid = 2*lid + 0;
	dlong vid = 2*lid + 1;

	dfloat fac = 2; 
	
        if (mapB[n+m*mesh.Nq+e*mesh.Np]!=1) {

	  dlong vbase = e*mesh.Np*mesh.Nvgeo;
          A[uid] = 0;
	  A[vid] = 0;

          for (int k=0;k<mesh.Nq;k++) {
            int id = k+m*mesh.Nq;
	    dfloat Dkn = mesh.D[n+k*mesh.Nq];
	    dfloat rx = mesh.vgeo[vbase + id + mesh.RXID*mesh.Np];
	    dfloat ry = mesh.vgeo[vbase + id + mesh.RYID*mesh.Np];
	    dfloat wJ = mesh.vgeo[vbase + id + mesh.JWID*mesh.Np];	    
	    dfloat nut_km = nut[e*mesh.Np+id];

	    dfloat uGrr = (fac*rx*rx + ry*ry)*nut_km*wJ;
            A[uid] += uGrr*Dkn*Dkn; // strided for gather
	    
	    dfloat vGrr = (fac*ry*ry + rx*rx)*nut_km*wJ;
	    A[vid] += vGrr*Dkn*Dkn;
          }


          for (int k=0;k<mesh.Nq;k++) {
            int id = n+k*mesh.Nq;
	    dfloat Dkm = mesh.D[m+k*mesh.Nq];
	    dfloat sx = mesh.vgeo[vbase + id + mesh.SXID*mesh.Np];
	    dfloat sy = mesh.vgeo[vbase + id + mesh.SYID*mesh.Np];
	    dfloat wJ = mesh.vgeo[vbase + id + mesh.JWID*mesh.Np];	    
	    dfloat nut_nk = nut[e*mesh.Np+id];
	    
	    dfloat uGss = (fac*sx*sx + sy*sy)*nut_nk*wJ;
            A[uid] += uGss*Dkm*Dkm; // strided for gather
	    
	    dfloat vGss = (fac*sy*sy + sx*sx)*nut_nk*wJ;
	    A[vid] += vGss*Dkm*Dkm;
          }

	  {
            int id = n+m*mesh.Nq;
	    dfloat Dnn = mesh.D[n+n*mesh.Nq];
	    dfloat Dmm = mesh.D[m+m*mesh.Nq];
	    dfloat rx = mesh.vgeo[vbase + id + mesh.RXID*mesh.Np];
	    dfloat ry = mesh.vgeo[vbase + id + mesh.RYID*mesh.Np];
	    dfloat sx = mesh.vgeo[vbase + id + mesh.SXID*mesh.Np];
	    dfloat sy = mesh.vgeo[vbase + id + mesh.SYID*mesh.Np];
	    dfloat wJ = mesh.vgeo[vbase + id + mesh.JWID*mesh.Np];	    
	    dfloat nut_nm = nut[e*mesh.Np+id];
	    
	    dfloat uGrs = 2.*(fac*rx*sx + ry*sy)*nut_nm*wJ;
            A[uid] += uGrs*Dnn*Dmm; // strided for gather
	    
	    dfloat vGrs = 2.*(fac*ry*sy + rx*sx)*nut_nm*wJ;
	    A[vid] += vGrs*Dnn*Dmm;

	    // do not need off diagonal blocks
	  }
	  
	  dfloat JW = mesh.vgeo[vbase + n + m*mesh.Nq + mesh.JWID*mesh.Np];
          A[uid] += JW*lambda;
	  A[vid] += JW*lambda;

        } else {
          A[uid] = 1; //just put a 1 so A is invertable
	  A[vid] = 1;
        }
      }
    }

    //add the rank boost for the allNeumann Poisson problem
    if (allNeumann) {
      for(int n=0;n<mesh.Np;++n){
        if (mapB[n+e*mesh.Np]!=1) { //dont fill rows for masked nodes
	  dlong id = e*mesh.Np+n;
	  dfloat fac = allNeumannPenalty*allNeumannScale*allNeumannScale;
          A[Nfields*id+0] += fac;
	  A[Nfields*id+1] += fac;
        }
      }
    }
  }
}


void stress_t::BuildOperatorDiagonalContinuousTet3D(memory<dfloat>& A) {

  std::cout << "BuildOperatorDiagonalContinuousTet3D not implemented" << std::endl;
  exit(-1);
}

void stress_t::BuildOperatorDiagonalContinuousHex3D(memory<dfloat>& A) {

  std::cout << "BuildOperatorDiagonalContinuousHex3D not implemented" << std::endl;
  exit(-1);
}


void stress_t::Setup(platform_t& _platform, mesh_t& _mesh,
		     settings_t& _settings, dfloat viscosity, dfloat _lambda,
                       const int _NBCTypes, const memory<int> _BCType){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;
  lambda = _lambda;

  Nfields = mesh.dim;

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
  std::string boundaryHeaderFileName;
  if (mesh.dim==2)
    boundaryHeaderFileName = std::string(DINS "/data/stressBoundary2D.h");
  else if (mesh.dim==3)
    boundaryHeaderFileName = std::string(DINS "/data/stressBoundary3D.h");
  kernelInfo["includes"] += boundaryHeaderFileName;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  kernelInfo["defines/" "p_Nfields"]= Nfields;
  
  int NblockV = std::max(1,blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  properties_t kernelInfoDouble = kernelInfo;
  kernelInfoDouble["defines/dfloat"] = "double";
  kernelInfoDouble["defines/dfloat4"] = "double4";

  properties_t kernelInfoFloat = kernelInfo;
  kernelInfoFloat["defines/dfloat"] = "float";
  kernelInfoFloat["defines/dfloat4"] = "float4";

  // Ax kernel (assume C0)
  fileName   = oklFilePrefix + "stressKernels" + suffix + oklFileSuffix;

  kernelName = "stressPartialAx" + suffix;
  partialAxKernel = platform.buildKernel(fileName, kernelName,
					 kernelInfoDouble);
  floatPartialAxKernel = platform.buildKernel(fileName, kernelName,
					      kernelInfoFloat);
  
  /* Preconditioner Setup */
  Ndofs = ogsMasked.Ngather*Nfields;
  Nhalo = gHalo.Nhalo*Nfields;

  nut.malloc(mesh.Np*mesh.Nelements);
  for(int e=0;e<mesh.Nelements;++e){
    for(int n=0;n<mesh.Np;++n){
      dlong id = e*mesh.Np + n;
      dfloat xn = mesh.x[id];
      dfloat yn = mesh.y[id];
      nut[id] = 1 + 0.5*(1+tanh(40.*(xn-1.2)))*(0.01/viscosity);  // 1 + excess scaled by 1/viscosity
    }
  }
  o_nut  = platform.malloc<dfloat>(mesh.Np*mesh.Nelements, nut);

  kernelName = "stressBuildOperatorDiagonal" + suffix;
  buildOperatorDiagonalKernel = platform.buildKernel(fileName, kernelName,
						     kernelInfo);

  // diagonal inverse (dfloat=>(pfloat)(1/float))
  kernelName = "stressReciprocal";
  reciprocalKernel = platform.buildKernel(fileName, kernelName,
					  kernelInfo);
  
  // assume Jacobi
  precon.Setup<StressJacobiPrecon>(*this);
  
}

void stress_t::ZeroMean(deviceMemory<double> &o_q){

  double qmean = platform.linAlg().sum(Ndofs, o_q, mesh.comm);

  // normalize
  qmean *= allNeumannScale*allNeumannScale;
  // q[n] = q[n] - qmean
  platform.linAlg().add(Ndofs, -qmean, o_q);
}

void stress_t::ZeroMean(deviceMemory<float> &o_q){

  float qmean = platform.linAlg().sum(Ndofs, o_q, mesh.comm);

  // normalize
  qmean *= allNeumannScale*allNeumannScale;
  // q[n] = q[n] - qmean
  platform.linAlg().add(Ndofs, -qmean, o_q);
}



void StressJacobiPrecon::Update(){
  stress.BuildOperatorDiagonal(o_invDiagA);
}

// Jacobi preconditioner
StressJacobiPrecon::StressJacobiPrecon(stress_t& _stress):
  stress(_stress) {

  o_invDiagA = stress.platform.malloc<pfloat>(stress.Ndofs);
  
  stress.BuildOperatorDiagonal(o_invDiagA);
  
}


void StressJacobiPrecon::Operator(deviceMemory<pfloat>& o_r, deviceMemory<pfloat>& o_Mr) {

  linAlg_t& linAlg = stress.platform.linAlg();

  // Mr = invDiag.*r
  linAlg.amxpy(stress.Ndofs, (pfloat)1.0, o_invDiagA, o_r, (pfloat)0.0, o_Mr);

  // zero mean of RHS
  if(stress.allNeumann){
    stress.ZeroMean(o_Mr);
  }
}
