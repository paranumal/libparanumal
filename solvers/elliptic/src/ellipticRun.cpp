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

#include "elliptic.hpp"
#include "timer.hpp"

void elliptic_t::Run(){

  //setup linear algebra module
  platform.linAlg().InitKernels({"set"});

  //setup linear solver
  hlong NglobalDofs;
  if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    NglobalDofs = ogsMasked.NgatherGlobal*Nfields;
  } else {
    NglobalDofs = mesh.NelementsGlobal*mesh.Np*Nfields;
  }

  linearSolver_t<dfloat> linearSolver;
  if (settings.compareSetting("LINEAR SOLVER","NBPCG")){
    linearSolver.Setup<LinearSolver::nbpcg<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","NBFPCG")){
    linearSolver.Setup<LinearSolver::nbfpcg<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","PCG")){
    linearSolver.Setup<LinearSolver::pcg<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","PGMRES")){
    linearSolver.Setup<LinearSolver::pgmres<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","PMINRES")){
    linearSolver.Setup<LinearSolver::pminres<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  }

  properties_t kernelInfo = mesh.props; //copy base occa properties

  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  //add standard boundary functions
  std::string boundaryHeaderFileName;
  if (mesh.dim==1)
    boundaryHeaderFileName = std::string(DELLIPTIC "/data/ellipticBoundary1D.h");
  else if (mesh.dim==2)
    boundaryHeaderFileName = std::string(DELLIPTIC "/data/ellipticBoundary2D.h");
  else if (mesh.dim==3)
    boundaryHeaderFileName = std::string(DELLIPTIC "/data/ellipticBoundary3D.h");
  kernelInfo["includes"] += boundaryHeaderFileName;

  int Nmax = std::max(mesh.Np, mesh.Nfaces*mesh.Nfp);
  kernelInfo["defines/" "p_Nmax"]= Nmax;

  kernelInfo["defines/" "p_Nfields"]= Nfields;

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();
  
  std::string oklFilePrefix = DELLIPTIC "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  fileName   = oklFilePrefix + "ellipticRhs" + suffix + oklFileSuffix;
  kernelName = "ellipticRhs" + suffix;
  kernel_t forcingKernel = platform.buildKernel(fileName, kernelName,
                                                    kernelInfo);

  kernel_t rhsBCKernel, addBCKernel;
  if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    fileName   = oklFilePrefix + "ellipticRhsBCIpdg" + suffix + oklFileSuffix;
    kernelName = "ellipticRhsBCIpdg" + suffix;

    rhsBCKernel = platform.buildKernel(fileName,kernelName, kernelInfo);
  } else if (settings.compareSetting("DISCRETIZATION","CONTINUOUS")) {
    fileName   = oklFilePrefix + "ellipticRhsBC" + suffix + oklFileSuffix;
    kernelName = "ellipticRhsBC" + suffix;

    rhsBCKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

    fileName   = oklFilePrefix + "ellipticAddBC" + suffix + oklFileSuffix;
    kernelName = "ellipticAddBC" + suffix;

    addBCKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
  }

  //create occa buffers
  dlong Nall = mesh.Np*(mesh.Nelements+mesh.totalHaloPairs);
  memory<dfloat> rL(Nall);
  memory<dfloat> xL(Nall);
  deviceMemory<dfloat> o_rL = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_xL = platform.malloc<dfloat>(Nall);

  deviceMemory<dfloat> o_r, o_x;
  if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    o_r = o_rL;
    o_x = o_xL;
  } else {
    dlong Ng = ogsMasked.Ngather;
    dlong Nghalo = gHalo.Nhalo;
    dlong Ngall = Ng + Nghalo;
    o_r = platform.malloc<dfloat>(Ngall);
    o_x = platform.malloc<dfloat>(Ngall);
  }

  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  //populate rhs forcing
  forcingKernel(mesh.Nelements,
                mesh.o_wJ,
                mesh.o_MM,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                lambda,
                o_rL);

  //Set x to zero
  platform.linAlg().set(mesh.Nelements*mesh.Np*Nfields, (dfloat)0.0, o_xL);

  //add boundary condition contribution to rhs
  if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    rhsBCKernel(mesh.Nelements,
                mesh.o_vmapM,
                tau,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                mesh.o_vgeo,
                mesh.o_sgeo,
                o_EToB,
                mesh.o_D,
                mesh.o_LIFT,
                mesh.o_MM,
                o_rL);
  } else if (settings.compareSetting("DISCRETIZATION","CONTINUOUS")) {
    rhsBCKernel(mesh.Nelements,
                mesh.o_wJ,
                mesh.o_ggeo,
                mesh.o_sgeo,
                mesh.o_D,
                mesh.o_S,
                mesh.o_MM,
                mesh.o_vmapM,
                mesh.o_sM,
                lambda,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                o_mapB,
                o_rL);
  }

  // gather rhs to globalDofs if c0
  if(settings.compareSetting("DISCRETIZATION","CONTINUOUS")){
    ogsMasked.Gather(o_r, o_rL, 1, ogs::Add, ogs::Trans);
    ogsMasked.Gather(o_x, o_xL, 1, ogs::Add, ogs::NoTrans);
  }

  int maxIter = 5000;
  int verbose = settings.compareSetting("VERBOSE", "TRUE") ? 1 : 0;

  timePoint_t start = GlobalPlatformTime(platform);

  //call the solver
  dfloat tol = (sizeof(dfloat)==sizeof(double)) ? 1.0e-18 : 1.0e-5;
  int iter = Solve(linearSolver, o_x, o_r, tol, maxIter, verbose);

  //add the boundary data to the masked nodes
  if(settings.compareSetting("DISCRETIZATION","CONTINUOUS")){
    // scatter x to LocalDofs if c0
    ogsMasked.Scatter(o_xL, o_x, 1, ogs::NoTrans);
    //fill masked nodes with BC data
    addBCKernel(mesh.Nelements,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                o_mapB,
                o_xL);
  }

  timePoint_t end = GlobalPlatformTime(platform);
  double elapsedTime = ElapsedTime(start, end);

  if ((mesh.rank==0) && verbose){
    printf("%d, " hlongFormat ", %g, %d, %g, %g; global: N, dofs, elapsed, iterations, time per node, nodes*iterations/time %s\n",
           mesh.N,
           NglobalDofs,
           elapsedTime,
           iter,
           elapsedTime/(NglobalDofs),
           NglobalDofs*((dfloat)iter/elapsedTime),
           (char*) settings.getSetting("PRECONDITIONER").c_str());
  }


#if 1

  dfloat normr = platform.linAlg().norm2(o_r.length(), o_r, comm);
  
  // make sure cubature is set up
  mesh.CubatureSetup();
  mesh.CubaturePhysicalNodes();
  
  occa::properties kernelInfoCubNp = kernelInfo;
  
  dlong p_NblockC = 1; 
  dlong NblocksC = (mesh.Nelements+p_NblockC-1)/p_NblockC;
  deviceMemory<dfloat> o_errH1 = platform.malloc<dfloat>(NblocksC);
  
  kernelInfoCubNp["defines/" "p_NblockC"] = p_NblockC;
  kernelInfoCubNp["defines/" "p_cubNp"]   = mesh.cubNp;
  kernelInfoCubNp["defines/" "p_cubNq"]   = mesh.cubNq;

  fileName   = oklFilePrefix + "/ellipticCubatureH1Error" + suffix + oklFileSuffix;
  kernelName = "ellipticCubatureH1Error" + suffix;
  kernel_t cubatureH1ErrorKernel = platform.buildKernel(fileName, kernelName, kernelInfoCubNp);
  
  NblocksC = (mesh.Nelements+p_NblockC-1)/p_NblockC;
  platform.linAlg().set(NblocksC, (dfloat)0.0, o_errH1);
  
  cubatureH1ErrorKernel(mesh.Nelements,
			mesh.o_cubx, mesh.o_cuby, mesh.o_cubz, mesh.o_cubwJ,
			mesh.o_cubvgeo, mesh.o_cubD, mesh.o_cubInterp, lambda,
			o_xL, o_errH1);
  
  dfloat normH1 = platform.linAlg().sum(NblocksC, o_errH1, mesh.comm);
  normH1 = sqrt(normH1);
  
  printf("%d &  %d & %d &  %.2e &  %.2e & %.2e; %% CG (REPORT): N, Ndofs, iter, normH1, ||r||, solveTime\n",
	 mesh.N, (int)  NglobalDofs,  iter, normH1, normr, elapsedTime);

#if 1
  {
    o_xL.copyTo(xL);
    FILE *fp = fopen("foo.dat", "w");
    fprintf(fp,"%%%% element, node, x, soln, exactSoln, error, Jacobian, drdx\n");
    for(int e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.Np;++n){
	int id = e*mesh.Np+n;
	dfloat xn = mesh.x[id];
	dfloat exu = sin(M_PI*xn)*exp(-10*sin(M_PI*xn)*sin(M_PI*xn));
	fprintf(fp,"%d, %d, % .10e, % .10e, % .10e, % .10e, % .10e, % .10e\n",
		e,
		n,
		mesh.x[id],
		xL[id],
		exu,
		xL[id]-exu,
		mesh.vgeo[mesh.Nvgeo*mesh.Np*e+n+mesh.Np*mesh.JID],
		mesh.vgeo[mesh.Nvgeo*mesh.Np*e+n+mesh.Np*mesh.RXID]);
      }
    }
    fclose(fp);
  }
#endif	     
      
  
#endif


  
  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

    // copy data back to host
    o_xL.copyTo(xL);

    // output field files
    std::string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d.vtu", name.c_str(), mesh.rank);

    PlotFields(xL, fname);
  }

  // output norm of final solution
  {
    //compute q.M*q
    dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
    deviceMemory<dfloat> o_MxL = platform.reserve<dfloat>(Nentries);
    mesh.MassMatrixApply(o_xL, o_MxL);

    dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_xL, o_MxL, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);
  }
}
