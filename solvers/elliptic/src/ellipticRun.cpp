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

  //setup linear solver
  hlong NglobalDofs;
  if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    NglobalDofs = ogsMasked.NgatherGlobal*Nfields;
  } else {
    NglobalDofs = mesh.NelementsGlobal*mesh.Np*Nfields;
  }

  linearSolver_t linearSolver;
  if (settings.compareSetting("LINEAR SOLVER","NBPCG")){
    linearSolver.Setup<LinearSolver::nbpcg>(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","NBFPCG")){
    linearSolver.Setup<LinearSolver::nbfpcg>(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","PCG")){
    linearSolver.Setup<LinearSolver::pcg>(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","PGMRES")){
    linearSolver.Setup<LinearSolver::pgmres>(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","PMINRES")){
    linearSolver.Setup<LinearSolver::pminres>(Ndofs, Nhalo, platform, settings, comm);
  }

  properties_t kernelInfo = mesh.props; //copy base occa properties

  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  //add standard boundary functions
  std::string boundaryHeaderFileName;
  if (mesh.dim==2)
    boundaryHeaderFileName = std::string(DELLIPTIC "/data/ellipticBoundary2D.h");
  else if (mesh.dim==3)
    boundaryHeaderFileName = std::string(DELLIPTIC "/data/ellipticBoundary3D.h");
  kernelInfo["includes"] += boundaryHeaderFileName;

  int Nmax = std::max(mesh.Np, mesh.Nfaces*mesh.Nfp);
  kernelInfo["defines/" "p_Nmax"]= Nmax;

  kernelInfo["defines/" "p_Nfields"]= Nfields;

  // set kernel name suffix
  std::string suffix;
  if(mesh.elementType==Mesh::TRIANGLES) {
    suffix = "Tri2D";
  } else if(mesh.elementType==Mesh::QUADRILATERALS) {
    if(mesh.dim==2)
      suffix = "Quad2D";
    else
      suffix = "Quad3D";
  } else if(mesh.elementType==Mesh::TETRAHEDRA) {
    suffix = "Tet3D";
  } else { //mesh.elementType==Mesh::HEXAHEDRA)
    suffix = "Hex3D";
  }

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
  memory<dfloat> rL(Nall, 0.0);
  memory<dfloat> xL(Nall, 0.0);
  deviceMemory<dfloat> o_rL = platform.malloc<dfloat>(rL);
  deviceMemory<dfloat> o_xL = platform.malloc<dfloat>(xL);

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

  //storage for M*q during reporting
  deviceMemory<dfloat> o_MxL = platform.malloc<dfloat>(xL);
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
  dfloat tol = (sizeof(dfloat)==sizeof(double)) ? 1.0e-8 : 1.0e-5;
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
    mesh.MassMatrixApply(o_xL, o_MxL);

    dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
    dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_xL, o_MxL, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);
  }
}
