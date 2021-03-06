/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

void elliptic_t::Run(){

  //setup linear solver
  hlong NglobalDofs;
  if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    NglobalDofs = ogsMasked->NgatherGlobal*Nfields;
  } else {
    NglobalDofs = mesh.NelementsGlobal*mesh.Np*Nfields;
  }
  linearSolver_t *linearSolver = linearSolver_t::Setup(Ndofs, Nhalo,
                                                       platform, settings, mesh.comm);

  occa::properties kernelInfo = mesh.props; //copy base occa properties

  string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  //add standard boundary functions
  char *boundaryHeaderFileName;
  if (mesh.dim==2)
    boundaryHeaderFileName = strdup(DELLIPTIC "/data/ellipticBoundary2D.h");
  else if (mesh.dim==3)
    boundaryHeaderFileName = strdup(DELLIPTIC "/data/ellipticBoundary3D.h");
  kernelInfo["includes"] += boundaryHeaderFileName;

  int Nmax = mymax(mesh.Np, mesh.Nfaces*mesh.Nfp);
  kernelInfo["defines/" "p_Nmax"]= Nmax;

  kernelInfo["defines/" "p_Nfields"]= Nfields;

  // set kernel name suffix
  char *suffix;
  if(mesh.elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(mesh.elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");

  if(mesh.elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(mesh.elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  if(mesh.elementType==QUADRILATERALS){
    if(mesh.dim==2)
      suffix = strdup("Quad2D");
    else
      suffix = strdup("Quad3D");
  }

  sprintf(fileName, DELLIPTIC "/okl/ellipticRhs%s.okl", suffix);
  sprintf(kernelName, "ellipticRhs%s", suffix);
  occa::kernel forcingKernel = platform.buildKernel(fileName, kernelName,
                                                    kernelInfo);

  occa::kernel rhsBCKernel, addBCKernel;
  if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    sprintf(fileName, DELLIPTIC "/okl/ellipticRhsBCIpdg%s.okl", suffix);
    sprintf(kernelName, "ellipticRhsBCIpdg%s", suffix);

    rhsBCKernel = platform.buildKernel(fileName,kernelName, kernelInfo);
  } else if (settings.compareSetting("DISCRETIZATION","CONTINUOUS")) {
    sprintf(fileName, DELLIPTIC "/okl/ellipticRhsBC%s.okl", suffix);
    sprintf(kernelName, "ellipticRhsBC%s", suffix);

    rhsBCKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

    sprintf(fileName, DELLIPTIC "/okl/ellipticAddBC%s.okl", suffix);
    sprintf(kernelName, "ellipticAddBC%s", suffix);

    addBCKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
  }

  //create occa buffers
  dlong Nall = mesh.Np*(mesh.Nelements+mesh.totalHaloPairs);
  dfloat *rL = (dfloat*) calloc(Nall, sizeof(dfloat));
  dfloat *xL = (dfloat*) calloc(Nall, sizeof(dfloat));
  occa::memory o_rL = platform.malloc(Nall*sizeof(dfloat), rL);
  occa::memory o_xL = platform.malloc(Nall*sizeof(dfloat), xL);

  occa::memory o_r, o_x;
  if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    o_r = o_rL;
    o_x = o_xL;
  } else {
    dlong Ng = ogsMasked->Ngather;
    dlong Nghalo = ogsMasked->NgatherHalo;
    dlong Ngall = Ng + Nghalo;
    o_r = platform.malloc(Ngall*sizeof(dfloat));
    o_x = platform.malloc(Ngall*sizeof(dfloat));
  }

  //storage for M*q during reporting
  occa::memory o_MxL = platform.malloc(Nall*sizeof(dfloat), xL);
  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  //populate rhs forcing
  forcingKernel(mesh.Nelements,
                mesh.o_ggeo,
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
    ogsMasked->Gather(o_r, o_rL, ogs_dfloat, ogs_add, ogs_trans);
    ogsMasked->Gather(o_x, o_xL, ogs_dfloat, ogs_add, ogs_notrans);
  }

  int maxIter = 5000;
  int verbose = settings.compareSetting("VERBOSE", "TRUE") ? 1 : 0;

  MPI_Barrier(mesh.comm);
  double startTime = MPI_Wtime();

  //call the solver
  dfloat tol = 1e-8;
  int iter = Solve(*linearSolver, o_x, o_r, tol, maxIter, verbose);


  //add the boundary data to the masked nodes
  if(settings.compareSetting("DISCRETIZATION","CONTINUOUS")){
    // scatter x to LocalDofs if c0
    ogsMasked->Scatter(o_xL, o_x, ogs_dfloat, ogs_add, ogs_notrans);
    //fill masked nodes with BC data
    addBCKernel(mesh.Nelements,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                o_mapB,
                o_xL);
  }

  MPI_Barrier(mesh.comm);
  double endTime = MPI_Wtime();
  double elapsedTime = endTime - startTime;

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
    string name;
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
    dfloat norm2 = sqrt(linAlg.innerProd(Nentries, o_xL, o_MxL, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);
  }

  free(rL); free(xL);
  o_rL.free(); o_xL.free();
  o_r.free(); o_x.free();
  o_MxL.free();
  delete linearSolver;
}
