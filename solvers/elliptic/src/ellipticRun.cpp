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
  linearSolver_t *linearSolver = linearSolver_t::Setup(*this);

  int weighted = settings.compareSetting("DISCRETIZATION", "CONTINUOUS") ? 1 : 0;
  linearSolver->Init(weighted, o_weight);

  occa::properties kernelInfo = props; //copy base occa properties

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
  sprintf(fileName, DELLIPTIC "/okl/ellipticRhs%s.okl", suffix);
  sprintf(kernelName, "ellipticRhs%s", suffix);

  occa::kernel forcingKernel = buildKernel(device, fileName, kernelName,
                                                    kernelInfo, comm);

  occa::kernel rhsBCKernel, addBCKernel;
  if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    sprintf(fileName, DELLIPTIC "/okl/ellipticRhsBCIpdg%s.okl", suffix);
    sprintf(kernelName, "ellipticRhsBCIpdg%s", suffix);

    rhsBCKernel = buildKernel(device, fileName,kernelName, kernelInfo, comm);
  } else if (settings.compareSetting("DISCRETIZATION","CONTINUOUS")) {
    sprintf(fileName, DELLIPTIC "/okl/ellipticRhsBC%s.okl", suffix);
    sprintf(kernelName, "ellipticRhsBC%s", suffix);

    rhsBCKernel = buildKernel(device, fileName, kernelName, kernelInfo, comm);

    sprintf(fileName, DELLIPTIC "/okl/ellipticAddBC%s.okl", suffix);
    sprintf(kernelName, "ellipticAddBC%s", suffix);

    addBCKernel = buildKernel(device, fileName, kernelName, kernelInfo, comm);
  }

  //create occa buffers
  dlong Nall = mesh.Np*(mesh.Nelements+mesh.totalHaloPairs);
  dfloat *r = (dfloat*) calloc(Nall, sizeof(dfloat));
  dfloat *x = (dfloat*) calloc(Nall, sizeof(dfloat));
  occa::memory o_r = device.malloc(Nall*sizeof(dfloat));
  occa::memory o_x = device.malloc(Nall*sizeof(dfloat), x);

  //populate rhs forcing
  forcingKernel(mesh.Nelements,
                mesh.o_ggeo,
                mesh.o_MM,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                lambda,
                o_r);

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
                mesh.o_Dmatrices,
                mesh.o_LIFTT,
                mesh.o_MM,
                o_r);
  } else if (settings.compareSetting("DISCRETIZATION","CONTINUOUS")) {
    rhsBCKernel(mesh.Nelements,
                mesh.o_ggeo,
                mesh.o_sgeo,
                mesh.o_Dmatrices,
                mesh.o_Smatrices,
                mesh.o_MM,
                mesh.o_vmapM,
                mesh.o_sMT,
                lambda,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                o_mapB,
                o_r);
  }

  // gather-scatter and mask rhs if c0
  if(settings.compareSetting("DISCRETIZATION","CONTINUOUS")){
    mesh.ogs->GatherScatter(o_r, ogs_dfloat, ogs_add, ogs_sym);
    if (Nmasked) maskKernel(Nmasked, o_maskIds, o_r);
  }

  int maxIter = 5000;
  int verbose = settings.compareSetting("VERBOSE", "TRUE") ? 1 : 0;

  MPI_Barrier(comm);
  double startTime = MPI_Wtime();

  //call the solver
  dfloat tol = 1e-10;
  int iter = Solve(*linearSolver, o_x, o_r, tol, maxIter, verbose);

  //add the boundary data to the masked nodes
  if(settings.compareSetting("DISCRETIZATION","CONTINUOUS")){
    addBCKernel(mesh.Nelements,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                o_mapB,
                o_x);
  }

  MPI_Barrier(comm);
  double endTime = MPI_Wtime();
  double elapsedTime = endTime - startTime;

#if 1
  o_x.copyTo(x);
  double err = 0.f; 
  for(int e=0; e<mesh.Nelements; e++){
    for(int n=0; n<mesh.Np; n++){
      const int id       = n + e*mesh.Np; 
       // const double exact = sin(M_PI*mesh.x[id])*sin(M_PI*mesh.y[id]);
      const double exact = sin(M_PI*mesh.x[id])*sin(M_PI*mesh.y[id])*sin(M_PI*mesh.z[id]); ;
      err = mymax(fabs(x[id]-exact), err);
    }
  }
  printf("Error = %.8e \n", err); 
#endif

  if ((mesh.rank==0) && verbose){
    printf("%d, " hlongFormat ", %g, %d, %g, %g; global: N, dofs, elapsed, iterations, time per node, nodes*iterations/time %s\n",
           mesh.N,
           mesh.NelementsGlobal*mesh.Np,
           elapsedTime,
           iter,
           elapsedTime/(mesh.Np*mesh.NelementsGlobal),
           mesh.NelementsGlobal*((dfloat)iter*mesh.Np/elapsedTime),
           (char*) settings.getSetting("PRECONDITIONER").c_str());
  }

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

    // copy data back to host
    o_x.copyTo(x);


#if 1
  for(int e=0; e<mesh.Nelements; e++){
    for(int n=0; n<mesh.Np; n++){
      const int id       = n + e*mesh.Np; 
       // const double exact = sin(M_PI*mesh.x[id])*sin(M_PI*mesh.y[id]);
     const double exact = sin(M_PI*mesh.x[id])*sin(M_PI*mesh.y[id])*sin(M_PI*mesh.z[id]); ;
       x[id] = x[id]-exact;
    }
  }
#endif

    // output field files
    string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d.vtu", name.c_str(), mesh.rank);

    PlotFields(x, fname);
  }

  free(r); free(x);
  o_r.free(); o_x.free();
}
