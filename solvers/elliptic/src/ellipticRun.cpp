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
#include "mesh/meshDefines2D.h"
#include "mesh/meshDefines3D.h"

void elliptic_t::Run(){

  //setup linear solver
  int weighted = settings.compareSetting("DISCRETIZATION", "CONTINUOUS") ? 1 : 0;

  dlong Nlocal = mesh.Nelements*mesh.Np*Nfields;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np*Nfields;
  linearSolver_t *linearSolver = linearSolver_t::Setup(Nlocal, Nhalo,
                                                       platform, settings, mesh.comm,
                                                       weighted, o_weight);

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
  dfloat *r = (dfloat*) calloc(Nall, sizeof(dfloat));
  dfloat *x = (dfloat*) calloc(Nall, sizeof(dfloat));
  occa::memory o_r = platform.malloc(Nall*sizeof(dfloat), r);
  occa::memory o_x = platform.malloc(Nall*sizeof(dfloat), x);

  //storage for M*q during reporting
  occa::memory o_Mx = platform.malloc(Nall*sizeof(dfloat), x);
  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  //populate rhs forcing
  forcingKernel(mesh.Nelements,
                mesh.o_ggeo,
                mesh.o_MM,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                lambda,
                o_r);

#if 0
  dfloat *tmpr = (dfloat*) calloc(mesh.Np*mesh.Nelements, sizeof(dfloat));
  for(dlong e=0;e<mesh.Nelements;++e){
    for(dlong n=0;n<mesh.Np;++n){
      tmpr[e*mesh.Np+n] = drand48()*mesh.ggeo[e*mesh.Np*mesh.Nggeo + n + GWJID*mesh.Np];
    }
  }
  o_r.copyFrom(tmpr);
#endif
  
#if 0
  dfloat *h_r = (dfloat*) calloc(mesh.Np*mesh.Nelements, sizeof(dfloat));
  for(dlong e=0;e<mesh.Nelements;++e){
    for(dlong n=0;n<mesh.Np;++n){
      h_r[e*mesh.Np+n] =
	drand48()*mesh.ggeo[mesh.Np*mesh.Nggeo*e+n+GWJID*mesh.Np];
    }
  }
  o_r.copyFrom(h_r);
#endif
  
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
                o_r);
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
                o_r);
  }

  // gather-scatter and mask rhs if c0
  if(settings.compareSetting("DISCRETIZATION","CONTINUOUS")){
    mesh.ogs->GatherScatter(o_r, ogs_dfloat, ogs_add, ogs_sym);
    if (Nmasked) maskKernel(Nmasked, o_maskIds, o_r);
  }

  int maxIter = 5000;
  int verbose = settings.compareSetting("VERBOSE", "TRUE") ? 1 : 0;

  MPI_Barrier(mesh.comm);
  double startTime = MPI_Wtime();

  //call the solver
  dfloat tol = 1e-8;
  settings.getSetting("LINEAR SOLVER CONVERGENCE TOLERANCE", tol);
  
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

  MPI_Barrier(mesh.comm);
  double endTime = MPI_Wtime();
  double elapsedTime = endTime - startTime;

  if ((mesh.rank==0) && verbose){
    dfloat epsy = 1.;
    mesh.settings.getSetting("BOX COORDINATE MAP PARAMETER Y", epsy);
    
    printf("%d, " hlongFormat ", %g, %d, %g, %g, %g; global: N, dofs, elapsed, iterations, time per node, nodes*iterations/time, epsy, %s\n",
           mesh.N,
           mesh.NelementsGlobal*mesh.Np,
           elapsedTime,
           iter,
           elapsedTime/(mesh.Np*mesh.NelementsGlobal),
           mesh.NelementsGlobal*((dfloat)iter*mesh.Np/elapsedTime),
	   epsy,
           (char*) settings.getSetting("PRECONDITIONER").c_str());
  }

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

    // copy data back to host
    o_x.copyTo(x);

    // output field files
    string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d.vtu", name.c_str(), mesh.rank);

    PlotFields(x, fname);
  }

  // output norm of final solution
  {
    //compute q.M*q
    mesh.MassMatrixApply(o_x, o_Mx);

    dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
    dfloat norm2 = sqrt(linAlg.innerProd(Nentries, o_x, o_Mx, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);
  }

  free(r); free(x);
  o_r.free(); o_x.free();
  o_Mx.free();
  delete linearSolver;
}
