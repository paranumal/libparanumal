#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"
#include "ins2D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  // SET OPTIONS
  // method = ALGEBRAIC, STIFFLYSTABLE (default for now)
  // grad-div   = WEAK, NODAL (default nodal)
  // out  = REPORT, REPORT+VTU
  // adv  = CUBATURE, COLLOCATION
  // disc = DISCONT_GALERKIN, CONT_GALERKIN 
  //char *options = strdup("method = ALGEBRAIC, grad-div= BROKEN, out=REPORT, adv=CUBATURE, disc = DISCONT_GALERKIN"); // SUBCYCLING
  


  char *velSolverOptions =
    strdup("solver=PCG method=IPDG basis=NODAL preconditioner=MASSMATRIX");
  char *velParAlmondOptions =
    strdup("solver= smoother= partition=");

  char *prSolverOptions =
    //strdup("solver=PCG,FLEXIBLE method=IPDG basis=NODAL preconditioner=MULTIGRID,HALFDOFS smoother=DAMPEDJACOBI,CHEBYSHEV");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG basis=NODAL preconditioner=NONE");
    strdup("solver=PCG,FLEXIBLE, method=IPDG basis=NODAL preconditioner=FULLALMOND");
    //strdup("solver=PCG,FLEXIBLE, method=IPDG basis=NODAL preconditioner=OMS,APPROXPATCH coarse=COARSEGRID,ALMOND");

  char *prParAlmondOptions =
    strdup("solver=KCYCLE smoother=CHEBYSHEV partition=STRONGNODES");

  if(argc!=3 && argc!=4 && argc!=5){
    printf("usage 1: ./main meshes/cavityH005.msh N\n");
    printf("usage 2: ./main meshes/cavityH005.msh N insUniformFlowBoundaryConditions.h\n");
    printf("usage 3: ./main meshes/cavityH005.msh N insUniformFlowBoundaryConditions.h Nsubstep\n");
    exit(-1);
  }
  // int specify polynomial degree
  int N = atoi(argv[2]);

  // set up mesh stuff
  mesh2D *mesh = meshSetupTri2D(argv[1], N); 

  #if 1

  // capture header file
  char *boundaryHeaderFileName;
  if(argc==3)
    boundaryHeaderFileName = strdup(DHOLMES "/examples/insTri2D/insUniform2D.h"); // default
  else
    boundaryHeaderFileName = strdup(argv[3]);

  int Ns = 0; // Default no-subcycling 
  if(argc==5)
   Ns = atoi(argv[4]); // Number of substeps
  
  
  char *options; 
  if(Ns==0)
      options = strdup("method = ALGEBRAIC, grad-div= BROKEN, out=VTU, adv=CUBATURE, disc = DISCONT_GALERKIN"); // SUBCYCLING
  else
      options = strdup("method = ALGEBRAIC, grad-div= BROKEN, SUBCYCLING, out=VTU, adv=CUBATURE, disc = DISCONT_GALERKIN"); // SUBCYCLING

    #if 1
     printf("Setup INS Solver: \n");
    ins_t *ins = insSetup2D(mesh,Ns,options, velSolverOptions,   velParAlmondOptions,
                            prSolverOptions, prParAlmondOptions, boundaryHeaderFileName);
    printf("Running INS solver\n");
    insRun2D(ins,options);
    #else
    printf("OCCA Run Timer: \n");
    insRunTimer2D(mesh,options,boundaryHeaderFileName);
    #endif
    
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
  #endif
}
