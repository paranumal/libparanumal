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
  char *options = strdup("method = ALGEBRAIC, grad-div= BROKEN, out=REPORT+VTU, adv=CUBATURE, disc = DISCONT_GALERKIN"); // SUBCYCLING
  //  char *options = strdup("out=REPORT+VTU, adv=COLLOCATION, disc = DISCONT_GALERKIN");
 
  char *velSolverOptions =
    strdup("solver=PCG method=IPDG preconditioner=MASSMATRIX");
  char *velParAlmondOptions =
    strdup("solver= smoother= partition=");

  char *prSolverOptions =
   strdup("solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=DAMPEDJACOBI,CHEBYSHEV");
    //strdup("solver=PCG,FLEXIBLE,method=IPDG  preconditioner=FULLALMOND");
    //strdup("solver=PCG,FLEXIBLE, method=IPDG preconditioner=OMS,APPROXPATCH coarse=COARSEGRID,ALMOND");

  char *prParAlmondOptions =
    strdup("solver=KCYCLE smoother=CHEBYSHEV partition=STRONGNODES");

  if(argc!=3 && argc!=4){
    printf("usage 1: ./main meshes/cavityH005.msh N\n");
    printf("usage 2: ./main meshes/cavityH005.msh N insUniformFlowBoundaryConditions.h\n");
    exit(-1);
  }
  // int specify polynomial degree
  int N = atoi(argv[2]);

  // set up mesh stuff
  mesh2D *mesh = meshSetupTri2D(argv[1], N); 

  // capture header file
  char *boundaryHeaderFileName;
  if(argc==3)
    boundaryHeaderFileName = strdup(DHOLMES "/examples/insTri2D/insUniform2D.h"); // default
  else
    boundaryHeaderFileName = strdup(argv[3]);

  for(iint i=0; i<1; i++){
    //iint i=0; 
    printf("Setup INS Solver: \n");
    ins_t *ins = insSetup2D(mesh,i,options, velSolverOptions,velParAlmondOptions,
                            prSolverOptions, prParAlmondOptions, boundaryHeaderFileName);

    printf("OCCA Run: \n");
    insRun2D(ins,options);
  }

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
