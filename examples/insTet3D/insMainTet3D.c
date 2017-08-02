#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh3D.h"
#include "ins3D.h"

int main(int argc, char **argv){
  // start up MPI
  MPI_Init(&argc, &argv);

  char *options = strdup("out=REPORT+VTU, adv=CUBATURE, disc = DISCONT_GALERKIN"); // SUBCYCLING
  //  char *options = strdup("out=REPORT+VTU, adv=COLLOCATION, disc = DISCONT_GALERKIN");
 
  char *velSolverOptions = 
      strdup("solver=PCG method=IPDG preconditioner=BLOCKJACOBI");

  char *prSolverOptions =
      strdup("solver=PCG,FLEXIBLE method=IPDG,preconditioner=FULLALMOND"); // ,FORCESYMMETRY"); // ,FORCESYMMETRY");

  if(argc!=3 && argc!=4){
    printf("usage 1: ./main meshes/cavityH005.msh N\n");
    printf("usage 2: ./main meshes/cavityH005.msh N insUniformFlowBoundaryConditions.h\n");
    exit(-1);
  }
  // int specify polynomial degree
  int N = atoi(argv[2]);

  // set up mesh stuff
  mesh3D *mesh = meshSetupTet3D(argv[1], N); 

  // capture header file
  char *boundaryHeaderFileName;
  if(argc==3)
    boundaryHeaderFileName = strdup(DHOLMES "/examples/insTet3D/insUniform3D.h"); // default
  else
    boundaryHeaderFileName = strdup(argv[3]);

  printf("Setup INS Solver: \n");
  ins_t *ins = insSetup3D(mesh,options,velSolverOptions,prSolverOptions,boundaryHeaderFileName);

  printf("OCCA Run: \n");
  insRun3D(ins,options);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
