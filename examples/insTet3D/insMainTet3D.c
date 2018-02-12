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

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  char *velSolverOptions = 
    strdup("solver=PCG method=IPDG preconditioner=MULTIGRID smoother=CHEBYSHEV");
  char *velParAlmondOptions = 
    strdup("solver=KCYCLE smoother=CHEBYSHEV partition=STRONGNODES");

  char *prSolverOptions =
    //strdup("solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID, HALFDOFS smoother=DAMPEDJACOBI,CHEBYSHEV");
    strdup("solver=PCG,FLEXIBLE method=CONTINUOUS preconditioner=FULLALMOND");

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
  mesh3D *mesh = meshSetupTet3D(argv[1], N); 

  // capture header file
  char *boundaryHeaderFileName;
  if(argc==3)
    boundaryHeaderFileName = strdup(DHOLMES "/examples/insTet3D/insUniform3D.h"); // default
  else
    boundaryHeaderFileName = strdup(argv[3]);

  //int Ns = 0; // Default no-subcycling 
  int Ns =4;
  if(argc==5)
   Ns = atoi(argv[4]); // Number of substeps
  
  
  char *options; 
 if(Ns==0)
      options = strdup("method = ALGEBRAIC, grad-div= BROKEN, out=ADAPTIVECONTOUR, adv=CUBATURE, disc = DISCONT_GALERKIN"); // SUBCYCLING
  else
      options = strdup("method = ALGEBRAIC, grad-div= BROKEN, SUBCYCLING, out=ADAPTIVECONTOUR, adv=CUBATURE, disc = DISCONT_GALERKIN"); // SUBCYCLING

  if (rank==0) printf("Setup INS Solver: \n");
  ins_t *ins = insSetup3D(mesh, Ns, options,
                          velSolverOptions,velParAlmondOptions,
                          prSolverOptions, prParAlmondOptions,
                          boundaryHeaderFileName);

  if (rank==0) printf("OCCA Run: \n");
  insRun3D(ins,options);


   // printf("OCCA Run Timer: \n");
   // insRunTimer3D(mesh,options,boundaryHeaderFileName);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
