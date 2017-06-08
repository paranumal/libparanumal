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
  // out  = REPORT, REPORT+VTU
  // adv  = CUBATURE, COLLOCATION
  // disc = DISCONT_GALERKIN, CONT_GALERKIN  

  char *options = strdup("out=REPORT+VTU, adv=CUBATURE, disc = DISCONT_GALERKIN");
  
  char *velSolverOptions = 
   strdup("type=VELOCITY solver=PCG,FLEXIBLE method=IPDG preconditioner=BLOCKJACOBI"); 
  char *prSolverOptions = 
   strdup("type=PRESSURE solver=PCG,FLEXIBLE method=IPDG preconditioner=FULLALMOND");

  if(argc!=3){
    printf("usage: ./main meshes/cavityH005.msh N\n");
    exit(-1);
  }
  // int specify polynomial degree 
  int N = atoi(argv[2]);
  // set up mesh stuff
  mesh2D *mesh = meshSetupQuad2D(argv[1], N);  
  //
  printf("Setup INS Solver: \n");   
  ins_t *ins = insSetupQuad2D(mesh,options,velSolverOptions,prSolverOptions); 

  // printf("OCCA Run: \n");  
  //insRun2D(ins,options); 
  

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
