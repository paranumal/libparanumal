#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"
#include "boltzmann2D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  
  // SET OPTIONS
  // mode        = TEST, SOLVER // do not use test mode, for developing purposes
  // relaxation  = CUBATURE, COLLOCATION, 
  // time        = LSERK, LSIMEX, SARK3, SAAB3
  // out         = REPORT, REPORT-VTU, NO  
  // bc          = UNSPLITPML, SPLITPML, NONE
  // pmlprofile  = CONSTANT, QUADRATIC
  
  char *options =strdup("out = REPORT-VTU,relaxation = CUBATURE,time = LSIMEX, bc=PML, pmlprofile=CONSTANT");
  
    if(argc!=3){
          printf("usage: ./main meshes/cavityH005.msh N\n");
          exit(-1);
    }
    // int specify polynomial degree 
    int N = atoi(argv[2]);
    // set up mesh stuff
    mesh2D *mesh = meshSetupTri2D(argv[1], N);  

    printf("Setup Boltzmann Solver: \n");   
    boltzmannSetup2D(mesh,options);  
    
    printf("Occa Run: \n");  
    boltzmannRun2D(mesh,options); 
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
