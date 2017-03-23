#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=4){
    // to run cavity test case with degree N elements
    printf("usage: ./main meshes/cavityH005.msh N dtfactor\n");
    exit(-1);
  }

  // int specify polynomial degree 
  int N = atoi(argv[2]);

  int dtfactor = atoi(argv[3]);

  // set up mesh stuff
  mesh2D *mesh = meshSetupTri2D(argv[1], N);
  

  // set up boltzmann stuff
  void boltzmannSplitPmlSetup2D(mesh2D *mesh);
  void boltzmannSplitPmlRun2D(mesh2D *mesh);
  

  // run
  //  boltzmannRun2D(mesh);
  printf("occa run: \n");

  



  FILE *fp;
  char fname[BUFSIZ];
  sprintf(fname, "TimeStepSize.txt");
  

  mesh->dtfactor = dtfactor; 
  
  boltzmannSplitPmlSetup2D(mesh);   
  boltzmannSplitPmlRun2D(mesh);
   

  fp = fopen(fname, "a");
  fprintf(fp, "Time steping method, TauInv, dt , Error\n");
  fprintf(fp, "%d %.6e %.6e  %.6e \n",TIME_DISC, mesh->tauInv, mesh->dt, mesh->maxErrorBoltzmann);
  fclose(fp);

 
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
