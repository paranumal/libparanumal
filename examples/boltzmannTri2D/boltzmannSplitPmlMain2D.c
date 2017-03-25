#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  //




  #if(TEST_MODE==1)

  if(argc!=4){
    // to run cavity test case with degree N elements
    printf("usage: ./main meshes/cavityH005.msh N dtfactor\n");
    exit(-1);
  }

  // int specify polynomial degree 
  int N = atoi(argv[2]);

  dfloat dtfactor = atoi(argv[3]);
  

  // set up mesh stuff
  mesh2D *mesh = meshSetupTri2D(argv[1], N);

  mesh->dtfactor  = dtfactor; 
  
  // set up boltzmann stuff
  void boltzmannSplitPmlSetup2D(mesh2D *mesh);
  void boltzmannSplitPmlRun2D(mesh2D *mesh);
  

  // run
  printf("occa run: \n");

  FILE *fp;
  char fname[BUFSIZ];
  sprintf(fname, "TimeStepSize.txt");
  
  boltzmannSplitPmlSetup2D(mesh);   
  boltzmannSplitPmlRun2D(mesh);
   

  fp = fopen(fname, "a");
  fprintf(fp, "Time steping method, Order, Element Number, TauInv, dt , Error\n");
  fprintf(fp, "%d %d %d %.6e %.6e  %.6e \n",TIME_DISC, N, mesh->Nelements, mesh->tauInv, mesh->dt, mesh->maxErrorBoltzmann);
  fclose(fp);


  #else
      if(argc!=3){
        // to run test case with degree N elements
        printf("usage: ./main meshes/cavityH005.msh N\n");
        exit(-1);
      }

      // int specify polynomial degree 
      int N = atoi(argv[2]);

      dfloat  dtfactor = 1;

      // set up mesh stuff
      mesh2D *mesh = meshSetupTri2D(argv[1], N);


      // set up boltzmann stuff
      void boltzmannSplitPmlSetup2D(mesh2D *mesh);
      void boltzmannSplitPmlRun2D(mesh2D *mesh);


      // run
      //  boltzmannRun2D(mesh);
      printf("occa run: \n");

      mesh->dtfactor = dtfactor; 

      boltzmannSplitPmlSetup2D(mesh);   
      boltzmannSplitPmlRun2D(mesh);
       
  #endif

 
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
