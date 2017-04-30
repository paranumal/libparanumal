#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  
  // SET OPTIONS
  // mode        = TEST, SOLVER // do not use test mode, for developing purposes
  // relaxation  = CUBATURE, COLLOCATION, 
  // time        = LSERK, LSIMEX, SARK3, SAAB3
  // out         = REPORT, REPORT-VTU, NO  
  // bc          = UNSPLITPML, SPLITPML, NONE
  
   char *options =strdup("mode = SOLVER , out = REPORT-VTU, relaxation = CUBATURE, stab = NO, time = SARK3, bc = UNSPLITPML");
  
  if(strstr(options, "SOLVER")){
    if(argc!=3){
          // to run test case with degree N elements
          printf("usage: ./main meshes/cavityH005.msh N\n");
          exit(-1);
    }

    // int specify polynomial degree 
    int N            = atoi(argv[2]);
    dfloat  dtfactor = 1;

    // set up mesh stuff
    mesh2D *mesh = meshSetupTri2D(argv[1], N);


    // set up boltzmann stuff
    void boltzmannSetup2D(mesh2D *mesh, char *options);
    void boltzmannRun2D(mesh2D *mesh, char *options);

    printf("occa run: \n");

    mesh->dtfactor = dtfactor; 

    boltzmannSetup2D(mesh,options);   
    boltzmannRun2D(mesh,options);
   

  }
   else if(strstr(options, "TEST")){
     if(argc!=4){
      // To find stable time-step sizes, give dtfactor
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
    void boltzmannSetup2D(mesh2D *mesh, char *options);
    void boltzmannRun2D(mesh2D *mesh, char *options);




    // run
    printf("occa run: \n");

    FILE *fp;
    char fname[BUFSIZ];
    sprintf(fname, "TimeStepSize.txt");

    boltzmannSetup2D(mesh,options);   
    boltzmannRun2D(mesh,options);


    // fp = fopen(fname, "a");
    // fprintf(fp, "Time steping method, Order, Element Number, TauInv, dt , Error\n");
    // fprintf(fp, "%d %d %d %.6e %.6e  %.6e \n",TIME_DISC, N, mesh->Nelements, mesh->tauInv, mesh->dt, mesh->maxErrorBoltzmann);
    // fclose(fp);

      
       
  
}
 
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
