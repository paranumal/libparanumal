#include "boltzmann3D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  // SET OPTIONS
  // mode        = TEST, SOLVER // do not use test mode, for developing purposes
  // relaxation  = CUBATURE, COLLOCATION, 
  // time        = LSERK, LSIMEX, SARK3, SAAB3
  // out         = REPORT, REPORT-VTU, NO  
  // bc          = PML, NONE
  
   char *options =strdup("mode = SOLVER,out = REPORT-VTU,relaxation=CUBATURE,stab=NO,time=LSERK,bc=PML");
  

    if(argc!=3){
          printf("usage: ./main meshes/cavityH005.msh N\n");
          exit(-1);
    }

    // int specify polynomial degree 
    int N = atoi(argv[2]);
    // set up mesh stuff
    mesh3D *mesh = meshSetupTet3D(argv[1], N);
        
    boltzmannSetup3D(mesh,options);   
    boltzmannRun3D(mesh,options);
   
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
