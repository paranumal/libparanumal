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
  // time        = LSERK, LSIMEX, SARK3, SAAB3, MRAB, MRSAAB
  // out         = REPORT, REPORT-VTU, NO  
  // bc          = UNSPLITPML, SPLITPML, NONE
  // pmlprofile  = CONSTANT, QUADRATIC
  
  char options[BUFSIZ];
  strcpy(options,"out = REPORT + VTU+PROBE, MR_GROUPS, relaxation = COLLOCATION, bc=PML, pmlprofile=QUADRATIC");
  
    int N, time_disc;
    char meshfile[BUFSIZ]; 
    char time_str[BUFSIZ];
    mesh2D *mesh; 

    if(argc==1){
      // int specify polynomial degree 
      N = 3;
      printf("Warning!!! N is not speciefied, defaulting N=%d\n",N);
      // set up mesh stuff
      sprintf(meshfile, "../../meshes/boltzmannSquareCylinderPML2D.msh"); 

      printf("Warning!!! meshfile is not speciefied, defaulting %s\n",meshfile);  
      mesh = meshSetupTri2D(meshfile, N);  
      
      strcpy(time_str,", time= LSERK"); 
      printf("Warning!!! time discretization is not speciefied, defaulting %s\n",time_str); 
      strcat(options,time_str);
    }  
    
    if(argc==2){
      // int specify polynomial degree 
      N = 3;
      printf("Warning!!!!!! N is not speciefied, defaulting N=%d\n",N);
      // set up mesh stuff   
      mesh = meshSetupTri2D(argv[1], N); 

      strcpy(time_str,", time= LSERK"); 
      printf("Warning!!! time discretization is not speciefied, defaulting %s\n",time_str); 
      strcat(options,time_str); 
    }  

    if(argc==3){
      // int specify polynomial degree 
      N = atoi(argv[2]);
      // set up mesh stuff   
      mesh = meshSetupTri2D(argv[1], N); 

      strcpy(time_str,", time= LSERK"); 
      printf("Warning!!! time discretization is not speciefied, defaulting %s\n",time_str); 
      strcat(options,time_str); 
    }  

     if(argc==4){
      // int specify polynomial degree 
      N = atoi(argv[2]);
      // set up mesh stuff   
      mesh = meshSetupTri2D(argv[1], N); 

      strcpy(time_str,", time="); 
      strcat(time_str,argv[3]); 
      strcat(options,time_str);
    }  

    if(argc>4){
       printf("usage: ./main meshes/sample.msh N time_disc\n");
       exit(-1);
    }  
    
    for(iint i=0; i<1;i++){

      mesh->Ntscale=i;

      printf("Setup Boltzmann Solver: \n");
      boltzmannSetup2D(mesh,options); 

      printf("Boltzmann Run: \n");
      boltzmannRun2D(mesh,options);  
    }      
    
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
