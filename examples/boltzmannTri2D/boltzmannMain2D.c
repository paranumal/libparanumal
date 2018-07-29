/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

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

  if(argc!=2){
    printf("usage2: ./boltzmannMainTri2D setupfile\n");
    exit(-1);
  }

  setupAide options(argv[1]);

  string fileName;
  int N;

  options.getArgs("MESH FILE", fileName);
  options.getArgs("POLYNOMIAL DEGREE", N);
  mesh2D *mesh = meshSetupTri2D((char*)fileName.c_str(), N);


  bns_t *bns = boltzmannSetup2D(mesh,options); 


  boltzmannRun2D(bns,options);



  #if 0
  char options[BUFSIZ];
  strcpy(options,"out = REPORT+ VTU+PROBE, MR_GROUPS, relaxation = CUBATURE, bc=PML, pmlprofile=FORTHORDER");
  
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

    int SCALE[4]; 

    SCALE[0] = 150; SCALE[1] = 300;  SCALE[2] = 10;   SCALE[3] = 20;   
    
    for(int i=0; i<1;i++){

      // mesh->Ntscale=SCALE[i];
      printf("Setup Boltzmann Solver: \n");
      bns_t *bns = boltzmannSetup2D(mesh,options); 

      bns->Ntscale=i;

      printf("Boltzmann Run: \n");
      boltzmannRun2D(bns,options);  
    } 

  #endif     
    
  // close down MPI
  MPI_Finalize();
  exit(0);
  return 0;
}
