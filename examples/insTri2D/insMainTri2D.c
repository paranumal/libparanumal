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
  char *options = strdup("out=REPORT+VTU, adv=CUBATURE, disc = DISCONT_GALERKIN"); // SUBCYCLING
  //  char *options = strdup("out=REPORT+VTU, adv=COLLOCATION, disc = DISCONT_GALERKIN");
 
  char *velSolverOptions =
    strdup("solver=PCG method=IPDG preconditioner=BLOCKJACOBI");

  char *prSolverOptions =
    strdup("solver=PCG,FLEXIBLE method=IPDG,PROJECT preconditioner=FULLALMOND,MATRIXFREE"); // ,FORCESYMMETRY"); // ,FORCESYMMETRY");

  if(argc!=3 && argc!=4){
    printf("usage 1: ./main meshes/cavityH005.msh N\n");
    printf("usage 2: ./main meshes/cavityH005.msh N insUniformFlowBoundaryConditions.h\n");
    exit(-1);
  }
  // int specify polynomial degree
  int N = atoi(argv[2]);

  // set up mesh stuff
  mesh2D *mesh = meshSetupTri2D(argv[1], N); 

  // capture header file
  char *boundaryHeaderFileName;
  if(argc==3)
    boundaryHeaderFileName = strdup(DHOLMES "/examples/insTri2D/insUniform2D.h"); // default
  else
    boundaryHeaderFileName = strdup(argv[3]);
 
  printf("Setup INS Solver: \n");
  ins_t *ins = insSetup2D(mesh,options,velSolverOptions,prSolverOptions,boundaryHeaderFileName);

  printf("OCCA Run: \n");
  insRun2D(ins,options);



  // char fname[BUFSIZ];
  //  sprintf(fname, "Int.txt");
  //  FILE *fp;
  //  fp = fopen(fname, "w");
  // mesh->o_intx.copyTo(mesh->intx);// = (dfloat*) calloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  // mesh->o_inty.copyTo(mesh->inty);// = (dfloat*) calloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  // for(iint e=0;e<mesh->Nelements;++e){
  //   for(iint f=0;f<mesh->Nfaces;++f){
  //     for(iint n=0;n<mesh->intNfp;++n){
  //       iint id = n + f*mesh->intNfp + e*mesh->Nfaces*mesh->intNfp;
  //       fprintf(fp,"%.15e %.15e\n",  mesh->intx[id],mesh->inty[id]);
        
  //     }
  //   }
  // }


  // for(iint e=0;e<mesh->Nelements;++e){
  //   for(iint f=0;f<mesh->Nfaces;++f){
  //       for(iint m=0;m<mesh->Nfp;++m){
  //         iint vid = mesh->vmapM[m+f*mesh->Nfp+e*mesh->Nfp*mesh->Nfaces];
  //         dfloat xm = mesh->x[vid];
  //         dfloat ym = mesh->y[vid];
  //          fprintf(fp,"%.15e %.15e\n",xm,ym);
  //       }
  //     }
  //   }
  
  // fprintf(fp,"\n");
  //  //
  // for(int n=0;n<mesh->Np;++n){
  // for(int m=0;m<mesh->intNfp*mesh->Nfaces;++m){
  //     dfloat opt = mesh->intLIFT[n+m*mesh->Np];
  //     fprintf(fp,"%.15e ", opt);
  //   }
  //  fprintf(fp,"\n");
  // }
  
  // fprintf(fp,"\n");

 
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
