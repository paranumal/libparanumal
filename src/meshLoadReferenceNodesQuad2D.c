#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void readDfloatArray(FILE *fp, const char *label, dfloat **A, int *Nrows, int* Ncols);
void readIntArray   (FILE *fp, const char *label, int **A   , int *Nrows, int* Ncols);

void meshLoadReferenceNodesQuad2D(mesh2D *mesh, int N){

  char fname[BUFSIZ];
  sprintf(fname, DHOLMES "/nodes/quadrilateralN%02d.dat", N);

  FILE *fp = fopen(fname, "r");

  char buf[BUFSIZ];
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  int Ncheck;
  sscanf(buf, "%d", &Ncheck);
  if(Ncheck != N) printf("bu55er - wrong data file\n");
  mesh->N = N;
  mesh->Nq = N+1;
  mesh->Nfp = N+1;

  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  int Npcheck;
  sscanf(buf, "%d", &Npcheck);
  mesh->Np = Npcheck;

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->r = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  mesh->s = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, dfloatFormat dfloatFormat, mesh->r+n, mesh->s+n);
  }

  // find node indices of vertex nodes
  dfloat NODETOL = 1e-6;
  mesh->vertexNodes = (int*) calloc(mesh->Nverts, sizeof(int));
  for(int n=0;n<mesh->Np;++n){
    if( (mesh->r[n]+1)*(mesh->r[n]+1)+(mesh->s[n]+1)*(mesh->s[n]+1)<NODETOL)
      mesh->vertexNodes[0] = n;
    if( (mesh->r[n]-1)*(mesh->r[n]-1)+(mesh->s[n]+1)*(mesh->s[n]+1)<NODETOL)
      mesh->vertexNodes[1] = n;
    if( (mesh->r[n]-1)*(mesh->r[n]-1)+(mesh->s[n]-1)*(mesh->s[n]-1)<NODETOL)
      mesh->vertexNodes[2] = n;
    if( (mesh->r[n]+1)*(mesh->r[n]+1)+(mesh->s[n]-1)*(mesh->s[n]-1)<NODETOL)
      mesh->vertexNodes[3] = n;
  }
  
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->Dr = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->Dr+n);
  }
  fgets(buf, BUFSIZ, fp); // read comment

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->Ds = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->Ds+n);
  }
  fgets(buf, BUFSIZ, fp); // read EOL

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->faceNodes = (int*) calloc(mesh->Nfp*mesh->Nfaces, sizeof(int));
  for(int f=0;f<mesh->Nfaces;++f){
    for(int n=0;n<mesh->Nfp;++n){
      fscanf(fp, "%d", mesh->faceNodes+n + f*mesh->Nfp);
    }
  }
  fgets(buf, BUFSIZ, fp); // read EOL

    
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->LIFT = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Nfaces*mesh->Nfp*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->LIFT+n);
  }
  fgets(buf, BUFSIZ, fp);

  /* 1D collocation differentiation matrix on GLL nodes */
  fgets(buf, BUFSIZ, fp); // read comment

  mesh->D = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->N+1;++n){
    for(int m=0;m<mesh->N+1;++m){
      fscanf(fp, dfloatFormat, mesh->D+m+n*(mesh->N+1));
    }
  }
  fgets(buf, BUFSIZ, fp); // read comment

  /* 1D GLL node coordinates */
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->gllz = (dfloat*) calloc(mesh->N+1, sizeof(dfloat));
  for(int n=0;n<mesh->N+1;++n){
    fscanf(fp, dfloatFormat, mesh->gllz+n);
  }
  fgets(buf, BUFSIZ, fp); // read comment
  
  /* 1D GLL node coordinates */
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->gllw = (dfloat*) calloc(mesh->N+1, sizeof(dfloat));
  for(int n=0;n<mesh->N+1;++n){
    fscanf(fp, dfloatFormat, mesh->gllw+n);
  }
  fgets(buf, BUFSIZ, fp); // read comment

  // read number of plot nodes
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); 
  sscanf(buf, "%d", &(mesh->plotNp));

  // read plot node coordinates (hard code triangles)
  mesh->plotR = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  mesh->plotS = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->plotNp;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, dfloatFormat dfloatFormat, mesh->plotR+n, mesh->plotS+n);
  }
  
  // read plot interpolation matrix
  mesh->plotInterp = (dfloat*) calloc(mesh->plotNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->plotNp;++n){
    for(int m=0;m<mesh->Np;++m){
      fscanf(fp, dfloatFormat, mesh->plotInterp+n*mesh->Np+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // read number of elements in plot node triangulation
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); 
  sscanf(buf, "%d", &(mesh->plotNelements));

  // read number of vertices per plot element
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, "%d", &(mesh->plotNverts));
  
  // build and read in plot node triangulation
  mesh->plotEToV = (int*) calloc(mesh->plotNelements*mesh->plotNverts, sizeof(int));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->plotNelements;++n){
    for(int m=0;m<mesh->plotNverts;++m){
      fscanf(fp, "%d", mesh->plotEToV+m + mesh->plotNverts*n);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // projection info for OAS precon (one node overlap)
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, "%d", &(mesh->NqP));
  fgets(buf, BUFSIZ, fp);
  mesh->oasForward = (dfloat*) calloc(mesh->NqP*mesh->NqP, sizeof(dfloat));
  for(int n=0;n<mesh->NqP;++n){
    for(int m=0;m<mesh->NqP;++m){
      fscanf(fp, dfloatFormat, mesh->oasForward+n*mesh->NqP+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  fgets(buf, BUFSIZ, fp);
  mesh->oasDiagOp = (dfloat*) calloc(mesh->NqP, sizeof(dfloat));
  for(int n=0;n<mesh->NqP;++n){
    fscanf(fp, dfloatFormat, mesh->oasDiagOp+n);
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  fgets(buf,BUFSIZ,fp); // rest of line
  mesh->oasBack = (dfloat*) calloc(mesh->NqP*mesh->NqP, sizeof(dfloat));
  for(int n=0;n<mesh->NqP;++n){
    for(int m=0;m<mesh->NqP;++m){
      fscanf(fp, dfloatFormat, mesh->oasBack+n*mesh->NqP+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // projection info for OAS precon (one node overlap)
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, "%d", &(mesh->NqP));
  fgets(buf, BUFSIZ, fp);
  mesh->oasForwardDg = (dfloat*) calloc(mesh->NqP*mesh->NqP, sizeof(dfloat));
  for(int n=0;n<mesh->NqP;++n){
    for(int m=0;m<mesh->NqP;++m){
      fscanf(fp, dfloatFormat, mesh->oasForwardDg+n*mesh->NqP+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  fgets(buf, BUFSIZ, fp);
  mesh->oasDiagOpDg = (dfloat*) calloc(mesh->NqP, sizeof(dfloat));
  for(int n=0;n<mesh->NqP;++n){
    fscanf(fp, dfloatFormat, mesh->oasDiagOpDg+n);
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  fgets(buf,BUFSIZ,fp); // rest of line
  mesh->oasBackDg = (dfloat*) calloc(mesh->NqP*mesh->NqP, sizeof(dfloat));
  for(int n=0;n<mesh->NqP;++n){
    for(int m=0;m<mesh->NqP;++m){
      fscanf(fp, dfloatFormat, mesh->oasBackDg+n*mesh->NqP+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  
  // read number of volume cubature nodes
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); 
  sscanf(buf, "%d", &(mesh->cubNp));

  // read volume cubature interpolation matrix
  mesh->cubInterp = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->cubNp;++n){
    for(int m=0;m<mesh->Np;++m){
      fscanf(fp, dfloatFormat, mesh->cubInterp+n*mesh->Np+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // read cubature weak 'r' differentiation matrix
  mesh->cubDrW = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->cubNp;++m){
      fscanf(fp, dfloatFormat, mesh->cubDrW+n*mesh->cubNp+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }
  // read cubature weak 's' differentiation matrix
  mesh->cubDsW = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->cubNp;++m){
      fscanf(fp, dfloatFormat, mesh->cubDsW+n*mesh->cubNp+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

    // read cubature projection matrix
  mesh->cubProject = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->cubNp;++m){
      fscanf(fp, dfloatFormat, mesh->cubProject+n*mesh->cubNp+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }


  // read number of surface integration nodes
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); 
  sscanf(buf, "%d", &(mesh->intNfp));

  // read surface intergration node interpolation matrix
  mesh->intInterp 
    = (dfloat*) calloc(mesh->intNfp*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->intNfp*mesh->Nfaces;++n){
    for(int m=0;m<mesh->Nfp;++m){
      fscanf(fp, dfloatFormat, mesh->intInterp+n*mesh->Nfp+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // read lift matrix from surface integration to volume nodes
  mesh->intLIFT = (dfloat*) calloc(mesh->intNfp*mesh->Nfaces*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment

  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->intNfp*mesh->Nfaces;++m){
      fscanf(fp, dfloatFormat, mesh->intLIFT+n*mesh->intNfp*mesh->Nfaces+m);
      //      printf("%g ", mesh->intLIFT[n*mesh->intNfp*mesh->Nfaces+m]);
    }
    //    printf("\n");
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  fclose(fp);
}

void readDfloatArray(FILE *fp, const char *label, dfloat **A, int *Nrows, int* Ncols){

  char buf[BUFSIZ];
  rewind(fp); // rewind to beginning

  int flag = 0;
  int status;
  char *stat;

  //search for label in file
  while(fgets(buf, BUFSIZ, fp)){
    if (strstr(buf, label)) {flag =1; break;};
  }
  if (flag==0) {
    printf("ERROR: Unable to find label: '%s' in node file.\n", label);
    exit(-1);
  }

  //if found read in number of rows and columns
  status = fscanf(fp, "%d %d",  Nrows, Ncols);
  stat = fgets(buf, BUFSIZ, fp); //read to end of line

  *A = (dfloat*) calloc((*Nrows)*(*Ncols), sizeof(dfloat)); //allocate space
  for(int n=0;n<(*Nrows)*(*Ncols);++n) //read matrix data
    status = fscanf(fp, dfloatFormat, (*A)+n);
}

void readIntArray(FILE *fp, const char *label, int **A, int *Nrows, int* Ncols){

  char buf[BUFSIZ];
  rewind(fp); // rewind to beginning

  int flag = 0;
  int status;
  char *stat;

  //search for label in file
  while(fgets(buf, BUFSIZ, fp)){
    if (strstr(buf, label)) {flag =1; break;};
  }
  if (flag==0) {
    printf("ERROR: Unable to find label: '%s' in node file.\n", label);
    exit(-1);
  }

  //if found read in number of rows and columns
  status = fscanf(fp, "%d %d",  Nrows, Ncols);
  stat = fgets(buf, BUFSIZ, fp); //read to end of line

  *A = (int*) calloc((*Nrows)*(*Ncols), sizeof(int)); //allocate space
  for(int n=0;n<(*Nrows)*(*Ncols);++n) //read matrix data
    status = fscanf(fp, "%d", (*A)+n);
}