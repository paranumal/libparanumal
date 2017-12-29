#include "ellipticTri2D.h"
void loadElementStiffnessMatricesTri2D(mesh2D *mesh, const char *options, int N){

  char fname[BUFSIZ];
  sprintf(fname, "sparseN%02d.dat", N);
  FILE *fp = fopen(fname, "r");
  printf("opened\n");

  char buf[BUFSIZ];
  printf("buffer created\n");
  if (fp == NULL) {
    printf("no file! %s\n", fname);
    exit(-1);
  }

  fgets(buf, BUFSIZ, fp); // read comment
  printf("Comment 1%s ", buf);

  fgets(buf, BUFSIZ, fp);
  int maxNnzPerRow;
  sscanf(buf, "%d", &maxNnzPerRow);
  printf("maxNNz = %d\n", maxNnzPerRow);

  mesh->maxNnzPerRow = maxNnzPerRow;
  fgets(buf, BUFSIZ, fp); // another  comment
  printf("Comment 2%s ", buf);
  mesh->Ind = (iint*) calloc(mesh->maxNnzPerRow*mesh->maxNnzPerRow, sizeof(iint));                     
  for(int n=0;n<mesh->maxNnzPerRow*mesh->maxNnzPerRow;++n){                                           
    fscanf(fp, "%d ", mesh->Ind+n);                                                
    //fscanf(fp, " %s ", buf);    
    printf(" %d ", mesh->Ind[n]);
    //printf("buf %s ", buf); 
    if ((n+1)%mesh->maxNnzPerRow == 0){printf(" \n ");}
  }   
  mesh->Srr =  (dfloat*) calloc(mesh->maxNnzPerRow*mesh->maxNnzPerRow, sizeof(dfloat));
  mesh->Srs =  (dfloat*) calloc(mesh->maxNnzPerRow*mesh->maxNnzPerRow, sizeof(dfloat));
  mesh->Sss =  (dfloat*) calloc(mesh->maxNnzPerRow*mesh->maxNnzPerRow, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp);
  //fscanf(fp, "%s", buf);
  printf("Comment 3 %s ", buf);
  for(int n=0;n<mesh->maxNnzPerRow*mesh->maxNnzPerRow;++n){                                           
    char test[1000];    
    int aa = fscanf(fp,  "%lf ", mesh->Srr+n);                                                
    printf(" %f ", mesh->Srr[n]);
    //mesh->Srr[n]); 
    if ((n+1)%mesh->maxNnzPerRow == 0){
      printf(" \n ");
    }

  }   
  fgets(buf, BUFSIZ, fp);

  for(int n=0;n<mesh->maxNnzPerRow*mesh->maxNnzPerRow;++n){                                           
    fscanf(fp, "%lf ", mesh->Srs+n);
    printf(" %f ", mesh->Srs[n]);
    //mesh->Srr[n]); 
    if ((n+1)%mesh->maxNnzPerRow == 0){
      printf(" \n ");
    }
                                                
  } 
  fgets(buf, BUFSIZ, fp);

  for(int n=0;n<mesh->maxNnzPerRow*mesh->maxNnzPerRow;++n){                                           
    fscanf(fp, "%lf ", mesh->Sss+n);                                                
    printf(" %f ", mesh->Sss[n]);
    //mesh->Srr[n]); 
    if ((n+1)%mesh->maxNnzPerRow == 0){
      printf(" \n ");
    }

  }   

  //now occa copy, transpose and stuff
  dfloat *SrrT, *SrsT, *SsrT, *SssT;
  iint *IndT;  
  SrrT = (dfloat *) calloc(mesh->maxNnzPerRow*mesh->maxNnzPerRow,sizeof(dfloat));
  SrsT = (dfloat *) calloc(mesh->maxNnzPerRow*mesh->maxNnzPerRow,sizeof(dfloat));
  SssT = (dfloat *) calloc(mesh->maxNnzPerRow*mesh->maxNnzPerRow,sizeof(dfloat));
  IndT = (iint*)   calloc(mesh->maxNnzPerRow*mesh->maxNnzPerRow,sizeof(iint));
  for (iint n=0;n<mesh->maxNnzPerRow;n++) {
    for (iint m=0;m<mesh->maxNnzPerRow;m++) {  
      SrrT[m+n*mesh->maxNnzPerRow] = mesh->Srr[n+m*mesh->maxNnzPerRow];
      SrsT[m+n*mesh->maxNnzPerRow] = mesh->Srs[n+m*mesh->maxNnzPerRow];
      SssT[m+n*mesh->maxNnzPerRow] = mesh->Sss[n+m*mesh->maxNnzPerRow];
      IndT[m+n*mesh->maxNnzPerRow] = mesh->Ind[n+m*mesh->maxNnzPerRow];
    }
  }
  mesh->o_SrrT = mesh->device.malloc(mesh->maxNnzPerRow*mesh->maxNnzPerRow*sizeof(dfloat), SrrT);
  mesh->o_SrsT = mesh->device.malloc(mesh->maxNnzPerRow*mesh->maxNnzPerRow*sizeof(dfloat), SrsT);
  mesh->o_SssT = mesh->device.malloc(mesh->maxNnzPerRow*mesh->maxNnzPerRow*sizeof(dfloat), SssT);
  mesh->o_IndT = mesh->device.malloc(mesh->maxNnzPerRow*mesh->maxNnzPerRow*sizeof(iint), IndT);

}
