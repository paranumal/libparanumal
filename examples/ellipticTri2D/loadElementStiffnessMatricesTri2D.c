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
  printf("maxNNz = %d Np = %d \n", maxNnzPerRow, mesh->Np);

  mesh->maxNnzPerRow = maxNnzPerRow;
  int paddedRowSize = mesh->maxNnzPerRow+ (4-(mesh->maxNnzPerRow%4));
  fgets(buf, BUFSIZ, fp); // another  comment
  printf("Comment 2%s ", buf);
  mesh->Ind = (iint*) calloc(mesh->maxNnzPerRow*mesh->Np, sizeof(iint));                     
  for(int n=0;n<mesh->maxNnzPerRow*mesh->Np;++n){                                           
    fscanf(fp, "%d ", mesh->Ind+n);                                                
    //fscanf(fp, " %s ", buf);    
    printf(" %d ", mesh->Ind[n]);
    //printf("buf %s ", buf); 
    if ((n+1)%mesh->Np == 0){printf(" \n ");}
  }   
  mesh->Srr =  (dfloat*) calloc(paddedRowSize*mesh->Np, sizeof(dfloat));
  mesh->Srs =  (dfloat*) calloc(paddedRowSize*mesh->Np, sizeof(dfloat));
  mesh->Sss =  (dfloat*) calloc(paddedRowSize*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp);
  //printf("allocating %d floats, mesh->Np = %d and paddedRowSize = %d \n", paddedRowSize*mesh->Np, mesh->Np, paddedRowSize);
  //fscanf(fp, "%s", buf);
  printf("Comment 3 %s ", buf);
  printf("Srr :\n");
  int count = -1;
  for(int n=0;n<mesh->maxNnzPerRow*mesh->Np;++n){                                           
    count++;
    int aa = fscanf(fp,  "%lf ", mesh->Srr+count);                                                

    printf(" %f ", mesh->Srr[count]);
    if ((n+1)%mesh->maxNnzPerRow == 0){
      while ((count+1)%paddedRowSize != 0){
        count++;
        mesh->Srr[count] = 0.0;
        printf(" %f ", mesh->Srr[count]);
      }
      printf(" \n ");
    }
  }   
  printf("\nSrs \n");


  fgets(buf, BUFSIZ, fp);
  count = -1;
  for(int n=0;n<mesh->maxNnzPerRow*mesh->Np;++n){                                           
    count++;
    int aa = fscanf(fp,  "%lf ", mesh->Srs+count);

    printf(" %f ", mesh->Srs[count]);
    if ((n+1)%mesh->maxNnzPerRow == 0){
      while ((count+1)%paddedRowSize != 0){
        count++;
        mesh->Srs[count] = 0.0;
        printf(" %f ", mesh->Srs[count]);
      }
      printf(" \n ");
    }  
  }


  fgets(buf, BUFSIZ, fp);
  printf("\n Sss \n");
  count = -1;
  for(int n=0;n<mesh->maxNnzPerRow*mesh->Np;++n){                                           
    count++;
    int aa = fscanf(fp,  "%lf ", mesh->Sss+count);

    printf(" %f ", mesh->Sss[count]);
    if ((n+1)%mesh->maxNnzPerRow == 0){
      while ((count+1)%paddedRowSize != 0){
        count++;
        mesh->Sss[count] = 0.0;
        printf(" %f ", mesh->Sss[count]);
      }
      printf(" \n ");
    }
  }   
  printf("\n count = %d \n", count);

  /*char4 packing*/

  //now occa copy, transpose and stuff
  dfloat *SrrT, *SrsT, *SsrT, *SssT;
  iint *IndT;  

  printf("before alloc \n");

  SrrT = (dfloat *) calloc(paddedRowSize*mesh->Np,sizeof(dfloat));
  SrsT = (dfloat *) calloc(paddedRowSize*mesh->Np,sizeof(dfloat));
  SssT = (dfloat *) calloc(paddedRowSize*mesh->Np,sizeof(dfloat));
  printf("allocated Srr etc \n");

  IndT = (iint*)   calloc(paddedRowSize*mesh->Np,sizeof(iint));
  printf("allocated Ind\n");


  printf("\nIND transose\n");
  for (iint n=0;n<mesh->Np;n++) {
    for (iint m=0;m<paddedRowSize;m++) {  
      SrrT[m+n*paddedRowSize] = mesh->Srr[n+m*mesh->Np];
      SrsT[m+n*paddedRowSize] = mesh->Srs[n+m*mesh->Np];
      SssT[m+n*paddedRowSize] = mesh->Sss[n+m*mesh->Np];
      if (m<mesh->maxNnzPerRow){      
        IndT[m+n*paddedRowSize] = mesh->Ind[n+m*mesh->Np];
      }
      //printf("IndT[%d] =  %d \n",m+n*paddedRowSize,  IndT[m+n*paddedRowSize]);   
    }
    //printf("\n");
  }
  printf("before char packing \n");

  char * IndTchar = (char*) calloc(paddedRowSize*mesh->Np,sizeof(char));
  for (iint n=0;n<mesh->Np;n++) {
    for (iint m=0;m<paddedRowSize;m++) {  
      IndTchar[m+n*paddedRowSize] = (char) IndT[m+n*paddedRowSize];
      printf("IndT[%d] =  %d \n",m+n*paddedRowSize,  IndT[m+n*paddedRowSize]);   
    }
    printf("\n");
  }
  printf("\n\nSrrT: ");
  for (iint n=0;n<mesh->Np;n++) {
    for (iint m=0;m<paddedRowSize;m++) {  
      printf("SrsT[%d] =  %f \n",m+n*paddedRowSize,  SrrT[m+n*paddedRowSize]);   
    }
    printf("\n");
  }


  printf("before occa alloc \n");

  //ACHTUNG!!! ACHTUNG!!! ACHTUNG!!! mesh->maxNnzPerRow changes to acchount for padding!!!
  mesh->maxNnzPerRow = paddedRowSize;
  mesh->o_SrrT = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), SrrT);
  mesh->o_SrsT = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), SrsT);
  mesh->o_SssT = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), SssT);
  mesh->o_IndT = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(iint), IndT);

  mesh->o_Srr= mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Srr);
  mesh->o_Srs = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Srs);
  mesh->o_Sss = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Sss);


  mesh->o_IndTchar = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(char), IndTchar);

  fclose(fp); 
  printf("matrices loaded \n");
}
