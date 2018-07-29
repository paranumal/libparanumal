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

#include "ellipticBenchmarkTet3D.h"

typedef struct {
  char x,y,z,w;
}char4;

void loadElementStiffnessMatricesTet3D(mesh_t *mesh, const char *options, int N){

  char fname[BUFSIZ];
  sprintf(fname, "sparseMatrices/sparseTet3DN%02d.dat", N);
  FILE *fp = fopen(fname, "r");

  char buf[BUFSIZ];
  if (fp == NULL) {
    printf("no file! %s\n", fname);
    exit(-1);
  }

  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  int maxNnzPerRow;
  sscanf(buf, "%d", &maxNnzPerRow);
  printf("maxNNz = %d Np = %d \n", maxNnzPerRow, mesh->Np);

  mesh->maxNnzPerRow = maxNnzPerRow;
  int paddedRowSize = 4*((mesh->maxNnzPerRow+3)/4);
  printf("maxNnzPerRow = %d, paddedNnzPerRow = %d\n", mesh->maxNnzPerRow, paddedRowSize);

  fgets(buf, BUFSIZ, fp); 
  mesh->Ind = (int*) calloc(paddedRowSize*mesh->Np, sizeof(int));                     
  for(int n=0;n<mesh->maxNnzPerRow*mesh->Np;++n){ 
    fscanf(fp, "%d ", mesh->Ind+n);                                             
  }   
  mesh->Srr =  (dfloat*) calloc(paddedRowSize*mesh->Np, sizeof(dfloat));
  mesh->Srs =  (dfloat*) calloc(paddedRowSize*mesh->Np, sizeof(dfloat));
  mesh->Srt =  (dfloat*) calloc(paddedRowSize*mesh->Np, sizeof(dfloat));
  mesh->Sss =  (dfloat*) calloc(paddedRowSize*mesh->Np, sizeof(dfloat));
  mesh->Sst =  (dfloat*) calloc(paddedRowSize*mesh->Np, sizeof(dfloat));
  mesh->Stt =  (dfloat*) calloc(paddedRowSize*mesh->Np, sizeof(dfloat));
  mesh->MM =  (dfloat*) calloc(paddedRowSize*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp);
  for(int r=0;r<mesh->maxNnzPerRow;++r){
    for(int n=0;n<mesh->Np;++n){                                           
      int count = n + r*mesh->Np;
      int aa = fscanf(fp,  "%lf ", mesh->Srr+count); 
    }
  }
  fgets(buf, BUFSIZ, fp);
  for(int r=0;r<mesh->maxNnzPerRow;++r){
    for(int n=0;n<mesh->Np;++n){                                           
      int count = n + r*mesh->Np;
      int aa = fscanf(fp,  "%lf ", mesh->Srs+count); 
    }
  }
  fgets(buf, BUFSIZ, fp);
  for(int r=0;r<mesh->maxNnzPerRow;++r){
    for(int n=0;n<mesh->Np;++n){                                           
      int count = n + r*mesh->Np;
      int aa = fscanf(fp,  "%lf ", mesh->Srt+count); 
    }
  }
  fgets(buf, BUFSIZ, fp);
  for(int r=0;r<mesh->maxNnzPerRow;++r){
    for(int n=0;n<mesh->Np;++n){                                           
      int count = n + r*mesh->Np;
      int aa = fscanf(fp,  "%lf ", mesh->Sss+count); 
    }
  }
  fgets(buf, BUFSIZ, fp);
  for(int r=0;r<mesh->maxNnzPerRow;++r){
    for(int n=0;n<mesh->Np;++n){                                           
      int count = n + r*mesh->Np;
      int aa = fscanf(fp,  "%lf ", mesh->Sst+count); 
    }
  }
  fgets(buf, BUFSIZ, fp);
  for(int r=0;r<mesh->maxNnzPerRow;++r){
    for(int n=0;n<mesh->Np;++n){                                           
      int count = n + r*mesh->Np;
      int aa = fscanf(fp,  "%lf ", mesh->Stt+count); 
    }
  }
  fgets(buf, BUFSIZ, fp);
  
  for(int r=0;r<mesh->maxNnzPerRow;++r){
    for(int n=0;n<mesh->Np;++n){                                           
      int count = n + r*mesh->Np;
      int aa = fscanf(fp,  "%lf ", mesh->MM+count); 
    }
  }
  fgets(buf, BUFSIZ, fp);
  
  /*char4 packing*/
  char4 * IndTchar = (char4*) calloc((paddedRowSize*mesh->Np)/4,sizeof(char4));
  for (int m=0;m<paddedRowSize;m+=4) {  
    for (int n=0;n<mesh->Np;n++) {
      char4 ind;
      ind.x = mesh->Ind[n + 0*mesh->Np + m*mesh->Np];
      ind.y = mesh->Ind[n + 1*mesh->Np + m*mesh->Np];
      ind.z = mesh->Ind[n + 2*mesh->Np + m*mesh->Np];
      ind.w = mesh->Ind[n + 3*mesh->Np + m*mesh->Np];
      
      IndTchar[n + mesh->Np*(m/4)] = ind;
    }
  }

  int nnz = 0;
  for(int n=0;n<paddedRowSize*mesh->Np;++n)
    nnz += (mesh->Ind[n]!=0);
  printf("NNZ = %d\n", nnz);
  
  mesh->maxNnzPerRow = paddedRowSize;
  mesh->o_SrrT = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Srr);
  mesh->o_SrsT = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Srs);
  mesh->o_SrtT = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Srt);
  mesh->o_SssT = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Sss);
  mesh->o_SstT = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Sst);
  mesh->o_SttT = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Stt);

  mesh->o_Srr = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Srr);
  mesh->o_Srs = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Srs);
  mesh->o_Srt = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Srt);
  mesh->o_Sss = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Sss);
  mesh->o_Sst = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Sst);
  mesh->o_Stt = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Stt);

  mesh->o_MM = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->MM);

  mesh->o_IndTchar = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(char), IndTchar);
  
  fclose(fp); 
}
