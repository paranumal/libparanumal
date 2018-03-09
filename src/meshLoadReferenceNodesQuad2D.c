#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void readDfloatArray(FILE *fp, const char *label, dfloat **A, int *Nrows, int* Ncols);
void readIntArray   (FILE *fp, const char *label, int **A   , int *Nrows, int* Ncols);

void meshLoadReferenceNodesQuad2D(mesh2D *mesh, int N){

  char fname[BUFSIZ];
  sprintf(fname, DHOLMES "/nodes/quadrilateralN%02d.dat", N);

  FILE *fp = fopen(fname, "r");

  if (!fp) {
    printf("ERROR: Cannot open file: '%s'\n", fname);
    exit(-1);
  }

  mesh->N = N;
  mesh->Nfp = N+1;
  mesh->Nq = (N+1);
  mesh->Np = (N+1)*(N+1);

  int Nrows, Ncols;

  /* Nodal Data */
  readDfloatArray(fp, "Nodal r-coordinates", &(mesh->r),&Nrows,&Ncols);
  readDfloatArray(fp, "Nodal s-coordinates", &(mesh->s),&Nrows,&Ncols);
  readDfloatArray(fp, "Nodal Dr differentiation matrix", &(mesh->Dr), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Ds differentiation matrix", &(mesh->Ds), &Nrows, &Ncols);
  readIntArray   (fp, "Nodal Face nodes", &(mesh->faceNodes), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Lift Matrix", &(mesh->LIFT), &Nrows, &Ncols);
  
  readDfloatArray(fp, "Nodal 1D GLL Nodes", &(mesh->gllz), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal 1D GLL Weights", &(mesh->gllw), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal 1D differentiation matrix", &(mesh->D), &Nrows, &Ncols);

  readDfloatArray(fp, "1D degree raise matrix", &(mesh->interpRaise), &Nrows, &Ncols);
  readDfloatArray(fp, "1D degree lower matrix", &(mesh->interpLower), &Nrows, &Ncols);

  /* Plotting data */ 
  readDfloatArray(fp, "Plotting r-coordinates", &(mesh->plotR),&Nrows,&Ncols);
  readDfloatArray(fp, "Plotting s-coordinates", &(mesh->plotS),&Nrows,&Ncols);
  mesh->plotNp = Nrows;

  readDfloatArray(fp, "Plotting Interpolation Matrix", &(mesh->plotInterp),&Nrows,&Ncols);
  readIntArray   (fp, "Plotting triangulation", &(mesh->plotEToV), &Nrows, &Ncols);
  mesh->plotNelements = Nrows;
  mesh->plotNverts = Ncols;

  /* Cubature data */ 
  readDfloatArray(fp, "Cubature r-coordinates", &(mesh->cubr),&Nrows,&Ncols);
  readDfloatArray(fp, "Cubature s-coordinates", &(mesh->cubs),&Nrows,&Ncols);
  readDfloatArray(fp, "Cubature weights", &(mesh->cubw),&Nrows,&Ncols);
  mesh->cubNp = Nrows;

  readDfloatArray(fp, "Cubature Interpolation Matrix", &(mesh->cubInterp),&Nrows,&Ncols);
  readDfloatArray(fp, "Cubature Weak Dr Differentiation Matrix", &(mesh->cubDrW),&Nrows,&Ncols);
  readDfloatArray(fp, "Cubature Weak Ds Differentiation Matrix", &(mesh->cubDsW),&Nrows,&Ncols);
  readDfloatArray(fp, "Cubature Projection Matrix", &(mesh->cubProject),&Nrows,&Ncols);
  readDfloatArray(fp, "Cubature Surface Interpolation Matrix", &(mesh->intInterp),&Nrows,&Ncols);
  mesh->intNfp = Nrows/mesh->Nfaces; //number of interpolation points per face

  readDfloatArray(fp, "Cubature Surface Lift Matrix", &(mesh->intLIFT),&Nrows,&Ncols);

  /* C0 patch data */ 
  readDfloatArray(fp, "C0 overlapping patch forward matrix", &(mesh->oasForward), &Nrows, &Ncols);   
  readDfloatArray(fp, "C0 overlapping patch diagonal scaling", &(mesh->oasDiagOp), &Nrows, &Ncols);   
  readDfloatArray(fp, "C0 overlapping patch backward matrix", &(mesh->oasBack), &Nrows, &Ncols);   
  /* IPDG patch data */ 
  readDfloatArray(fp, "IPDG overlapping patch forward matrix", &(mesh->oasForwardDg), &Nrows, &Ncols);   
  readDfloatArray(fp, "IPDG overlapping patch diagonal scaling", &(mesh->oasDiagOpDg), &Nrows, &Ncols);   
  readDfloatArray(fp, "IPDG overlapping patch backward matrix", &(mesh->oasBackDg), &Nrows, &Ncols);   
  mesh->NpP = Nrows; //overlapping patch size

  fclose(fp);

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