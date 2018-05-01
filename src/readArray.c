#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

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
