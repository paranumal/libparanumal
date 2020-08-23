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

#include "mesh.hpp"

void readDfloatArray(MPI_Comm comm, FILE *fp, const char *label, dfloat **A, int *Nrows, int* Ncols){

  int rank;
  MPI_Comm_rank(comm, &rank);

  char buf[BUFSIZ];
  int flag = 0;

  //only root rank reads
  if (!rank) {
    rewind(fp); // rewind to beginning

    //search for label in file
    while(fgets(buf, BUFSIZ, fp)){
      if (strstr(buf, label)) {flag =1; break;};
    }
    if (flag==0) {
      stringstream ss;
      ss << "Unable to find label: " << label;
      LIBP_ABORT(ss.str())
    }

    //if found read in number of rows and columns
    if (!fscanf(fp, "%d %d",  Nrows, Ncols)) {
      stringstream ss;
      ss << "Unable to find matrix dimensions for entry: " << label;
      LIBP_ABORT(ss.str())
    }
    if (!fgets(buf, BUFSIZ, fp)) { //read to end of line
      stringstream ss;
      ss << "Error reading label: " << label;
      LIBP_ABORT(ss.str())
    }
  }

  //broadcast to the comm
  MPI_Bcast(Nrows, 1, MPI_INT, 0, comm);
  MPI_Bcast(Ncols, 1, MPI_INT, 0, comm);

  *A = (dfloat*) calloc((*Nrows)*(*Ncols), sizeof(dfloat)); //allocate space

  if (!rank) {
    for(int n=0;n<(*Nrows)*(*Ncols);++n) {//read matrix data
      if (!fscanf(fp, dfloatFormat, (*A)+n)) {
        stringstream ss;
        ss << "Error reading label: " << label;
        LIBP_ABORT(ss.str())
      }
    }
  }

  //broadcast to the comm
  MPI_Bcast(*A, (*Nrows)*(*Ncols), MPI_DFLOAT, 0, comm);
}

void readIntArray(MPI_Comm comm, FILE *fp, const char *label, int **A, int *Nrows, int* Ncols){

  int rank;
  MPI_Comm_rank(comm, &rank);

  char buf[BUFSIZ];
  int flag = 0;

  //only root rank reads
  if (!rank) {
    rewind(fp); // rewind to beginning

    //search for label in file
    while(fgets(buf, BUFSIZ, fp)){
      if (strstr(buf, label)) {flag =1; break;};
    }
    if (flag==0) {
      stringstream ss;
      ss << "Unable to find label: " << label;
      LIBP_ABORT(ss.str())
    }

    //if found read in number of rows and columns
    if (!fscanf(fp, "%d %d",  Nrows, Ncols)) {
      stringstream ss;
      ss << "Unable to find matrix dimensions for entry: " << label;
      LIBP_ABORT(ss.str())
    }
    if (!fgets(buf, BUFSIZ, fp)) { //read to end of line
      stringstream ss;
      ss << "Error reading label: " << label;
      LIBP_ABORT(ss.str())
    }
  }

  //broadcast to the comm
  MPI_Bcast(Nrows, 1, MPI_INT, 0, comm);
  MPI_Bcast(Ncols, 1, MPI_INT, 0, comm);

  *A = (int*) calloc((*Nrows)*(*Ncols), sizeof(int)); //allocate space

  if (!rank) {
    for(int n=0;n<(*Nrows)*(*Ncols);++n) {//read matrix data
      if (!fscanf(fp, "%d", (*A)+n)) {
        stringstream ss;
        ss << "Error reading label: " << label;
        LIBP_ABORT(ss.str())
      }
    }
  }

  //broadcast to the comm
  MPI_Bcast(*A, (*Nrows)*(*Ncols), MPI_INT, 0, comm);
}
