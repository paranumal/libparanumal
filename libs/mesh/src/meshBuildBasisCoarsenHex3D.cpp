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
#include "mesh3D.hpp"

void meshHex3D::BuildBasisCoarsen(dfloat**R, occa::memory& o_R, int Nf, int Nc) {

  int NqFine   = Nf+1;
  int NqCoarse = Nc+1;
  dfloat *P    = (dfloat *) calloc(NqFine*NqCoarse,sizeof(dfloat));
  dfloat *Ptmp = (dfloat *) calloc(NqFine*NqCoarse,sizeof(dfloat));

  //initialize P as identity
  for (int i=0;i<NqCoarse;i++) P[i*NqCoarse+i] = 1.0;

  for (int n=Nc;n<Nf;n++) {

    int Nqp1 = n+2;
    int Nqp  = n+1;

    //copy P
    for (int i=0;i<Nqp*NqCoarse;i++) Ptmp[i] = P[i];

    //get the raise op from the node file
    char fname[BUFSIZ];
    sprintf(fname, LIBP_DIR "/nodes/hexN%02d.dat", n);

    FILE *fp = fopen(fname, "r");

    if (!fp) {
      stringstream ss;
      ss << "Cannot open file: " << fname;
      LIBP_ABORT(ss.str())
    }

    int Nrows, Ncols;
    dfloat *InterpRaise;
    readDfloatArray(comm, fp, "1D degree raise matrix", &(InterpRaise), &Nrows, &Ncols);

    //Multiply by the raise op
    for (int i=0;i<Nqp1;i++) {
      for (int j=0;j<NqCoarse;j++) {
        P[i*NqCoarse + j] = 0.;
        for (int k=0;k<Nqp;k++) {
          P[i*NqCoarse + j] += InterpRaise[i*Nqp+k]*Ptmp[k*NqCoarse + j];
        }
      }
    }

    fclose(fp);
    free(InterpRaise);
  }

  //the coarsen matrix is P^T
  *R = (dfloat *) calloc(NqFine*NqCoarse,sizeof(dfloat));
  for (int i=0;i<NqCoarse;i++) {
    for (int j=0;j<NqFine;j++) {
      (*R)[i*NqFine+j] = P[j*NqCoarse+i];
    }
  }
  o_R = device.malloc(NqFine*NqCoarse*sizeof(dfloat), *R);

  free(P); free(Ptmp);
}


void meshHex3D::BuildInterpolation(dfloat**R, occa::memory& o_R, int Nf, int Nc) {

  int NqFine   = Nf+1;
  int NqCoarse = Nc+1;

  dfloat *P    = (dfloat *) calloc(NqFine*NqFine,sizeof(dfloat));
  dfloat *Ptmp = (dfloat *) calloc(NqFine*NqFine,sizeof(dfloat));

  //initialize P as Nf-1 identity
  for (int i=0;i<NqFine;i++) P[i*NqFine+i] = 1.0;

  for (int n=Nf;n>Nc;n--) {

    int Nqp1 = n;
    int Nqp  = n+1;

    //copy P
    for (int i=0;i<Nqp*NqFine;i++) Ptmp[i] = P[i];

    //get the raise op from the node file
    char fname[BUFSIZ];
    sprintf(fname, LIBP_DIR "/nodes/quadrilateralN%02d.dat", n);

    FILE *fp = fopen(fname, "r");

    if (!fp) {
      stringstream ss;
      ss << "Cannot open file: " << fname;
      LIBP_ABORT(ss.str())
    }

    int Nrows, Ncols;
    dfloat *InterpLower;
    readDfloatArray(comm, fp, "1D degree lower matrix", &(InterpLower), &Nrows, &Ncols);

    //Multiply by the raise op
    for (int i=0;i<Nqp1;i++) {
      for (int j=0;j<NqFine;j++) {
        P[i*NqFine + j] = 0.;
        for (int k=0;k<Nqp;k++) {
          P[i*NqFine + j] += InterpLower[i*Nqp+k]*Ptmp[k*NqFine + j];
        }
      }
    }

    fclose(fp);
    free(InterpLower);
  }

  //copy non-zero part of the interpolation matrix
  *R = (dfloat *) calloc(NqFine*NqCoarse,sizeof(dfloat));
  for (int i=0;i<NqCoarse;i++) {
    for (int j=0;j<NqFine;j++) {
      (*R)[i*NqFine+j] = P[i*NqFine+j];
    }
  }
  o_R = device.malloc(NqFine*NqCoarse*sizeof(dfloat), *R);

  free(P); free(Ptmp);
}
