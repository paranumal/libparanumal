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

void meshHex3D::LoadReferenceNodes(int N_){

  char fname[BUFSIZ];
  sprintf(fname, LIBP_DIR "/nodes/hexN%02d.dat", N);

  FILE *fp = fopen(fname, "r");

  if (!fp) {
    stringstream ss;
    ss << "Cannot open file: " << fname;
    LIBP_ABORT(ss.str())
  }

  N = N_;
  Nq = N+1;
  Nfp = (N+1)*(N+1);
  Np = (N+1)*(N+1)*(N+1);
  Nverts = 8;

  int Nrows, Ncols;

  /* Nodal Data */
  readDfloatArray(fp, "Nodal r-coordinates", &(r),&Nrows,&Ncols);
  readDfloatArray(fp, "Nodal s-coordinates", &(s),&Nrows,&Ncols);
  readDfloatArray(fp, "Nodal t-coordinates", &(t),&Nrows,&Ncols);

  if (0) { //may not be present in the node file
    readDfloatArray(fp, "Nodal Dr differentiation matrix", &(Dr), &Nrows, &Ncols);
    readDfloatArray(fp, "Nodal Ds differentiation matrix", &(Ds), &Nrows, &Ncols);
    readDfloatArray(fp, "Nodal Dt differentiation matrix", &(Dt), &Nrows, &Ncols);
    readDfloatArray(fp, "Nodal Lift Matrix", &(LIFT), &Nrows, &Ncols);
  }

  readIntArray   (fp, "Nodal Face nodes", &(faceNodes), &Nrows, &Ncols);

  readDfloatArray(fp, "Nodal 1D GLL Nodes", &(gllz), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal 1D GLL Weights", &(gllw), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal 1D differentiation matrix", &(D), &Nrows, &Ncols);

  readDfloatArray(fp, "1D degree raise matrix", &(interpRaise), &Nrows, &Ncols);
  readDfloatArray(fp, "1D degree lower matrix", &(interpLower), &Nrows, &Ncols);

  /* Plotting data */
  readDfloatArray(fp, "Plotting r-coordinates", &(plotR),&Nrows,&Ncols);
  readDfloatArray(fp, "Plotting s-coordinates", &(plotS),&Nrows,&Ncols);
  readDfloatArray(fp, "Plotting t-coordinates", &(plotT),&Nrows,&Ncols);
  plotNp = Nrows;

  readDfloatArray(fp, "Plotting Interpolation Matrix", &(plotInterp),&Nrows,&Ncols);
  readIntArray   (fp, "Plotting triangulation", &(plotEToV), &Nrows, &Ncols);
  plotNelements = Nrows;
  plotNverts = Ncols;

  /* Quadrature data */
  readDfloatArray(fp, "Quadrature r-coordinates", &(cubr),&Nrows,&Ncols);
  readDfloatArray(fp, "Quadrature weights", &(cubw),&Nrows,&Ncols);
  cubNq = Nrows;
  cubNfp = cubNq*cubNq;
  cubNp = cubNq*cubNq*cubNq;

  readDfloatArray(fp, "Quadrature Interpolation Matrix", &(cubInterp),&Nrows,&Ncols);
  readDfloatArray(fp, "Quadrature Weak D Differentiation Matrix", &(cubDW),&Nrows,&Ncols);
  readDfloatArray(fp, "Quadrature Differentiation Matrix", &(cubD),&Nrows,&Ncols);
  readDfloatArray(fp, "Quadrature Projection Matrix", &(cubProject),&Nrows,&Ncols);


  if (0) { //may not be present in the node file
    readDfloatArray(fp, "Cubature r-coordinates", &(cubr),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature s-coordinates", &(cubs),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature t-coordinates", &(cubt),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature weights", &(cubw),&Nrows,&Ncols);
    cubNp = Nrows;

    readDfloatArray(fp, "Cubature Interpolation Matrix", &(cubInterp),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature Weak Dr Differentiation Matrix", &(cubDrW),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature Weak Ds Differentiation Matrix", &(cubDsW),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature Projection Matrix", &(cubProject),&Nrows,&Ncols);
  }

  if (0) { //may not be present in the node file
    readDfloatArray(fp, "Cubature Surface Interpolation Matrix", &(intInterp),&Nrows,&Ncols);
    intNfp = Nrows/Nfaces; //number of interpolation points per face

    readDfloatArray(fp, "Cubature Surface Lift Matrix", &(intLIFT),&Nrows,&Ncols);
  }
  intNfp = 0;
  intLIFT = NULL;

  /* C0 patch data */
  readDfloatArray(fp, "C0 overlapping patch forward matrix", &(oasForward), &Nrows, &Ncols);
  readDfloatArray(fp, "C0 overlapping patch diagonal scaling", &(oasDiagOp), &Nrows, &Ncols);
  readDfloatArray(fp, "C0 overlapping patch backward matrix", &(oasBack), &Nrows, &Ncols);
  /* IPDG patch data */
  readDfloatArray(fp, "IPDG overlapping patch forward matrix", &(oasForwardDg), &Nrows, &Ncols);
  readDfloatArray(fp, "IPDG overlapping patch diagonal scaling", &(oasDiagOpDg), &Nrows, &Ncols);
  readDfloatArray(fp, "IPDG overlapping patch backward matrix", &(oasBackDg), &Nrows, &Ncols);
  NpP = Nrows; //overlapping patch size

  readIntArray   (fp, "SEMFEM reference mesh", &(FEMEToV), &Nrows, &Ncols);
  NelFEM = Nrows;
  NpFEM = Np;

  readDfloatArray(fp, "Gauss Legendre 1D quadrature nodes", &(gjr), &Nrows, &Ncols);
  readDfloatArray(fp, "Gauss Legendre 1D quadrature weights", &(gjw), &Nrows, &Ncols);
  readDfloatArray(fp, "GLL to Gauss Legendre interpolation matrix", &(gjI), &Nrows, &Ncols);
  readDfloatArray(fp, "GLL to Gauss Legendre differentiation matrix", &(gjD), &Nrows, &Ncols);
  readDfloatArray(fp, "Gauss Legendre to Gauss Legendre differentiation matrix", &(gjD2), &Nrows, &Ncols);

  fclose(fp);

  // find node indices of vertex nodes
  dfloat NODETOL = 1e-6;
  vertexNodes = (int*) calloc(Nverts, sizeof(int));
  for(int n=0;n<Np;++n){
    if( (r[n]+1)*(r[n]+1)+(s[n]+1)*(s[n]+1)+(t[n]+1)*(t[n]+1)<NODETOL)
      vertexNodes[0] = n;
    if( (r[n]-1)*(r[n]-1)+(s[n]+1)*(s[n]+1)+(t[n]+1)*(t[n]+1)<NODETOL)
      vertexNodes[1] = n;
    if( (r[n]-1)*(r[n]-1)+(s[n]-1)*(s[n]-1)+(t[n]+1)*(t[n]+1)<NODETOL)
      vertexNodes[2] = n;
    if( (r[n]+1)*(r[n]+1)+(s[n]-1)*(s[n]-1)+(t[n]+1)*(t[n]+1)<NODETOL)
      vertexNodes[3] = n;
    if( (r[n]+1)*(r[n]+1)+(s[n]+1)*(s[n]+1)+(t[n]-1)*(t[n]-1)<NODETOL)
      vertexNodes[4] = n;
    if( (r[n]-1)*(r[n]-1)+(s[n]+1)*(s[n]+1)+(t[n]-1)*(t[n]-1)<NODETOL)
      vertexNodes[5] = n;
    if( (r[n]-1)*(r[n]-1)+(s[n]-1)*(s[n]-1)+(t[n]-1)*(t[n]-1)<NODETOL)
      vertexNodes[6] = n;
    if( (r[n]+1)*(r[n]+1)+(s[n]-1)*(s[n]-1)+(t[n]-1)*(t[n]-1)<NODETOL)
      vertexNodes[7] = n;
  }

  max_EL_nnz = 0;
}

