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

void meshTet3D::LoadReferenceNodes(int N_){

  char fname[BUFSIZ];
  sprintf(fname, LIBP_DIR "/nodes/tetN%02d.dat", N);

  FILE *fp = fopen(fname, "r");

  if (!fp) {
    stringstream ss;
    ss << "Cannot open file: " << fname;
    LIBP_ABORT(ss.str())
  }

  N = N_;
  Np = ((N+1)*(N+2)*(N+3))/6;
  Nfp = ((N+1)*(N+2))/2;

  int Nrows, Ncols;

  /* Nodal Data */
  readDfloatArray(fp, "Nodal r-coordinates", &(r),&Nrows,&Ncols);
  readDfloatArray(fp, "Nodal s-coordinates", &(s),&Nrows,&Ncols);
  readDfloatArray(fp, "Nodal t-coordinates", &(t),&Nrows,&Ncols);
  readDfloatArray(fp, "Nodal Dr differentiation matrix", &(Dr), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Ds differentiation matrix", &(Ds), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Dt differentiation matrix", &(Dt), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Mass Matrix", &(MM), &Nrows, &Ncols);
  readIntArray   (fp, "Nodal Face nodes", &(faceNodes), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Lift Matrix", &(LIFT), &Nrows, &Ncols);
  //readIntArray   (fp, "Nodal rotation permutations", &(rmapP), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal degree raise matrix", &(interpRaise), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal degree lower matrix", &(interpLower), &Nrows, &Ncols);

  /* Plotting data */
  readDfloatArray(fp, "Plotting r-coordinates", &(plotR),&Nrows,&Ncols);
  readDfloatArray(fp, "Plotting s-coordinates", &(plotS),&Nrows,&Ncols);
  readDfloatArray(fp, "Plotting t-coordinates", &(plotT),&Nrows,&Ncols);
  plotNp = Nrows;

  readDfloatArray(fp, "Plotting Interpolation Matrix", &(plotInterp),&Nrows,&Ncols);
  readIntArray   (fp, "Plotting triangulation", &(plotEToV), &Nrows, &Ncols);
  plotNelements = Nrows;
  plotNverts = Ncols;

  readIntArray(fp,"Contour plot EToV", &(contourEToV), &Nrows, &Ncols);
  readDfloatArray(fp,"Contour plot VX", &(contourVX), &Nrows, &Ncols);
  readDfloatArray(fp,"Contour plot VY", &(contourVY), &Nrows, &Ncols);
  readDfloatArray(fp,"Contour plot VZ", &(contourVZ), &Nrows, &Ncols);

  readDfloatArray(fp, "Contour plot Interpolation",&(contourInterp), &Nrows, &Ncols);
  readDfloatArray(fp, "Contour plot Linear Interpolation",&(contourInterp1), &Nrows, &Ncols);
  readDfloatArray(fp, "Contour plot Filter",&(contourFilter), &Nrows, &Ncols);

  /* Cubature data */
  if (N<7) {
    readDfloatArray(fp, "Cubature r-coordinates", &(cubr),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature s-coordinates", &(cubs),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature t-coordinates", &(cubt),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature weights", &(cubw),&Nrows,&Ncols);
    cubNp = Nrows;

    readDfloatArray(fp, "Cubature Interpolation Matrix", &(cubInterp),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature Weak Dr Differentiation Matrix", &(cubDrW),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature Weak Ds Differentiation Matrix", &(cubDsW),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature Weak Dt Differentiation Matrix", &(cubDtW),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature Projection Matrix", &(cubProject),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature Surface Interpolation Matrix", &(intInterp),&Nrows,&Ncols);
    intNfp = Nrows/Nfaces; //number of interpolation points per face

    readDfloatArray(fp, "Cubature Surface Lift Matrix", &(intLIFT),&Nrows,&Ncols);
  }


  /* Bernstein-Bezier data */
  readDfloatArray(fp, "Bernstein-Bezier Vandermonde Matrix", &(VB),&Nrows,&Ncols);
  readDfloatArray(fp, "Bernstein-Bezier Inverse Vandermonde Matrix", &(invVB),&Nrows,&Ncols);
  readIntArray   (fp, "Bernstein-Bezier sparse D0 differentiation ids", &(D0ids), &Nrows, &Ncols);  //Ncols should be 4
  readIntArray   (fp, "Bernstein-Bezier sparse D1 differentiation ids", &(D1ids), &Nrows, &Ncols);  //Ncols should be 4
  readIntArray   (fp, "Bernstein-Bezier sparse D2 differentiation ids", &(D2ids), &Nrows, &Ncols);  //Ncols should be 4
  readIntArray   (fp, "Bernstein-Bezier sparse D3 differentiation ids", &(D3ids), &Nrows, &Ncols);  //Ncols should be 4
  readDfloatArray(fp, "Bernstein-Bezier sparse D differentiation values", &(Dvals), &Nrows, &Ncols);//Ncols should be 4

  readIntArray   (fp, "Bernstein-Bezier sparse D0T transpose differentiation ids", &(D0Tids), &Nrows, &Ncols);  //Ncols should be 4
  readIntArray   (fp, "Bernstein-Bezier sparse D1T transpose differentiation ids", &(D1Tids), &Nrows, &Ncols);  //Ncols should be 4
  readIntArray   (fp, "Bernstein-Bezier sparse D2T transpose differentiation ids", &(D2Tids), &Nrows, &Ncols);  //Ncols should be 4
  readIntArray   (fp, "Bernstein-Bezier sparse D3T transpose differentiation ids", &(D3Tids), &Nrows, &Ncols);  //Ncols should be 4
  readDfloatArray(fp, "Bernstein-Bezier sparse DT transpose differentiation values", &(DTvals), &Nrows, &Ncols);//Ncols should be 4

  readIntArray   (fp, "Bernstein-Bezier L0 Matrix ids", &(L0ids), &Nrows, &Ncols);
  readDfloatArray(fp, "Bernstein-Bezier L0 Matrix values", &(L0vals), &Nrows, &Ncols); //Ncols should be 7
  readIntArray   (fp, "Bernstein-Bezier EL lift ids", &(ELids), &Nrows, &Ncols);
  readDfloatArray(fp, "Bernstein-Bezier EL lift values", &(ELvals), &Nrows, &Ncols);
  max_EL_nnz = Ncols;

  readIntArray   (fp, "Bernstein-Bezier sparse 2D degree raise ids", &(BBRaiseids), &Nrows, &Ncols);     //Ncols should be 3
  readDfloatArray(fp, "Bernstein-Bezier sparse 2D degree raise values", &(BBRaiseVals), &Nrows, &Ncols); //Ncols should be 3
  readDfloatArray(fp, "Bernstein-Bezier sparse 2D degree lower matrix", &(BBLower), &Nrows, &Ncols);

  /* IPDG patch data */
  readDfloatArray(fp, "IPDG overlapping patch forward matrix", &(oasForwardDg), &Nrows, &Ncols);
  readDfloatArray(fp, "IPDG overlapping patch diagonal scaling", &(oasDiagOpDg), &Nrows, &Ncols);
  readDfloatArray(fp, "IPDG overlapping patch backward matrix", &(oasBackDg), &Nrows, &Ncols);
  NpP = Nrows; //overlapping patch size


  /* SEMFEM data */
  readDfloatArray(fp, "SEMFEM r-coordinates", &(rFEM),&Nrows,&Ncols);
  readDfloatArray(fp, "SEMFEM s-coordinates", &(sFEM),&Nrows,&Ncols);
  readDfloatArray(fp, "SEMFEM t-coordinates", &(tFEM),&Nrows,&Ncols);
  NpFEM = Nrows;

  readIntArray   (fp, "SEMFEM reference mesh", &(FEMEToV), &Nrows, &Ncols);
  NelFEM = Nrows;

  readDfloatArray(fp, "SEMFEM interpolation matrix", &(SEMFEMInterp),&Nrows,&Ncols);


  fclose(fp);

  // find node indices of vertex nodes
  dfloat NODETOL = 1e-6;
  vertexNodes = (int*) calloc(Nverts, sizeof(int));
  for(int n=0;n<Np;++n){
    if( (r[n]+1)*(r[n]+1)+(s[n]+1)*(s[n]+1)+(t[n]+1)*(t[n]+1)<NODETOL)
      vertexNodes[0] = n;
    if( (r[n]-1)*(r[n]-1)+(s[n]+1)*(s[n]+1)+(t[n]+1)*(t[n]+1)<NODETOL)
      vertexNodes[1] = n;
    if( (r[n]+1)*(r[n]+1)+(s[n]-1)*(s[n]-1)+(t[n]+1)*(t[n]+1)<NODETOL)
      vertexNodes[2] = n;
    if( (r[n]+1)*(r[n]+1)+(s[n]+1)*(s[n]+1)+(t[n]-1)*(t[n]-1)<NODETOL)
      vertexNodes[3] = n;
  }
}


