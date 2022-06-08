/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"
#include "parAlmond/parAlmondAMGLevel.hpp"

namespace libp {

namespace parAlmond {

amgLevel::amgLevel(parCSR& _A, settings_t& _settings):
  multigridLevel(_A.Nrows, _A.Ncols, _A.platform, _settings, _A.comm) {

  A = _A;

  //determine smoother
  if (settings.compareSetting("PARALMOND SMOOTHER", "CHEBYSHEV")) {
    stype = CHEBYSHEV;
    settings.getSetting("PARALMOND CHEBYSHEV DEGREE", ChebyshevIterations);
  } else { //default to DAMPED_JACOBI
    stype = DAMPED_JACOBI;
  }
}

void amgLevel::Operator(deviceMemory<dfloat>& o_X, deviceMemory<dfloat>& o_Ax){
  A.SpMV(1.0, o_X, 0.0, o_Ax);
}

void amgLevel::coarsen   (deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_Rr){
  R.SpMV(1.0, o_r, 0.0, o_Rr);
}

void amgLevel::prolongate(deviceMemory<dfloat>& o_X, deviceMemory<dfloat>& o_Px){
  P.SpMV(1.0, o_X, 1.0, o_Px);
}

void amgLevel::residual  (deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_X,
                          deviceMemory<dfloat>& o_RES) {
  A.SpMV(-1.0, o_X, 1.0, o_RHS, o_RES);
}

void amgLevel::smooth(deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_X, bool x_is_zero){
  if(stype == DAMPED_JACOBI){
    A.smoothDampedJacobi(o_RHS, o_X, lambda,
                          x_is_zero, o_scratch);
  } else if(stype == CHEBYSHEV){
    A.smoothChebyshev(o_RHS, o_X, lambda0, lambda1,
                       x_is_zero, o_scratch,
                       ChebyshevIterations);
  }
}

void amgLevel::setupSmoother(){

  if (stype == DAMPED_JACOBI) {
    lambda = (4./3.)/A.rho;
  } else if (stype == CHEBYSHEV) {
    lambda1 = A.rho;
    lambda0 = A.rho/10.;
  }
}

void amgLevel::syncToDevice(){
  if (A.Nrows>0) A.syncToDevice();
  if (P.Nrows>0) P.syncToDevice();
  if (R.Nrows>0) R.syncToDevice();
}

void amgLevel::Report() {

  //This setup can be called by many subcommunicators, so only
  // print on the global root.
  int totalActive=(Nrows>0) ? 1:0;
  A.comm.Allreduce(totalActive);

  dlong minNrows=Nrows, maxNrows=Nrows;
  hlong totalNrows=Nrows;
  A.comm.Allreduce(maxNrows, Comm::Max);
  A.comm.Allreduce(totalNrows, Comm::Sum);
  dfloat avgNrows = (dfloat) totalNrows/totalActive;

  if (Nrows==0) minNrows=maxNrows; //set this so it's ignored for the global min
  A.comm.Allreduce(minNrows, Comm::Min);

  long long int nnz = A.diag.nnz+A.offd.nnz;
  long long int minNnz=nnz, maxNnz=nnz, totalNnz=nnz;
  A.comm.Allreduce(maxNnz, Comm::Max);
  A.comm.Allreduce(totalNnz, Comm::Sum);

  if (nnz==0) minNnz = maxNnz; //set this so it's ignored for the global min
  A.comm.Allreduce(minNnz, Comm::Min);

  dfloat nnzPerRow = (Nrows==0) ? 0 : (dfloat) nnz/Nrows;
  dfloat minNnzPerRow=nnzPerRow, maxNnzPerRow=nnzPerRow, avgNnzPerRow=nnzPerRow;
  A.comm.Allreduce(maxNnzPerRow, Comm::Max);
  A.comm.Allreduce(avgNnzPerRow, Comm::Sum);
  avgNnzPerRow /= totalActive;

  if (Nrows==0) minNnzPerRow = maxNnzPerRow;
  A.comm.Allreduce(minNnzPerRow, Comm::Min);

  char smootherString[BUFSIZ];
  if (stype==DAMPED_JACOBI)
    strcpy(smootherString, "Damped Jacobi   ");
  else if (stype==CHEBYSHEV)
    strcpy(smootherString, "Chebyshev       ");

  if (comm.rank()==0){
    printf(      "|  parAlmond |  %12lld  |  %12d  | %13d   |   %s|\n", (long long int) totalNrows, minNrows, (int)minNnzPerRow, smootherString);
    printf("      |            |                |  %12d  | %13d   |                   |\n", maxNrows, (int)maxNnzPerRow);
    printf("      |            |                |  %12d  | %13d   |                   |\n", (int)avgNrows, (int)avgNnzPerRow);
  }
}

} //namespace parAlmond

} //namespace libp
