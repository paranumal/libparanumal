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

#include "elliptic.hpp"
#include "ellipticPrecon.hpp"

void MGLevel::Ax(occa::memory &o_X, occa::memory &o_Ax) {
  elliptic.Operator(o_X,o_Ax);
}

void MGLevel::residual(occa::memory &o_RHS, occa::memory &o_X, occa::memory &o_RES) {
  elliptic.Operator(o_X,o_RES);

  // subtract res = rhs - A*x
  linAlg.axpy(mesh.Np*mesh.Nelements, 1.f, o_RHS, -1.f, o_RES);
}

void MGLevel::coarsen(occa::memory &o_X, occa::memory &o_Rx) {
  if (elliptic.disc_c0)
    linAlg.amx(mesh.Nelements*NpF, 1.0, o_weightF, o_X);

  coarsenKernel(mesh.Nelements, o_R, o_X, o_Rx);

  if (elliptic.disc_c0) {
    elliptic.ogsMasked->GatherScatter(o_Rx, ogs_dfloat, ogs_add, ogs_sym);
    if (elliptic.Nmasked)
      elliptic.maskKernel(elliptic.Nmasked, elliptic.o_maskIds, o_Rx);
  }
}

void MGLevel::prolongate(occa::memory &o_X, occa::memory &o_Px) {
  prolongateKernel(mesh.Nelements, o_R, o_X, o_Px);
}

void MGLevel::smooth(occa::memory &o_RHS, occa::memory &o_X, bool x_is_zero) {
  if (stype==JACOBI) {
    smoothJacobi(o_RHS, o_X, x_is_zero);
  } else if (stype==CHEBYSHEV) {
    smoothChebyshev(o_RHS, o_X, x_is_zero);
  }
}

void MGLevel::smoothJacobi(occa::memory &o_r, occa::memory &o_X, bool xIsZero) {

  dlong Ntotal = mesh.Np*mesh.Nelements;
  occa::memory &o_RES = o_smootherResidual;

  if (xIsZero) {
    linAlg.amxpy(Ntotal, 1.0, o_invDiagA, o_r, 0.0, o_X);
    return;
  }

  //res = r-Ax
  Ax(o_X,o_RES);
  linAlg.axpy(Ntotal, 1.f, o_r, -1.f, o_RES);

  //smooth the fine problem x = x + S(r-Ax)
  linAlg.amxpy(Ntotal, 1.0, o_invDiagA, o_RES, 1.0, o_X);
}

void MGLevel::smoothChebyshev (occa::memory &o_r, occa::memory &o_X, bool xIsZero) {

  const dfloat theta = 0.5*(lambda1+lambda0);
  const dfloat delta = 0.5*(lambda1-lambda0);
  const dfloat invTheta = 1.0/theta;
  const dfloat sigma = theta/delta;
  dfloat rho_n = 1./sigma;
  dfloat rho_np1;

  dlong Ntotal = mesh.Np*mesh.Nelements;
  occa::memory &o_RES = o_smootherResidual;
  occa::memory &o_Ad  = o_smootherResidual2;
  occa::memory &o_d   = o_smootherUpdate;

  if(xIsZero){ //skip the Ax if x is zero
    //res = S*r
    linAlg.amxpy(Ntotal, 1.0, o_invDiagA, o_r, 0.f, o_RES);

    //d = invTheta*res
    linAlg.axpy(Ntotal, invTheta, o_RES, 0.f, o_d);
  } else {
    //res = S*(r-Ax)
    Ax(o_X,o_RES);
    linAlg.axpy(Ntotal, 1.f, o_r, -1.f, o_RES);
    linAlg.amx(Ntotal, 1.f, o_invDiagA, o_RES);

    //d = invTheta*res
    linAlg.axpy(Ntotal, invTheta, o_RES, 0.f, o_d);
  }

  for (int k=0;k<ChebyshevIterations;k++) {
    //x_k+1 = x_k + d_k
    if (xIsZero&&(k==0))
      linAlg.axpy(Ntotal, 1.f, o_d, 0.f, o_X);
    else
      linAlg.axpy(Ntotal, 1.f, o_d, 1.f, o_X);

    //r_k+1 = r_k - SAd_k
    Ax(o_d,o_Ad);
    linAlg.amxpy(Ntotal, -1.f, o_invDiagA, o_Ad, 1.f, o_RES);

    rho_np1 = 1.0/(2.*sigma-rho_n);
    dfloat rhoDivDelta = 2.0*rho_np1/delta;

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    linAlg.axpy(Ntotal, rhoDivDelta, o_RES, rho_np1*rho_n, o_d);

    rho_n = rho_np1;
  }
  //x_k+1 = x_k + d_k
  linAlg.axpy(Ntotal, 1.f, o_d, 1.0, o_X);
}


/******************************************
*
* MG Level Setup
*
*******************************************/

size_t  MGLevel::smootherResidualBytes;
dfloat* MGLevel::smootherResidual;
occa::memory MGLevel::o_smootherResidual;
occa::memory MGLevel::o_smootherResidual2;
occa::memory MGLevel::o_smootherUpdate;

//build a level and connect it to the previous one
MGLevel::MGLevel(elliptic_t& _elliptic, int k,
                 int Nf, int Npf, occa::memory o_weightF_,
                 parAlmond::KrylovType ktype_, parAlmond::CycleType ctype):
  multigridLevel(_elliptic.mesh.Nelements*_elliptic.mesh.Np,
                (_elliptic.mesh.Nelements+_elliptic.mesh.totalHaloPairs)*_elliptic.mesh.Np,
                ktype_,
                _elliptic.mesh.comm),
  elliptic(_elliptic),
  mesh(_elliptic.mesh),
  linAlg(_elliptic.linAlg) {

  weighted = false;

  printf("Creating level order : Coarse %d and Fine:%d \n", elliptic.mesh.N, Nf);

  //save size and weight of previous level
  NpF = Npf;
  o_weightF = o_weightF_;

  //check for weighted inner products
  if (elliptic.settings.compareSetting("DISCRETIZATION","CONTINUOUS")) {
    weighted = true;
    o_weight = elliptic.o_weight;
  }

  SetupSmoother();
  AllocateStorage(k, ctype);

  if (mesh.N<Nf) {
    mesh.BuildBasisCoarsen(&R, o_R, Nf, mesh.N);

    //build kernels
    occa::properties kernelInfo = elliptic.props;

    // set kernel name suffix
    char *suffix;
    if(mesh.elementType==TRIANGLES)
      suffix = strdup("Tri2D");
    else if(mesh.elementType==QUADRILATERALS)
      suffix = strdup("Quad2D");
    else if(mesh.elementType==TETRAHEDRA)
      suffix = strdup("Tet3D");
    else if(mesh.elementType==HEXAHEDRA)
      suffix = strdup("Hex3D");

    char fileName[BUFSIZ], kernelName[BUFSIZ];

    kernelInfo["defines/" "p_NqFine"]= Nf+1;
    kernelInfo["defines/" "p_NqCoarse"]= mesh.N+1;

    kernelInfo["defines/" "p_NpFine"]= Npf;
    kernelInfo["defines/" "p_NpCoarse"]= mesh.Np;

    int NblockVFine = 512/Npf;
    int NblockVCoarse = 512/mesh.Np;
    kernelInfo["defines/" "p_NblockVFine"]= NblockVFine;
    kernelInfo["defines/" "p_NblockVCoarse"]= NblockVCoarse;

    sprintf(fileName, DELLIPTIC "/okl/ellipticPreconCoarsen%s.okl", suffix);
    sprintf(kernelName, "ellipticPreconCoarsen%s", suffix);
    coarsenKernel = buildKernel(mesh.device, fileName, kernelName, kernelInfo, mesh.comm);

    sprintf(fileName, DELLIPTIC "/okl/ellipticPreconProlongate%s.okl", suffix);
    sprintf(kernelName, "ellipticPreconProlongate%s", suffix);
    prolongateKernel = buildKernel(mesh.device, fileName, kernelName, kernelInfo, mesh.comm);
  }
}

void MGLevel::AllocateStorage(int k, parAlmond::CycleType ctype) {
  // extra storage for smoothing op
  size_t Nbytes = Ncols*sizeof(dfloat);
  if (smootherResidualBytes < Nbytes) {
    if (o_smootherResidual.size()) {
      free(smootherResidual);
      o_smootherResidual.free();
      o_smootherResidual2.free();
      o_smootherUpdate.free();
    }

    smootherResidual = (dfloat *) calloc(Ncols,sizeof(dfloat));
    o_smootherResidual = mesh.device.malloc(Nbytes,smootherResidual);
    o_smootherResidual2 = mesh.device.malloc(Nbytes,smootherResidual);
    o_smootherUpdate = mesh.device.malloc(Nbytes,smootherResidual);
    smootherResidualBytes = Nbytes;
  }

  if (k) x     = (dfloat *) calloc(Ncols,sizeof(dfloat));
  if (k) rhs   = (dfloat *) calloc(Nrows,sizeof(dfloat));
  if (k) o_x   = mesh.device.malloc(Ncols*sizeof(dfloat),x);
  if (k) o_rhs = mesh.device.malloc(Nrows*sizeof(dfloat),rhs);
 

  res  = (dfloat *) calloc(Ncols,sizeof(dfloat));
  o_res = mesh.device.malloc(Ncols*sizeof(dfloat),res);

  //kcycle vectors
  if (ctype==parAlmond::KCYCLE) {
    if ((k>0) && (k<NUMKCYCLES+1)) {
      ck = (dfloat *) calloc(Ncols,sizeof(dfloat));
      vk = (dfloat *) calloc(Nrows,sizeof(dfloat));
      wk = (dfloat *) calloc(Nrows,sizeof(dfloat));
      o_ck = mesh.device.malloc(Ncols*sizeof(dfloat),ck);
      o_vk = mesh.device.malloc(Nrows*sizeof(dfloat),vk);
      o_wk = mesh.device.malloc(Nrows*sizeof(dfloat),wk);
    }
  }
}

void MGLevel::Report() {

  hlong hNrows = (hlong) Nrows;

  dlong minNrows=0, maxNrows=0;
  hlong totalNrows=0;
  dfloat avgNrows;

  MPI_Allreduce(&Nrows, &maxNrows, 1, MPI_DLONG, MPI_MAX, mesh.comm);
  MPI_Allreduce(&hNrows, &totalNrows, 1, MPI_HLONG, MPI_SUM, mesh.comm);
  avgNrows = (dfloat) totalNrows/mesh.size;

  if (Nrows==0) Nrows=maxNrows; //set this so it's ignored for the global min
  MPI_Allreduce(&Nrows, &minNrows, 1, MPI_DLONG, MPI_MIN, mesh.comm);

  char smootherString[BUFSIZ];
  if (stype==JACOBI)
    strcpy(smootherString, "Damped Jacobi   ");
  else if (stype==CHEBYSHEV)
    strcpy(smootherString, "Chebyshev       ");

  if (mesh.rank==0){
    printf(     "|    pMG     |    %10d  |   Matrix-free   |   %s|\n",minNrows, smootherString);
    printf("     |            |    %10d  |     Degree %2d   |                   |\n", maxNrows, mesh.N);
    printf("     |            |    %10d  |                 |                   |\n", (int) avgNrows);
  }
}


void MGLevel::SetupSmoother() {

  //set up the fine problem smoothing
  dfloat *diagA    = (dfloat*) calloc(mesh.Np*mesh.Nelements, sizeof(dfloat));
  dfloat *invDiagA = (dfloat*) calloc(mesh.Np*mesh.Nelements, sizeof(dfloat));
  elliptic.BuildOperatorDiagonal(diagA);

  for (dlong n=0;n<mesh.Nelements*mesh.Np;n++)
    invDiagA[n] = 1.0/diagA[n];

  o_invDiagA = mesh.device.malloc(mesh.Np*mesh.Nelements*sizeof(dfloat), invDiagA);

  if (elliptic.settings.compareSetting("MULTIGRID SMOOTHER","CHEBYSHEV")) {
    stype = CHEBYSHEV;

    ChebyshevIterations = 2; //default to degree 2
    elliptic.settings.getSetting("MULTIGRID CHEBYSHEV DEGREE", ChebyshevIterations);

    //estimate the max eigenvalue of S*A
    dfloat rho = maxEigSmoothAx();

    lambda1 = rho;
    lambda0 = rho/10.;
  } else {
    stype = JACOBI;

    //estimate the max eigenvalue of S*A
    dfloat rho = maxEigSmoothAx();

    //set the stabilty weight (jacobi-type interation)
    lambda0 = (4./3.)/rho;

    for (dlong n=0;n<mesh.Np*mesh.Nelements;n++)
      invDiagA[n] *= lambda0;

    //update diagonal with weight
    o_invDiagA.copyFrom(invDiagA);
  }
  free(diagA);
  free(invDiagA);
}

extern "C"
{
  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI,
              double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );
}

static void eig(const int Nrows, double *A, double *WR, double *WI){

  int NB  = 256;
  char JOBVL  = 'V';
  char JOBVR  = 'V';
  int     N = Nrows;
  int   LDA = Nrows;
  int  LWORK  = (NB+2)*N;

  double *WORK  = new double[LWORK];
  double *VL  = new double[Nrows*Nrows];
  double *VR  = new double[Nrows*Nrows];

  int INFO = -999;

  dgeev_ (&JOBVL, &JOBVR, &N, A, &LDA, WR, WI,
    VL, &LDA, VR, &LDA, WORK, &LWORK, &INFO);

  if(INFO) {
    stringstream ss;
    ss << "MGLevel: dgeev reports info = " << INFO;
    LIBP_WARNING(ss.str())
  }

  delete [] VL;
  delete [] VR;
  delete [] WORK;
}

dfloat MGLevel::maxEigSmoothAx(){

  const dlong N = Nrows;
  const dlong M = Ncols;

  int k = 10;

  hlong Nlocal = (hlong) Nrows;
  hlong Ntotal = 0;
  MPI_Allreduce(&Nlocal, &Ntotal, 1, MPI_HLONG, MPI_SUM, mesh.comm);
  if(k > Ntotal) k = (int) Ntotal;

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  double *H = (double *) calloc(k*k,sizeof(double));

  // allocate memory for basis
  dfloat *Vx = (dfloat*) calloc(M, sizeof(dfloat));
  //  occa::memory *o_V = (occa::memory *) calloc(k+1, sizeof(occa::memory));
  occa::memory *o_V = new occa::memory[k+1];

  occa::memory o_Vx  = mesh.device.malloc(M*sizeof(dfloat),Vx);
  occa::memory o_AVx = mesh.device.malloc(M*sizeof(dfloat),Vx);

  for(int i=0; i<=k; i++)
    o_V[i] = mesh.device.malloc(M*sizeof(dfloat),Vx);

  // generate a random vector for initial basis vector
  for (dlong i=0;i<N;i++) Vx[i] = (dfloat) drand48();

  //gather-scatter
  if (weighted) {
    elliptic.ogsMasked->GatherScatter(Vx, ogs_dfloat, ogs_add, ogs_sym);
    for (dlong i=0;i<elliptic.Nmasked;i++) Vx[elliptic.maskIds[i]] = 0.;
  }

  o_Vx.copyFrom(Vx); //copy to device

  dfloat norm_vo;
  if (weighted)
    norm_vo =  linAlg.weightedNorm2(N, o_weight, o_Vx, mesh.comm);
  else
    norm_vo =  linAlg.norm2(N, o_Vx, mesh.comm);

  linAlg.axpy(N, 1./norm_vo, o_Vx, 0.f, o_V[0]);

  for(int j=0; j<k; j++){
    // v[j+1] = invD*(A*v[j])
    Ax(o_V[j],o_AVx);
    linAlg.amxpy(N, 1.0, o_invDiagA, o_AVx, 0.0, o_V[j+1]);

    // modified Gram-Schmidth
    for(int i=0; i<=j; i++){
      // H(i,j) = v[i]'*A*v[j]
      dfloat hij;
      if (weighted)
        hij =  linAlg.weightedInnerProd(N, o_weight, o_V[i], o_V[j+1], mesh.comm);
      else
        hij =  linAlg.innerProd(N, o_V[i], o_V[j+1], mesh.comm);

      // v[j+1] = v[j+1] - hij*v[i]
      linAlg.axpy(N, -hij, o_V[i], 1.f, o_V[j+1]);

      H[i + j*k] = (double) hij;
    }

    if(j+1 < k){
      // v[j+1] = v[j+1]/||v[j+1]||
      dfloat norm_vj;
      if (weighted)
        norm_vj =  linAlg.weightedNorm2(N, o_weight, o_V[j+1], mesh.comm);
      else
        norm_vj =  linAlg.norm2(N, o_V[j+1], mesh.comm);
      linAlg.scale(N, 1.0/norm_vj, o_V[j+1]);

      H[j+1+ j*k] = (double) norm_vj;
    }
  }

  double *WR = (double *) calloc(k,sizeof(double));
  double *WI = (double *) calloc(k,sizeof(double));

  eig(k, H, WR, WI);

  double rho = 0.;

  for(int i=0; i<k; i++){
    double rho_i  = sqrt(WR[i]*WR[i] + WI[i]*WI[i]);

    if(rho < rho_i) {
      rho = rho_i;
    }
  }

  // free memory
  free(H);
  free(WR);
  free(WI);

  free(Vx);
  o_Vx.free();
  o_AVx.free();
  for(int i=0; i<=k; i++) o_V[i].free();
  delete[] o_V;

  // if((mesh.rank==0)&&(mesh.settings.compareSetting("VERBOSE","TRUE"))) printf("weight = %g \n", rho);

  return rho;
}
