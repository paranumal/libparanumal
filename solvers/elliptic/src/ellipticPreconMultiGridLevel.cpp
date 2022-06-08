/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

void MGLevel::Operator(deviceMemory<dfloat>& o_X, deviceMemory<dfloat>& o_Ax) {
  elliptic.Operator(o_X,o_Ax);
}

void MGLevel::residual(deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_X, deviceMemory<dfloat>& o_RES) {
  elliptic.Operator(o_X,o_RES);

  // subtract res = rhs - A*x
  platform.linAlg().axpy(elliptic.Ndofs, 1.f, o_RHS, -1.f, o_RES);
}

void MGLevel::coarsen(deviceMemory<dfloat>& o_X, deviceMemory<dfloat>& o_Rx) {

  linAlg_t& linAlg = platform.linAlg();

  if (elliptic.disc_c0) {
    //scratch spaces
    deviceMemory<dfloat>& o_wx = o_smootherResidual;
    deviceMemory<dfloat>& o_RxL = o_transferScratch;

    //pre-weight
    linAlg.amxpy(elliptic.Ndofs, 1.0, elliptic.o_weightG, o_X, 0.0, o_wx);

    elliptic.gHalo.ExchangeStart(o_wx, 1);

    if(mesh.NlocalGatherElements/2)
      partialCoarsenKernel(mesh.NlocalGatherElements/2,
                           mesh.o_localGatherElementList,
                           elliptic.o_GlobalToLocal,
                           o_P, o_wx, o_RxL);

    elliptic.gHalo.ExchangeFinish(o_wx, 1);

    if(mesh.NglobalGatherElements)
      partialCoarsenKernel(mesh.NglobalGatherElements,
                           mesh.o_globalGatherElementList,
                           elliptic.o_GlobalToLocal,
                           o_P, o_wx, o_RxL);

    ellipticC.ogsMasked.GatherStart(o_Rx, o_RxL, 1, ogs::Add, ogs::Trans);

    if((mesh.NlocalGatherElements+1)/2)
      partialCoarsenKernel((mesh.NlocalGatherElements+1)/2,
                           mesh.o_localGatherElementList + mesh.NlocalGatherElements/2,
                           elliptic.o_GlobalToLocal,
                           o_P, o_wx, o_RxL);

    ellipticC.ogsMasked.GatherFinish(o_Rx, o_RxL, 1, ogs::Add, ogs::Trans);

  } else {
    coarsenKernel(mesh.Nelements, o_P, o_X, o_Rx);
  }
}

void MGLevel::prolongate(deviceMemory<dfloat>& o_X, deviceMemory<dfloat>& o_Px) {

  linAlg_t& linAlg = platform.linAlg();

  if (elliptic.disc_c0) {
    //scratch spaces
    deviceMemory<dfloat>& o_PxG = o_smootherResidual;
    deviceMemory<dfloat>& o_PxL = o_transferScratch;

    ellipticC.gHalo.ExchangeStart(o_X, 1);

    if(meshC.NlocalGatherElements/2)
      partialProlongateKernel(meshC.NlocalGatherElements/2,
                              meshC.o_localGatherElementList,
                              ellipticC.o_GlobalToLocal,
                              o_P, o_X, o_PxL);

    ellipticC.gHalo.ExchangeFinish(o_X, 1);

    if(meshC.NglobalGatherElements)
      partialProlongateKernel(meshC.NglobalGatherElements,
                              meshC.o_globalGatherElementList,
                              ellipticC.o_GlobalToLocal,
                              o_P, o_X, o_PxL);

    //ogs_notrans -> no summation at repeated nodes, just one value
    elliptic.ogsMasked.GatherStart(o_PxG, o_PxL, 1, ogs::Add, ogs::NoTrans);

    if((meshC.NlocalGatherElements+1)/2)
      partialProlongateKernel((meshC.NlocalGatherElements+1)/2,
                              meshC.o_localGatherElementList + meshC.NlocalGatherElements/2,
                              ellipticC.o_GlobalToLocal,
                              o_P, o_X, o_PxL);

    elliptic.ogsMasked.GatherFinish(o_PxG, o_PxL, 1, ogs::Add, ogs::NoTrans);

    linAlg.axpy(elliptic.Ndofs, 1.f, o_PxG, 1.f, o_Px);

  } else {
    prolongateKernel(mesh.Nelements, o_P, o_X, o_Px);
  }
}

void MGLevel::smooth(deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_X, bool x_is_zero) {
  if (stype==JACOBI) {
    smoothJacobi(o_RHS, o_X, x_is_zero);
  } else if (stype==CHEBYSHEV) {
    smoothChebyshev(o_RHS, o_X, x_is_zero);
  }
}

void MGLevel::smoothJacobi(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_X, bool xIsZero) {

  linAlg_t& linAlg = platform.linAlg();

  deviceMemory<dfloat>& o_RES = o_smootherResidual;

  if (xIsZero) {
    linAlg.amxpy(elliptic.Ndofs, 1.0, o_invDiagA, o_r, 0.0, o_X);
    return;
  }

  //res = r-Ax
  Operator(o_X,o_RES);
  linAlg.axpy(elliptic.Ndofs, 1.f, o_r, -1.f, o_RES);

  //smooth the fine problem x = x + S(r-Ax)
  linAlg.amxpy(elliptic.Ndofs, 1.0, o_invDiagA, o_RES, 1.0, o_X);
}

void MGLevel::smoothChebyshev (deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_X, bool xIsZero) {

  const dfloat theta = 0.5*(lambda1+lambda0);
  const dfloat delta = 0.5*(lambda1-lambda0);
  const dfloat invTheta = 1.0/theta;
  const dfloat sigma = theta/delta;
  dfloat rho_n = 1./sigma;
  dfloat rho_np1;

  deviceMemory<dfloat>& o_RES = o_smootherResidual;
  deviceMemory<dfloat>& o_Ad  = o_smootherResidual2;
  deviceMemory<dfloat>& o_d   = o_smootherUpdate;

  linAlg_t& linAlg = platform.linAlg();

  if(xIsZero){ //skip the Ax if x is zero
    //res = S*r
    linAlg.amxpy(elliptic.Ndofs, 1.0, o_invDiagA, o_r, 0.f, o_RES);

    //d = invTheta*res
    linAlg.axpy(elliptic.Ndofs, invTheta, o_RES, 0.f, o_d);
  } else {
    //res = S*(r-Ax)
    Operator(o_X,o_RES);
    linAlg.axpy(elliptic.Ndofs, 1.f, o_r, -1.f, o_RES);
    linAlg.amx(elliptic.Ndofs, 1.f, o_invDiagA, o_RES);

    //d = invTheta*res
    linAlg.axpy(elliptic.Ndofs, invTheta, o_RES, 0.f, o_d);
  }

  for (int k=0;k<ChebyshevIterations;k++) {
    //x_k+1 = x_k + d_k
    if (xIsZero&&(k==0))
      linAlg.axpy(elliptic.Ndofs, 1.f, o_d, 0.f, o_X);
    else
      linAlg.axpy(elliptic.Ndofs, 1.f, o_d, 1.f, o_X);

    //r_k+1 = r_k - SAd_k
    Operator(o_d,o_Ad);
    linAlg.amxpy(elliptic.Ndofs, -1.f, o_invDiagA, o_Ad, 1.f, o_RES);

    rho_np1 = 1.0/(2.*sigma-rho_n);
    dfloat rhoDivDelta = 2.0*rho_np1/delta;

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    linAlg.axpy(elliptic.Ndofs, rhoDivDelta, o_RES, rho_np1*rho_n, o_d);

    rho_n = rho_np1;
  }
  //x_k+1 = x_k + d_k
  linAlg.axpy(elliptic.Ndofs, 1.f, o_d, 1.0, o_X);
}


/******************************************
*
* MG Level Setup
*
*******************************************/

dlong  MGLevel::NsmootherResidual=0;
dlong  MGLevel::Nscratch=0;
memory<dfloat> MGLevel::smootherResidual;
deviceMemory<dfloat> MGLevel::o_smootherResidual;
deviceMemory<dfloat> MGLevel::o_smootherResidual2;
deviceMemory<dfloat> MGLevel::o_smootherUpdate;
deviceMemory<dfloat> MGLevel::o_transferScratch;

//build a level and connect it to the next one
MGLevel::MGLevel(elliptic_t& _elliptic,
                 dlong _Nrows, dlong _Ncols,
                 int Nc, int NpCoarse):
  multigridLevel(_Nrows, _Ncols,
                 _elliptic.platform,
                 _elliptic.settings,
                 _elliptic.comm),
  elliptic(_elliptic),
  mesh(_elliptic.mesh) {

  SetupSmoother();
  AllocateStorage();

  if (   mesh.elementType==Mesh::QUADRILATERALS
      || mesh.elementType==Mesh::HEXAHEDRA) {
    mesh.DegreeRaiseMatrix1D(Nc, mesh.N, P);
  } else if (mesh.elementType==Mesh::TRIANGLES) {
    mesh.DegreeRaiseMatrixTri2D(Nc, mesh.N, P);
  } else { //Mesh::TETRAHEDRA
    mesh.DegreeRaiseMatrixTet3D(Nc, mesh.N, P);
  }
  o_P = elliptic.platform.malloc<dfloat>(P);

  //build kernels
  properties_t kernelInfo = elliptic.platform.props();

  // set kernel name suffix
  std::string suffix;
  if(mesh.elementType==Mesh::TRIANGLES)
    suffix = "Tri2D";
  else if(mesh.elementType==Mesh::QUADRILATERALS)
    suffix = "Quad2D";
  else if(mesh.elementType==Mesh::TETRAHEDRA)
    suffix = "Tet3D";
  else if(mesh.elementType==Mesh::HEXAHEDRA)
    suffix = "Hex3D";

  std::string oklFilePrefix = DELLIPTIC "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  kernelInfo["defines/" "p_NqFine"]= mesh.N+1;
  kernelInfo["defines/" "p_NqCoarse"]= Nc+1;

  kernelInfo["defines/" "p_NpFine"]= mesh.Np;
  kernelInfo["defines/" "p_NpCoarse"]= NpCoarse;

  int blockMax = 256;
  if (elliptic.platform.device.mode() == "CUDA") blockMax = 512;

  int NblockVFine = std::max(1,blockMax/mesh.Np);
  int NblockVCoarse = std::max(1,blockMax/NpCoarse);
  kernelInfo["defines/" "p_NblockVFine"]= NblockVFine;
  kernelInfo["defines/" "p_NblockVCoarse"]= NblockVCoarse;

  if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    fileName   = oklFilePrefix + "ellipticPreconCoarsen" + suffix + oklFileSuffix;
    kernelName = "ellipticPartialPreconCoarsen" + suffix;
    partialCoarsenKernel = elliptic.platform.buildKernel(fileName, kernelName, kernelInfo);

    fileName   = oklFilePrefix + "ellipticPreconProlongate" + suffix + oklFileSuffix;
    kernelName = "ellipticPartialPreconProlongate" + suffix;
    partialProlongateKernel = elliptic.platform.buildKernel(fileName, kernelName, kernelInfo);
  } else { //IPDG
    fileName   = oklFilePrefix + "ellipticPreconCoarsen" + suffix + oklFileSuffix;
    kernelName = "ellipticPreconCoarsen" + suffix;
    coarsenKernel = elliptic.platform.buildKernel(fileName, kernelName, kernelInfo);

    fileName   = oklFilePrefix + "ellipticPreconProlongate" + suffix + oklFileSuffix;
    kernelName = "ellipticPreconProlongate" + suffix;
    prolongateKernel = elliptic.platform.buildKernel(fileName, kernelName, kernelInfo);
  }
}

void MGLevel::AllocateStorage() {
  // extra storage for smoothing op
  if (NsmootherResidual < Ncols) {
    smootherResidual.malloc(Ncols, 0);
    o_smootherResidual  = elliptic.platform.malloc<dfloat>(smootherResidual);
    o_smootherResidual2 = elliptic.platform.malloc<dfloat>(smootherResidual);
    o_smootherUpdate    = elliptic.platform.malloc<dfloat>(smootherResidual);
    NsmootherResidual = Ncols;
  }

  if (Nscratch < mesh.Nelements*mesh.Np) {
    memory<dfloat> dummy(mesh.Nelements*mesh.Np,0);
    o_transferScratch = elliptic.platform.malloc<dfloat>(dummy);
    Nscratch = mesh.Nelements*mesh.Np;
  }
}

void MGLevel::Report() {

  int totalActive=(Nrows>0) ? 1:0;
  mesh.comm.Allreduce(totalActive);

  dlong minNrows=Nrows, maxNrows=Nrows;
  hlong totalNrows=Nrows;
  mesh.comm.Allreduce(maxNrows, Comm::Max);
  mesh.comm.Allreduce(totalNrows, Comm::Sum);
  dfloat avgNrows = static_cast<dfloat>(totalNrows)/totalActive;

  if (Nrows==0) Nrows=maxNrows; //set this so it's ignored for the global min
  mesh.comm.Allreduce(minNrows, Comm::Min);

  char smootherString[BUFSIZ];
  if (stype==JACOBI)
    strcpy(smootherString, "Damped Jacobi   ");
  else if (stype==CHEBYSHEV)
    strcpy(smootherString, "Chebyshev       ");

  //This setup can be called by many subcommunicators, so only
  // print on the global root.
  if (mesh.rank==0){
    printf(      "|    pMG     |    %10lld  |    %10d  |   Matrix-free   |   %s|\n", (long long int)totalNrows, minNrows, smootherString);
    printf("      |            |                |    %10d  |     Degree %2d   |                   |\n", maxNrows, mesh.N);
    printf("      |            |                |    %10d  |                 |                   |\n", (int) avgNrows);
  }
}

void MGLevel::SetupSmoother() {

  //set up the fine problem smoothing
  memory<dfloat> diagA   (Nrows);
  memory<dfloat> invDiagA(Nrows);
  elliptic.BuildOperatorDiagonal(diagA);

  for (dlong n=0;n<Nrows;n++) {
    invDiagA[n] = 1.0/diagA[n];
  }

  o_invDiagA = elliptic.platform.malloc<dfloat>(invDiagA);

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

    for (dlong n=0;n<Nrows;n++)
      invDiagA[n] *= lambda0;

    //update diagonal with weight
    o_invDiagA.copyFrom(invDiagA);
  }
}


//------------------------------------------------------------------------
//
//  Estimate max Eigenvalue of diagA^{-1}*A
//
//------------------------------------------------------------------------

dfloat MGLevel::maxEigSmoothAx(){

  const dlong N = Nrows;
  const dlong M = Ncols;

  linAlg_t& linAlg = platform.linAlg();

  int k = 10;

  hlong Ntotal = Nrows;
  mesh.comm.Allreduce(Ntotal);
  if(k > Ntotal) k = static_cast<int>(Ntotal);

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  memory<double> H(k*k,0.0);

  // allocate memory for basis
  memory<dfloat> Vx(M);
  memory<deviceMemory<dfloat>> o_V(k+1);

  deviceMemory<dfloat> o_Vx  = elliptic.platform.malloc<dfloat>(Vx);
  deviceMemory<dfloat> o_AVx = elliptic.platform.malloc<dfloat>(Vx);

  for(int i=0; i<=k; i++)
    o_V[i] = elliptic.platform.malloc<dfloat>(Vx);

  // generate a random vector for initial basis vector
  for (dlong i=0;i<N;i++) Vx[i] = (dfloat) drand48();

  o_Vx.copyFrom(Vx); //copy to device

  dfloat norm_vo =  linAlg.norm2(N, o_Vx, mesh.comm);

  linAlg.axpy(N, 1./norm_vo, o_Vx, 0.f, o_V[0]);

  for(int j=0; j<k; j++){
    // v[j+1] = invD*(A*v[j])
    Operator(o_V[j],o_AVx);
    linAlg.amxpy(N, 1.0, o_invDiagA, o_AVx, 0.0, o_V[j+1]);

    // modified Gram-Schmidth
    for(int i=0; i<=j; i++){
      // H(i,j) = v[i]'*A*v[j]
      dfloat hij =  linAlg.innerProd(N, o_V[i], o_V[j+1], mesh.comm);

      // v[j+1] = v[j+1] - hij*v[i]
      linAlg.axpy(N, -hij, o_V[i], 1.f, o_V[j+1]);

      H[i + j*k] = static_cast<double>(hij);
    }

    if(j+1 < k){
      // v[j+1] = v[j+1]/||v[j+1]||
      dfloat norm_vj =  linAlg.norm2(N, o_V[j+1], mesh.comm);
      linAlg.scale(N, 1.0/norm_vj, o_V[j+1]);

      H[j+1+ j*k] = static_cast<double>(norm_vj);
    }
  }

  memory<double> WR(k);
  memory<double> WI(k);

  linAlg_t::matrixEigenValues(k, H, WR, WI);

  double rho = 0.;

  for(int i=0; i<k; i++){
    double rho_i  = sqrt(WR[i]*WR[i] + WI[i]*WI[i]);

    if(rho < rho_i) {
      rho = rho_i;
    }
  }

  // if((mesh.rank==0)) printf("weight = %g \n", rho);

  return rho;
}
