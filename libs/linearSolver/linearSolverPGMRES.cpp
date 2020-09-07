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

#include "linearSolver.hpp"

#define PGMRES_RESTART 20

pgmres::pgmres(dlong _N, dlong _Nhalo,
         platform_t& _platform, settings_t& _settings,
         int _weighted, occa::memory& _o_weight):
  linearSolver_t(_N, _Nhalo, _platform, _settings) {

  // Make sure LinAlg has the necessary kernels
  platform.linAlg.InitKernels({"axpy", "zaxpy",
                               "innerProd", "weightedInnerProd",
                               "norm2", "weightedNorm2"});

  dlong Ntotal = N + Nhalo;

  //Number of iterations between restarts
  //TODO make this modifyable via settings
  restart=PGMRES_RESTART;

  dfloat *dummy = (dfloat *) calloc(Ntotal,sizeof(dfloat)); //need this to avoid uninitialized memory warnings

  o_V = new occa::memory[restart];
  for(int i=0; i<restart; ++i){
    o_V[i] = platform.malloc(Ntotal*sizeof(dfloat), dummy);
  }

  H  = (dfloat *) calloc((restart+1)*(restart+1), sizeof(dfloat));
  sn = (dfloat *) calloc(restart, sizeof(dfloat));
  cs = (dfloat *) calloc(restart, sizeof(dfloat));
  s = (dfloat *) calloc(restart+1, sizeof(dfloat));
  y = (dfloat *) calloc(restart, sizeof(dfloat));

  /*aux variables */
  o_Ax = platform.malloc(Ntotal*sizeof(dfloat), dummy);
  o_z  = platform.malloc(Ntotal*sizeof(dfloat), dummy);
  o_r  = platform.malloc(Ntotal*sizeof(dfloat), dummy);
  free(dummy);

  weighted = _weighted;
  o_w = _o_weight;
}

int pgmres::Solve(solver_t& solver, precon_t& precon,
               occa::memory &o_x, occa::memory &o_b,
               const dfloat tol, const int MAXIT, const int verbose) {

  int rank = platform.rank;
  linAlg_t &linAlg = platform.linAlg;

  // compute A*x
  solver.Operator(o_x, o_Ax);

  // subtract z = b - A*x
  linAlg.zaxpy(N, -1.f, o_Ax, 1.f, o_b, o_z);

  // r = Precon^{-1} (r-A*x)
  precon.Operator(o_z, o_r);

  dfloat nr;
  if (weighted)
    nr = linAlg.weightedNorm2(N, o_w, o_r, platform.comm);
  else
    nr = linAlg.norm2(N, o_r, platform.comm);

  dfloat error = nr;
  const dfloat TOL = mymax(tol*nr,tol);

  if (verbose&&(rank==0))
    printf("PGMRES: initial res norm %12.12f \n", nr);

  int iter=0;

  //exit if tolerance is reached
  if(error<=TOL) return iter;


  for(iter=1;iter<MAXIT;){

    s[0] = nr;

    // V(:,0) = r/nr
    linAlg.axpy(N, (1./nr), o_r, 0., o_V[0]);

    //Construct orthonormal basis via Gram-Schmidt
    for(int i=0;i<restart;++i){
      // compute z = A*V(:,i)
      solver.Operator(o_V[i], o_z);

      // r = Precon^{-1} z
      precon.Operator(o_z, o_r);

      for(int k=0; k<=i; ++k){
        dfloat hki;
        if (weighted)
          hki = linAlg.weightedInnerProd(N, o_w, o_r, o_V[k], platform.comm);
        else
          hki = linAlg.innerProd(N, o_r, o_V[k], platform.comm);

        // r = r - hki*V[k]
        linAlg.axpy(N, -hki, o_V[k], 1.0, o_r);

        // H(k,i) = hki
        H[k + i*(restart+1)] = hki;
      }

      dfloat nw;
      if (weighted)
        nw = linAlg.weightedNorm2(N, o_w, o_r, platform.comm);
      else
        nw = linAlg.norm2(N, o_r, platform.comm);
      H[i+1 + i*(restart+1)] = nw;

      // V(:,i+1) = r/nw
      if (i<restart-1)
        linAlg.axpy(N, (1./nw), o_r, 0., o_V[i+1]);

      //apply Givens rotation
      for(int k=0; k<i; ++k){
        const dfloat h1 = H[k +     i*(restart+1)];
        const dfloat h2 = H[k + 1 + i*(restart+1)];

        H[k +     i*(restart+1)] =  cs[k]*h1 + sn[k]*h2;
        H[k + 1 + i*(restart+1)] = -sn[k]*h1 + cs[k]*h2;
      }

      // form i-th rotation matrix
      const dfloat h1 = H[i+    i*(restart+1)];
      const dfloat h2 = H[i+1 + i*(restart+1)];
      const dfloat hr = sqrt(h1*h1 + h2*h2);
      cs[i] = h1/hr;
      sn[i] = h2/hr;

      H[i   + i*(restart+1)] = cs[i]*h1 + sn[i]*h2;
      H[i+1 + i*(restart+1)] = 0;

      //approximate residual norm
      s[i+1] = -sn[i]*s[i];
      s[i]   =  cs[i]*s[i];

      iter++;
      error = fabs(s[i+1]);

      if (verbose&&(rank==0)) {
        printf("GMRES: it %d, approx residual norm %12.12le \n", iter, error);
      }

      if(error < TOL || iter==MAXIT) {
        //update approximation
        UpdateGMRES(o_x, i+1);
        break;
      }
    }

    //exit if tolerance is reached
    if(error < TOL || iter==MAXIT) break;

    //update approximation
    UpdateGMRES(o_x, restart);

    // compute A*x
    solver.Operator(o_x, o_Ax);

    // subtract z = b - A*x
    linAlg.zaxpy(N, -1.f, o_Ax, 1.f, o_b, o_z);

    // r = Precon^{-1} (r-A*x)
    precon.Operator(o_z, o_r);

    if (weighted)
      nr = linAlg.weightedNorm2(N, o_w, o_r, platform.comm);
    else
      nr = linAlg.norm2(N, o_r, platform.comm);

    error = nr;
    //exit if tolerance is reached
    if(error<=TOL) return iter;
  }

  return iter;
}

void pgmres::UpdateGMRES(occa::memory& o_x, const int I){

  for(int k=I-1; k>=0; --k){
    y[k] = s[k];

    for(int m=k+1; m<I; ++m)
      y[k] -= H[k + m*(restart+1)]*y[m];

    y[k] /= H[k + k*(restart+1)];
  }

  //TODO this is really a GEMM, should write it that way
  for(int j=0; j<I; ++j){
    platform.linAlg.axpy(N, y[j], o_V[j], 1.0, o_x);
  }
}

pgmres::~pgmres() {
  if(H) free(H);
  if(sn) free(sn);
  if(cs) free(cs);
  if(s) free(s);
  if(y) free(y);

  if (o_V) delete[] o_V;
}