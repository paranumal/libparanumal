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

#include "linearSolver.hpp"

namespace libp {

namespace LinearSolver {

#define PGMRES_RESTART 20

template<typename T> pgmres<T>::pgmres(dlong _N, dlong _Nhalo,
         platform_t& _platform, settings_t& _settings, comm_t _comm):
  linearSolverBase_t(_N, _Nhalo, _platform, _settings, _comm) {

  // Make sure LinAlg has the necessary kernels
  platform.linAlg().InitKernels({"axpy", "zaxpy",
	"innerProd", "norm2"});

  //Number of iterations between restarts
  //TODO make this modifyable via settings
  restart=PGMRES_RESTART;

  H .malloc((restart+1)*(restart+1), 0.0);
  sn.malloc(restart);
  cs.malloc(restart);
  s.malloc(restart+1);
  y.malloc(restart);
}

template<typename T> int pgmres<T>::Solve(operator_t& linearOperator, operator_t& precon,
               deviceMemory<T>& o_x, deviceMemory<T>& o_b,
               const T tol, const int MAXIT, const int verbose) {

  int rank = comm.rank();
  linAlg_t &linAlg = platform.linAlg();

  /*Pre-reserve memory pool space to avoid some unnecessary re-sizing*/
  dlong Ntotal = N + Nhalo;
  platform.reserve<T>((restart+3)*Ntotal
                          +(restart+3) * platform.memPoolAlignment<T>());

  /*aux variables */
  deviceMemory<T> o_Ax = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_z  = platform.reserve<T>(Ntotal);
  deviceMemory<T> o_r  = platform.reserve<T>(Ntotal);

  platform.reserve<pfloat>(2*Ntotal +
                           + 4 * platform.memPoolAlignment<T>());
  
  deviceMemory<pfloat> o_pfloat_r  = platform.reserve<pfloat>(Ntotal);
  deviceMemory<pfloat> o_pfloat_z  = platform.reserve<pfloat>(Ntotal);

  
  memory<deviceMemory<T>> o_V(restart);
  for(int i=0; i<restart; ++i){
    o_V[i] = platform.reserve<T>(Ntotal);
  }

  // compute A*x
  linearOperator.Operator(o_x, o_Ax);

  // subtract z = b - A*x
  linAlg.zaxpy(N, (T)-1.f, o_Ax, (T)1.f, o_b, o_z);

  // r = Precon^{-1} (r-A*x)
  //  precon.Operator(o_z, o_r);
  if(sizeof(pfloat)==sizeof(T)){
    precon.Operator(o_z, o_r);
  }
  else{
    linAlg.d2p(N, o_z, o_pfloat_z);
    precon.Operator(o_pfloat_z, o_pfloat_r);
    linAlg.p2d(N, o_pfloat_r, o_r);
  }


  T nr = linAlg.norm2(N, o_r, comm);

  T error = nr;
  const T TOL = std::max(tol*nr,tol);

  if (verbose&&(rank==0))
    printf("PGMRES: initial res norm %12.12f \n", nr);

  int iter=0;

  //exit if tolerance is reached
  if(error<=TOL) return iter;


  for(iter=1;iter<MAXIT;){

    s[0] = nr;

    // V(:,0) = r/nr
    linAlg.axpy(N, (T)(1./nr), o_r, (T) 0., o_V[0]);

    //Construct orthonormal basis via Gram-Schmidt
    for(int i=0;i<restart;++i){
      // compute z = A*V(:,i)
      linearOperator.Operator(o_V[i], o_z);

      // r = Precon^{-1} z
      //      precon.Operator(o_z, o_r);
      if(sizeof(pfloat)==sizeof(T)){
	precon.Operator(o_z, o_r);
      }
      else{
	linAlg.d2p(N, o_z, o_pfloat_z);
	precon.Operator(o_pfloat_z, o_pfloat_r);
	linAlg.p2d(N, o_pfloat_r, o_r);
      }

      for(int k=0; k<=i; ++k){
        T hki = linAlg.innerProd(N, o_r, o_V[k], comm);

        // r = r - hki*V[k]
        linAlg.axpy(N, -hki, o_V[k], (T)1.0, o_r);

        // H(k,i) = hki
        H[k + i*(restart+1)] = hki;
      }

      T nw = linAlg.norm2(N, o_r, comm);
      H[i+1 + i*(restart+1)] = nw;

      // V(:,i+1) = r/nw
      if (i<restart-1)
        linAlg.axpy(N, (T)(1./nw), o_r, (T)0., o_V[i+1]);

      //apply Givens rotation
      for(int k=0; k<i; ++k){
        const T h1 = H[k +     i*(restart+1)];
        const T h2 = H[k + 1 + i*(restart+1)];

        H[k +     i*(restart+1)] =  cs[k]*h1 + sn[k]*h2;
        H[k + 1 + i*(restart+1)] = -sn[k]*h1 + cs[k]*h2;
      }

      // form i-th rotation matrix
      const T h1 = H[i+    i*(restart+1)];
      const T h2 = H[i+1 + i*(restart+1)];
      const T hr = sqrt(h1*h1 + h2*h2);
      cs[i] = h1/hr;
      sn[i] = h2/hr;

      H[i   + i*(restart+1)] = cs[i]*h1 + sn[i]*h2;
      H[i+1 + i*(restart+1)] = 0;

      //approximate residual norm
      s[i+1] = -sn[i]*s[i];
      s[i]   =  cs[i]*s[i];

      iter++;
      error = std::abs(s[i+1]);

      if (verbose&&(rank==0)) {
        printf("GMRES: it %d, approx residual norm %12.12le \n", iter, error);
      }

      if(error < TOL || iter==MAXIT) {
        //update approximation
        UpdateGMRES(o_V, o_x, i+1);
        break;
      }
    }

    //exit if tolerance is reached
    if(error < TOL || iter==MAXIT) break;

    //update approximation
    UpdateGMRES(o_V, o_x, restart);

    // compute A*x
    linearOperator.Operator(o_x, o_Ax);

    // subtract z = b - A*x
    linAlg.zaxpy(N, (T)-1.f, o_Ax, (T)1.f, o_b, o_z);

    // r = Precon^{-1} (r-A*x)    
    //    precon.Operator(o_z, o_r);
    // double check direction
    if(sizeof(pfloat)==sizeof(T)){
      precon.Operator(o_z, o_r);
    }
    else{
      linAlg.d2p(N, o_z, o_pfloat_z);
      precon.Operator(o_pfloat_z, o_pfloat_r);
      linAlg.p2d(N, o_pfloat_r, o_r);
    }
    

    nr = linAlg.norm2(N, o_r, comm);

    error = nr;
    //exit if tolerance is reached
    if(error<=TOL) return iter;
  }

  return iter;
}

template<typename T> void pgmres<T>::UpdateGMRES(memory<deviceMemory<T>>& o_V,
                         deviceMemory<T>& o_x,
                         const int I){

  for(int k=I-1; k>=0; --k){
    y[k] = s[k];

    for(int m=k+1; m<I; ++m)
      y[k] -= H[k + m*(restart+1)]*y[m];

    y[k] /= H[k + k*(restart+1)];
  }

  //TODO this is really a GEMM, should write it that way
  for(int j=0; j<I; ++j){
    platform.linAlg().axpy(N, y[j], o_V[j], (T)1.0, o_x);
  }
}

} //namespace LinearSolver

} //namespace libp
