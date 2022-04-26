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

#include "parAdogs.hpp"
#include "parAdogs/parAdogsGraph.hpp"
#include "parAdogs/parAdogsPartition.hpp"

namespace libp {

namespace paradogs {

/****************************************/
/* Solve A_{l}*x = b                    */
/****************************************/
int graph_t::Solve(const int level, 
                   const dfloat TOL,
                   memory<dfloat>& r,
                   memory<dfloat>& x,
                   memory<dfloat>& scratch) {

  parCSR& A = L[level].A;
  const dlong N = A.Nrows;
  const dlong Ncols = L[level].Ncols;

  memory<dfloat> p  = scratch + 0*Ncols;
  memory<dfloat> Ap = scratch + 1*Ncols;
  memory<dfloat> z  = scratch + 2*Ncols;

  // register scalars
  dfloat rdotz1 = 0.0;
  dfloat rdotz2 = 0.0;
  dfloat alpha = 0.0, beta = 0.0, pAp = 0.0;
  dfloat rdotr = 1.0;
  const int MAXIT = 5000;

  /* We assume that x is initialized to some guess and
     r = b-A*x */

  /*Compute x = A^{-1} b*/
  int cg_iter;
  for(cg_iter=0;cg_iter<MAXIT;++cg_iter){

    // Exit if tolerance is reached, taking at least one step.
    if (((cg_iter == 0) && (rdotr == 0.0)) ||
        ((cg_iter > 0) && (sqrt(rdotr) <= TOL))) {
      break;
    }

    // z = Precon^{-1} r
    MultigridVcycle(level, r, z);

    // r.z
    rdotz2 = rdotz1;
    rdotz1 = 0.0;
    for (dlong n=0;n<N;++n) {
      rdotz1 += z[n]*r[n];
    }
    comm.Allreduce(rdotz1);

    beta = (cg_iter==0) ? 0.0 : rdotz1/rdotz2;

    // p = z + beta*p
    if (cg_iter==0) {
      #pragma omp parallel for
      for (dlong n=0;n<N;++n) {
        p[n] = z[n];
      }
    } else {
      #pragma omp parallel for
      for (dlong n=0;n<N;++n) {
        p[n] = z[n] + beta*p[n];
      }
    }

    // A*p
    A.SpMV(1.0, p, 0.0, Ap);

    // p.Ap
    pAp = 0.0;
    for (dlong n=0;n<N;++n) {
      pAp += p[n]*Ap[n];
    }
    comm.Allreduce(pAp);

    alpha = rdotz1/pAp;

    //  x <= x + alpha*p
    //  r <= r - alpha*A*p
    //  dot(r,r)
    rdotr = 0.0;
    for (dlong n=0;n<N;++n) {
      x[n] = x[n] + alpha*p[n];
      r[n] = r[n] - alpha*Ap[n];
      rdotr += r[n]*r[n];
    }
    comm.Allreduce(rdotr);

    if(rdotr<0) printf("WARNING CG: rdotr = %17.15lf\n", rdotr);

    // printf("CG: it %d, r norm %12.12le, alpha = %le \n", cg_iter+1, sqrt(rdotr), alpha);
  }

  return cg_iter;
}

} //namespace paradogs

} //namespace libp
