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
/* Refine Fiedler Vector               */
/****************************************/
void graph_t::Refine(const int level) {

  parCSR& A = L[level].A;
  memory<dfloat>& null = L[level].null;
  const dlong N = L[level].Nrows;
  const dlong Ncols = L[level].Ncols;

  memory<dfloat>& Fiedler = L[level].Fiedler;

  /*******************************************************/
  /*Improve fine Fiedler vector via Inverse Iteration    */
  /*******************************************************/

  const dfloat RELTOL = 3.0e-1;
  const dfloat CG_TOL = 1.0e-2;

  const int maxIters=1;

  memory<dfloat> x(Ncols);
  memory<dfloat> scratch(3*Ncols);
  memory<dfloat> AF = scratch;

  /*AF = A*F*/
  A.SpMV(1.0, Fiedler, 0.0, AF);

  /*theta = F^T * A * F */
  dfloat theta = 0.0;
  dfloat normAF = 0.0;
  for (dlong n=0;n<N;++n) {
    theta += Fiedler[n]*AF[n];
    normAF += AF[n]*AF[n];
  }
  comm.Allreduce(theta);
  comm.Allreduce(normAF);

  dfloat err = sqrt(std::abs(normAF - theta*theta))/theta;

  // if (rank==0) printf("Intial err = %f, theta = %f, ||AF|| = %f \n", err, theta, sqrt(normAF));

  for (int it=0;it<maxIters;++it) {

    if (err<RELTOL) break;

    #pragma omp parallel for
    for (dlong n=0;n<N;++n) {
      x[n] = Fiedler[n]/theta;
    }

    #pragma omp parallel for
    for (dlong n=0;n<N;++n) {
      Fiedler[n] = Fiedler[n] - AF[n]/theta;
    }

    /*Solve A_{l}*x = Fiedler*/
    (void) Solve(level, CG_TOL, Fiedler, x, scratch);
    // const int cg_iter = Solve(level, CG_TOL, Fiedler, x, scratch);

    /*Project out null vector*/
    dfloat dot=0.0;
    for (int n=0;n<N;++n) {
      dot += x[n]*null[n];
    }
    comm.Allreduce(dot);

    #pragma omp parallel for
    for (int n=0;n<N;++n) {
      x[n] -= dot*null[n];
    }

    dfloat normx = 0.0;
    for (dlong n=0;n<N;++n) {
      normx += x[n]*x[n];
    }
    comm.Allreduce(normx);
    normx = sqrt(normx);

    /*F = x /||x||*/
    #pragma omp parallel for
    for (dlong n=0;n<N;++n) {
      Fiedler[n] = x[n]/normx;
    }

    /*AF = A*F*/
    A.SpMV(1.0, Fiedler, 0.0, AF);

    /*theta = F^T * A * F */
    theta = 0.0;
    normAF = 0.0;
    for (dlong n=0;n<N;++n) {
      theta += Fiedler[n]*AF[n];
      normAF += AF[n]*AF[n];
    }
    comm.Allreduce(theta);
    comm.Allreduce(normAF);

    err = sqrt(std::abs(normAF - theta*theta))/theta;

    // if (rank==0)  printf("err = %f, theta = %f, ||AF|| = %f, cg_iter=%d\n", err, theta, sqrt(normAF), cg_iter);
  }
}

} //namespace paradogs

} //namespace libp
