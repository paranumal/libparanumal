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

extern "C" {
  void dsyev_ (char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO);
}

namespace libp {

namespace paradogs {

/****************************************/
/* Serial Inertial Bipartition          */
/****************************************/
void graph_t::InertialBipartition(const dfloat targetFraction[2]) {

  memory<int> partition(Nverts);

  memory<double> I;
  I.calloc(9);

  memory<dfloat> x, y, z;

  if (dim==2) {
    x.malloc(Nverts);
    y.malloc(Nverts);

    /*Compute center of mass of each element*/
    for (dlong e=0;e<Nverts;++e) {
      x[e]=0.0;
      y[e]=0.0;
      for (int v=0;v<NelementVerts;++v) {
        x[e] += elements[e].EX[v];
        y[e] += elements[e].EY[v];
      }
      x[e] /= NelementVerts;
      y[e] /= NelementVerts;
    }

    /*Compute center of mass of whole mesh*/
    memory<double> avg(2);

    avg[0]=0.0;
    avg[1]=0.0;
    for (dlong e=0;e<Nverts;++e) {
      avg[0] += x[e];
      avg[1] += y[e];
    }
    comm.Allreduce(avg);

    avg[0] /= NVertsGlobal;
    avg[1] /= NVertsGlobal;

    for (dlong e=0;e<Nverts;++e) {
      const dfloat X = x[e] - avg[0];
      const dfloat Y = y[e] - avg[1];

      I[0] += X*X; I[1] += X*Y;
      I[2] += Y*X; I[3] += Y*Y;
    }
    comm.Allreduce(I);

  } else {
    x.malloc(Nverts);
    y.malloc(Nverts);
    z.malloc(Nverts);

    /*Compute center of mass of each element*/
    for (dlong e=0;e<Nverts;++e) {
      x[e]=0.0;
      y[e]=0.0;
      z[e]=0.0;
      for (int v=0;v<NelementVerts;++v) {
        x[e] += elements[e].EX[v];
        y[e] += elements[e].EY[v];
        z[e] += elements[e].EZ[v];
      }
      x[e] /= NelementVerts;
      y[e] /= NelementVerts;
      z[e] /= NelementVerts;
    }

    /*Compute center of mass of whole mesh*/
    memory<double> avg(3);

    avg[0]=0.0;
    avg[1]=0.0;
    avg[2]=0.0;
    for (dlong e=0;e<Nverts;++e) {
      avg[0] += x[e];
      avg[1] += y[e];
      avg[2] += z[e];
    }
    comm.Allreduce(avg);

    avg[0] /= NVertsGlobal;
    avg[1] /= NVertsGlobal;
    avg[2] /= NVertsGlobal;

    for (dlong e=0;e<Nverts;++e) {
      const dfloat X = x[e] - avg[0];
      const dfloat Y = y[e] - avg[1];
      const dfloat Z = z[e] - avg[2];

      I[0] += X*X; I[1] += X*Y; I[2] += X*Z;
      I[3] += Y*X; I[4] += Y*Y; I[5] += Y*Z;
      I[6] += Z*X; I[7] += Z*Y; I[8] += Z*Z;
    }
    comm.Allreduce(I);
  }

  /*Find the principal axis of inertia*/
  int N = dim;
  int INFO = -999;
  char JOBZ='V';
  char UPLO='L';
  int LDA = N;
  double W[3];
  int LWORK = 8;
  double WORK[8];
  dsyev_(&JOBZ, &UPLO, &N, I.ptr(), &LDA, W, WORK, &LWORK, &INFO);

  LIBP_ABORT("Paradogs: dsyev_ reports info = " << INFO << " in InertialBipartition",
             INFO);

  /*Find the largest eigenvalue*/
  double max = W[0];
  int maxloc = 0;
  for (int i=1;i<dim;++i) {
    if (W[i]>max) {
      max = W[i];
      maxloc = i;
    }
  }
  // printf("max = %f, maxloc = %d \n", max, maxloc);

  /*Princial axis is the eigenvector with largest eigenvalue*/
  double a[3];
  memory<double> maxV = I + maxloc*N;
  for (int i=0;i<N;++i) {
    a[i] = maxV[i];
  }

  /*Use principal axis to bipartion graph*/
  memory<dfloat> F(Nverts);

  if (dim==2) {
    for (dlong e=0;e<Nverts;++e) {
      F[e] = x[e]*a[0] + y[e]*a[1];
    }
  } else {
    for (dlong e=0;e<Nverts;++e) {
      F[e] = x[e]*a[0] + y[e]*a[1] + z[e]*a[2];
    }
  }

  const hlong K = std::ceil(targetFraction[0]*NVertsGlobal);
  const dfloat pivot = ParallelPivot(Nverts, F, K, comm);

  for (dlong n=0;n<Nverts;++n) {
    if (F[n]<=pivot) {
      partition[n] = 0;
    } else {
      partition[n] = 1;
    }
  }

  /*Split the graph according to this partitioning*/
  Split(partition);
}

} //namespace paradogs

} //namespace libp
