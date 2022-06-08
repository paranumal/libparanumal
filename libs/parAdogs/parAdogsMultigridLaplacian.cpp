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
#include "parAdogs/parAdogsMultigrid.hpp"

namespace libp {

namespace paradogs {

/*Create graph Laplacian from mesh data*/
void graph_t::CreateLaplacian() {

  Nlevels=1;
  L[0].A = parCSR(Nverts, Nverts, platform, comm);
  parCSR& A = L[0].A;

  A.rowOffsetL = VoffsetL;
  A.rowOffsetU = VoffsetU;
  A.colOffsetL = VoffsetL;
  A.colOffsetU = VoffsetU;

  /*Create a graph Laplacian from mesh info*/
  A.diag.rowStarts.malloc(Nverts+1);
  A.offd.rowStarts.malloc(Nverts+1);

  #pragma omp parallel for
  for (dlong n=0;n<Nverts+1;++n) {
    A.diag.rowStarts[n] = 0;
    A.offd.rowStarts[n] = 0;
  }

  for (dlong e=0;e<Nverts;++e) {
    A.diag.rowStarts[e+1]++; /*Count diagonal*/

    for (int n=0;n<Nfaces;++n) {
      const hlong gE = elements[e].E[n];
      if (gE!=-1) {
        if (gE>=VoffsetL && gE<VoffsetU) {
          A.diag.rowStarts[e+1]++; /*count connections per vert*/
        } else {
          A.offd.rowStarts[e+1]++; /*count connections per vert*/
        }
      }
    }
  }

  // count how many rows are shared
  A.offd.nzRows=0;
  for(dlong e=0;e<Nverts; e++) {
    if (A.offd.rowStarts[e+1]>0) A.offd.nzRows++;
  }

  A.offd.rows.malloc(A.offd.nzRows);
  A.offd.mRowStarts.malloc(A.offd.nzRows+1);

  /*cumulative sum*/
  dlong cnt=0;
  A.offd.mRowStarts[0]=0;
  for (dlong e=0;e<Nverts;++e) {
    if (A.offd.rowStarts[e+1]>0) {
      A.offd.rows[cnt] = e; //record row id
      A.offd.mRowStarts[cnt+1] = A.offd.mRowStarts[cnt] + A.offd.rowStarts[e+1];
      cnt++;
    }
    A.diag.rowStarts[e+1] += A.diag.rowStarts[e];
    A.offd.rowStarts[e+1] += A.offd.rowStarts[e];
  }
  A.diag.nnz = A.diag.rowStarts[Nverts];
  A.offd.nnz = A.offd.rowStarts[Nverts];

  /*Halo setup*/
  cnt=0;
  colIds.malloc(A.offd.nnz);
  for (dlong e=0;e<Nverts;++e) {
    for (int n=0;n<Nfaces;++n) {
      const hlong gE = elements[e].E[n];
      if (gE!=-1) {
        if (gE<VoffsetL || gE>=VoffsetU) {
          colIds[cnt++] = gE;
        }
      }
    }
  }
  A.haloSetup(colIds); //setup halo, and transform colIds to a local indexing
  Nhalo = A.Ncols-A.Nrows; /*Record how big the halo region is*/

  /*Build connectivity*/
  A.diagA.malloc(A.Ncols);
  A.diagInv.malloc(A.Nrows);
  A.diag.cols.malloc(A.diag.nnz);
  A.offd.cols.malloc(A.offd.nnz);
  A.diag.vals.malloc(A.diag.nnz);
  A.offd.vals.malloc(A.offd.nnz);

  A.diag.nnz=0;
  A.offd.nnz=0;
  for (dlong e=0;e<Nverts;++e) {
    A.diag.cols[A.diag.nnz] = e;
    pfloat& Ann = A.diag.vals[A.diag.nnz];
    A.diag.nnz++;

    Ann = 0.0;

    for (int n=0;n<Nfaces;++n) {
      const hlong gE = elements[e].E[n];
      if (gE!=-1) {
        if (gE>=VoffsetL && gE<VoffsetU) {
          A.diag.cols[A.diag.nnz] = static_cast<dlong>(gE-VoffsetL);
          A.diag.vals[A.diag.nnz] = -1.0;
          A.diag.nnz++;
        } else {
          A.offd.cols[A.offd.nnz] = colIds[A.offd.nnz];
          A.offd.vals[A.offd.nnz] = -1.0;
          A.offd.nnz++;
        }
        Ann += 1.0;
      }
    }
    A.diagA[e] = Ann;
    A.diagInv[e] = 1.0/Ann;
  }

  //fill the halo region
  A.halo.Exchange(A.diagA, 1);

  L[0].Nrows = A.Nrows;
  L[0].Ncols = A.Ncols;
  L[0].Nglobal = NVertsGlobal;

  /*Construct fine null vector*/
  L[0].null.malloc(Nverts);

  #pragma omp parallel for
  for (dlong n=0;n<Nverts;++n) {
    L[0].null[n] = 1.0/sqrt(NVertsGlobal);
  }
}

} //namespace paradogs

} //namespace libp
