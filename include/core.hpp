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

#ifndef CORE_HPP
#define CORE_HPP

#include <mpi.h>
#include <occa.h>
#include <string>
#include "types.h"
#include "utils.hpp"
#include "settings.hpp"


void *occaHostMallocPinned(occa::device &device, size_t size, void *source, occa::memory &mem, occa::memory &h_mem);

// serial sort
void mysort(hlong *data, int N, const char *order);

// sort entries in an array in parallel
void parallelSort(int size, int rank, MPI_Comm comm,
      int N, void *vv, size_t sz,
      int (*compare)(const void *, const void *),
      void (*match)(void *, void *)
      );

void readDfloatArray(MPI_Comm comm, FILE *fp, const char *label, dfloat **A, int *Nrows, int* Ncols);
void readIntArray   (MPI_Comm comm, FILE *fp, const char *label, int **A   , int *Nrows, int* Ncols);

void matrixRightSolve(int NrowsA, int NcolsA, dfloat *A, int NrowsB, int NcolsB, dfloat *B, dfloat *C);
void matrixEig(int N, dfloat *A, dfloat *VR, dfloat *WR, dfloat *WI);
void matrixInverse(int N, dfloat *A);
dfloat matrixConditionNumber(int N, dfloat *A);

void matrixTranspose(const int M, const int N,
                     const dfloat  *A, const int LDA,
                           dfloat *AT, const int LDAT);

void matrixTranspose(const int M, const int N,
                     const dlong  *A, const int LDA,
                           dlong *AT, const int LDAT);

#endif