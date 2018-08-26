/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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


//creators
csr * newCSRfromCOO(dlong N, hlong* globalRowStarts,
            		dlong NNZ,   hlong *Ai, hlong *Aj, dfloat *Avals);
void freeCSR(csr *A);
dcoo *newDCOO(parAlmond_t *parAlmond, csr *B);
hyb * newHYB(parAlmond_t *parAlmond, csr *csrA);


void axpy(csr *A, dfloat alpha, dfloat *x, dfloat beta, dfloat *y, bool nullSpace, dfloat nullSpacePenalty);

void axpy(parAlmond_t *parAlmond, dcoo *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void axpy(parAlmond_t *parAlmond, hyb *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y, bool nullSpace, dfloat nullSpacePenalty);

void axpy(parAlmond_t *parAlmond, ell *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void ax(parAlmond_t *parAlmond, coo *C, dfloat alpha, occa::memory o_x, occa::memory o_y);


//smoothing
void smoothJacobi      (parAlmond_t *parAlmond, agmgLevel *level, csr *A, dfloat *r, dfloat *x, bool x_is_zero);
void smoothDampedJacobi(parAlmond_t *parAlmond, agmgLevel *level, csr *A, dfloat *r, dfloat *x, bool x_is_zero);
void smoothChebyshev   (parAlmond_t *parAlmond, agmgLevel *level, csr *A, dfloat *r, dfloat *x, bool x_is_zero);
void smoothJacobi      (parAlmond_t *parAlmond, agmgLevel *level, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero);
void smoothDampedJacobi(parAlmond_t *parAlmond, agmgLevel *level, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero);
void smoothChebyshev   (parAlmond_t *parAlmond, agmgLevel *level, hyb *A, occa::memory o_r, occa::memory o_x, bool x_is_zero);

//halo exchange
void csrHaloSetup(csr *A, hlong *globalColStarts);
void csrHaloExchange(csr *A, size_t Nbytes, void *sourceBuffer, void *sendBuffer, void *recvBuffer);
void csrHaloExchangeStart(csr *A, size_t Nbytes, void *sourceBuffer, void *sendBuffer, void *recvBuffer);
void csrHaloExchangeFinish(csr *A);
void dcooHaloExchangeStart(dcoo *A, size_t Nbytes, void *sendBuffer, void *recvBuffer);
void dcooHaloExchangeFinish(dcoo *A);
void hybHaloExchangeStart(hyb *A, size_t Nbytes, void *sendBuffer, void *recvBuffer);
void hybHaloExchangeFinish(hyb *A);
