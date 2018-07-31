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



dfloat innerProd(dlong n, dfloat *a, dfloat *b);

void doubleInnerProd(dlong n, dfloat *aDotbc, dfloat *a, dfloat *b, dfloat *c);

void kcycleCombinedOp1(dlong n, dfloat *aDotbc, dfloat *a, dfloat *b, dfloat *c, dfloat *w, bool weighted);

void kcycleCombinedOp2(dlong n, dfloat *aDotbcd, dfloat *a, dfloat *b, dfloat *c, dfloat* d, dfloat *w, bool weighted);

void vectorAdd(dlong n, dfloat alpha, dfloat *x, dfloat beta, dfloat *y);

dfloat vectorAddInnerProd(dlong n, dfloat alpha, dfloat *x, dfloat beta, dfloat *y, dfloat *w, bool weighted);

void dotStar(dlong m, dfloat *a, dfloat *b);

void scaleVector(dlong m, dfloat *a, dfloat alpha);

void setVector(dlong m, dfloat *a, dfloat alpha);

dfloat sumVector(dlong m, dfloat *a);

void addScalar(dlong m, dfloat alpha, dfloat *a);

void randomize(dlong m, dfloat *a);

dfloat maxEntry(dlong n, dfloat *a);

void scaleVector(parAlmond_t *parAlmond, dlong N, occa::memory o_a, dfloat alpha);

void setVector(parAlmond_t *parAlmond, dlong N, occa::memory o_a, dfloat alpha);

dfloat sumVector(parAlmond_t *parAlmond, dlong N, occa::memory o_a);

void addScalar(parAlmond_t *parAlmond, dlong N, dfloat alpha, occa::memory o_a);

void dotStar(parAlmond_t *parAlmond, dlong N, occa::memory o_a, occa::memory o_b);

void dotStar(parAlmond_t *parAlmond, dlong N, dfloat alpha, occa::memory o_a,
	           occa::memory o_b, dfloat beta, occa::memory o_c);

dfloat innerProd(parAlmond_t *parAlmond, dlong N, occa::memory o_x, occa::memory o_y);

// returns aDotbc[0] = a\dot b, aDotbc[1] = a\dot c, aDotbc[2] = b\dot b,
void kcycleCombinedOp1(parAlmond_t *parAlmond, dlong n, dfloat *aDotbc, occa::memory o_a,
                                        occa::memory o_b, occa::memory o_c, occa::memory o_w, bool weighted);

// returns aDotbcd[0] = a\dot b, aDotbcd[1] = a\dot c, aDotbcd[2] = a\dot d,
void kcycleCombinedOp2(parAlmond_t *parAlmond, dlong n, dfloat *aDotbcd, occa::memory o_a,
                                              occa::memory o_b, occa::memory o_c, occa::memory o_d,
                                              occa::memory o_w, bool weighted);

// y = beta*y + alpha*x, and return y\dot y
dfloat vectorAddInnerProd(parAlmond_t *parAlmond, dlong n, dfloat alpha, occa::memory o_x,
                                                          dfloat beta, occa::memory o_y,
                                                          occa::memory o_w, bool weighted);

void vectorAdd(parAlmond_t *parAlmond, dlong N, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y);

void vectorAdd(parAlmond_t *parAlmond, dlong N, dfloat alpha, occa::memory o_x,
	 dfloat beta, occa::memory o_y, occa::memory o_z);
