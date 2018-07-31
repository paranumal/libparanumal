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


void agmgAx        (void **args, dfloat *x, dfloat *Ax);
void agmgCoarsen   (void **args, dfloat *r, dfloat *Rr);
void agmgProlongate(void **args, dfloat *x, dfloat *Px);
void agmgSmooth    (void **args, dfloat *rhs, dfloat *x, bool x_is_zero);

void device_agmgAx        (void **args, occa::memory &o_x, occa::memory &o_Ax);
void device_agmgCoarsen   (void **args, occa::memory &o_r, occa::memory &o_Rr);
void device_agmgProlongate(void **args, occa::memory &o_x, occa::memory &o_Px);
void device_agmgSmooth    (void **args, occa::memory &o_r, occa::memory &o_x, bool x_is_zero);

void setupSmoother(parAlmond_t *parAlmond, agmgLevel *level, SmoothType s);
void setupExactSolve(parAlmond_t *parAlmond, agmgLevel *level, bool nullSpace, dfloat nullSpacePenalty);
void exactCoarseSolve(parAlmond_t *parAlmond, int N, dfloat *rhs, dfloat *x);
void device_exactCoarseSolve(parAlmond_t *parAlmond, int N, occa::memory o_rhs, occa::memory o_x);
