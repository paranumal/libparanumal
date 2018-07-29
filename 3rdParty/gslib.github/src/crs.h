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

#ifndef CRS_H
#define CRS_H

#if !defined(COMM_H)
#warning "crs.h" requires "comm.h"
#endif

#define crs_setup PREFIXED_NAME(crs_setup)
#define crs_solve PREFIXED_NAME(crs_solve)
#define crs_stats PREFIXED_NAME(crs_stats)
#define crs_free  PREFIXED_NAME(crs_free )

struct crs_data;

struct crs_data *crs_setup(
  uint n, const ulong *id,
  uint nz, const uint *Ai, const uint *Aj, const double *A,
  uint null_space, const struct comm *comm);
void crs_solve(double *x, struct crs_data *data, double *b);
void crs_stats(struct crs_data *data);
void crs_free(struct crs_data *data);

#endif

