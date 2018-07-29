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

#ifndef XXT_H
#define XXT_H

/* requires "types.h", and, when MPI is defined, "crystal.h" */
#if !defined(TYPES_H) || ( defined(MPI) && !defined(CRYSTAL_H) )
#warning "xxt.h" requires "types.h" and "crystal.h"
#endif

typedef struct xxt_data_ xxt_data;

#ifndef MPI
#  define crystal_data void
#endif

#define xxt_free xxt_jl_free
#define xxt_solve xxt_jl_solve
#define xxt_stats xxt_jl_stats

xxt_data *xxt_setup(uint n, const ulong *id,
                    uint nz, const uint *Ai, const uint *Aj, const real *A,
                    uint null_space, crystal_data *crystal);
void xxt_solve(real *x, xxt_data *data, const real *b);
void xxt_stats(xxt_data *data);
void xxt_free(xxt_data *data);

#ifndef MPI
#  undef crystal_data
#endif

#endif

