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

#ifndef SPARSE_CHOLESKY_H
#define SPARSE_CHOLESKY_H

#if !defined(TYPES_H) || !defined(MEM_H)
#warning "sparse_cholesky.h" requires "types.h" and "mem.h"
#endif

#define sparse_cholesky_factor PREFIXED_NAME(sparse_cholesky_factor)
#define sparse_cholesky_solve  PREFIXED_NAME(sparse_cholesky_solve )
#define sparse_cholesky_free   PREFIXED_NAME(sparse_cholesky_free  )

struct sparse_cholesky {
  uint n, *Lrp, *Lj;
  double *L, *D;
};

/* input data is the usual CSR
   matrix is n by n
   Arp has n+1 elements
   elements of row i are A [Arp[i]], ..., A [Arp[i+1]-1]
              in columns Aj[Arp[i]], ..., Aj[Arp[i+1]-1]
*/
void sparse_cholesky_factor(uint n, const uint *Arp, const uint *Aj,
                            const double *A,
                            struct sparse_cholesky *out, buffer *buf);
                            
/* x = A^(-1) b;  works when x and b alias */
void sparse_cholesky_solve(
  double *x, const struct sparse_cholesky *fac, double *b);

void sparse_cholesky_free(struct sparse_cholesky *fac);

#endif

