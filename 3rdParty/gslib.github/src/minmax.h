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

#ifndef MINMAX_H
#define MINMAX_H

/* requires <math.h> and "types.h" */

#ifndef TYPES_H
#warning "minmax.h" depends on "types.h"
#endif

/*--------------------------------------------------------------------------
   Min, Max, Norm
  --------------------------------------------------------------------------*/

#define DECLMINMAX(type, prefix) \
static type prefix##min_2(type a, type b) { return b<a?b:a; } \
static type prefix##max_2(type a, type b) { return a>b?a:b; } \
static void prefix##minmax_2(type *min, type *max, type a, type b) \
{ if(b<a) *min=b, *max=a; else *min=a, *max=b; } \
static type prefix##min_3(type a, type b, type c) \
{ return b<a?(c<b?c:b):(c<a?c:a); } \
static type prefix##max_3(type a, type b, type c) \
{ return a>b?(a>c?a:c):(b>c?b:c); } \
static void prefix##minmax_3(type *min, type *max, type a, type b, type c) \
{ if(b<a)  *min=prefix##min_2(b,c), *max=prefix##max_2(a,c); \
  else    *min=prefix##min_2(a,c), *max=prefix##max_2(b,c); }

DECLMINMAX(int, i)
DECLMINMAX(unsigned, u)
DECLMINMAX(real, r)
#undef DECLMINMAX

static real r1norm_1(real a) { return fabsr(a); }
static real r1norm_2(real a, real b) { return fabsr(a)+fabsr(b); }
static real r1norm_3(real a, real b, real c)
{ return fabsr(a)+fabsr(b)+fabsr(c); }
static real r2norm_1(real a) { return sqrtr(a*a); }
static real r2norm_2(real a, real b) { return sqrtr(a*a+b*b); }
static real r2norm_3(real a, real b, real c)
{ return sqrtr(a*a+b*b+c*c); }

#endif

