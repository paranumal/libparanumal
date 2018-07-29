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

#include <stdio.h>

#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "c99.h"
#include "types.h"
#include "name.h"
#include "fail.h"
#include "mem.h"
#include "obbox.h"
#include "poly.h"
#include "sort.h"
#include "sarray_sort.h"
#include "findpts_el.h"

struct uint_range { uint min, max; };
struct index_el { uint index, el; };

static struct dbl_range dbl_range_merge(struct dbl_range a, struct dbl_range b)
{
  struct dbl_range m;
  m.min = b.min<a.min?b.min:a.min,
  m.max = a.max>b.max?a.max:b.max;
  return m;
}

static sint ifloor(double x) { return floor(x); }
static sint iceil (double x) { return ceil (x); }

static uint hash_index_aux(double low, double fac, uint n, double x)
{
  const sint i = ifloor((x-low)*fac);
  return i<0 ? 0 : (n-1<(uint)i ? n-1 : (uint)i);
}

#define CODE_INTERNAL 0
#define CODE_BORDER 1
#define CODE_NOT_FOUND 2

#define D 2
#define WHEN_3D(a)
#include "findpts_local_imp.h"
#undef WHEN_3D
#undef D

#define D 3
#define WHEN_3D(a) a
#include "findpts_local_imp.h"
#undef WHEN_3D
#undef D
