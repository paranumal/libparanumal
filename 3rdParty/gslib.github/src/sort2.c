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

#include <limits.h>
#include <string.h>

#include "types.h"

typedef uint Index;

#define sort jl_sort
#define Value uint
#define Data sort_data
typedef struct { Value v; Index i; } Data;
#include "sort_imp.c"

#undef sort
#undef Value
#undef Data

#ifdef GLOBAL_INT
#  define Value ulong
#  define Data sort_data_long
   typedef struct { Value v; Index i; } Data;
#  define radix_count         radix_count_long
#  define radix_offsets       radix_offsets_long
#  define radix_zeros         radix_zeros_long
#  define radix_pass          radix_pass_long
#  define radix_sort          radix_sort_long
#  define radix_index_pass_b  radix_index_pass_b_long
#  define radix_index_pass_m  radix_index_pass_m_long
#  define radix_index_pass_e  radix_index_pass_e_long
#  define radix_index_pass_be radix_index_pass_be_long
#  define radix_index_sort    radix_index_sort_long
#  define merge_sort          merge_sort_long
#  define merge_index_sort    merge_index_sort_long
#  define sort                sort_long
#  define index_sort          index_sort_long
#  include "sort_imp.c"
#endif

