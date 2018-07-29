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

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include "c99.h"
#include "name.h"
#include "fail.h"
#include "types.h"
#include "mem.h"
#include "sort.h"
#include "rdtsc.h"

#if 1

DEFINE_HW_COUNTER()

#define N (1<<20)

ulong A[N], out[N];
uint P[N];

int main()
{
  buffer buf = null_buffer;
  uint i;
  unsigned long long tic, toc;
  unsigned r;
  #define TIME(t, repeat, what) do { \
    for(r=repeat;r;--r) { what; } \
    tic = getticks(); \
    for(r=repeat;r;--r) { what; } \
    toc = getticks(); \
    t = toc-tic; \
  } while(0)

  for(i=0;i!=N;++i) {
    A[i]=rand();
    A[i]<<=CHAR_BIT*sizeof(int)-1;
    A[i]^=rand();
    A[i]<<=CHAR_BIT*sizeof(int)-1;
    A[i]^=rand();
    if(0) A[i]&=0x000ff00;
  }

  for(i=N;i;i>>=1) {
    unsigned long long t;
    TIME(t, (N/i), 
      sortv_long(out, A,i,sizeof(ulong), &buf));
    printf("sortv %d : %g cycles per item\n",
      (int)i, t/(double)(N/i)/(double)i);
  }

  for(i=N;i;i>>=1) {
    unsigned long long t;
    TIME(t, (N/i), 
      sortp_long(&buf,0, A,i,sizeof(ulong)));
    printf("sortp %d : %g cycles per item\n",
      (int)i, t/(double)(N/i)/(double)i);
  }

  buffer_free(&buf);
  return 0;
}

#else

int main()
{
  return 0;
}

#endif

