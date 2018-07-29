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
#include <string.h>
#include <limits.h>
#include "c99.h"
#include "name.h"
#include "fail.h"
#include "types.h"
#include "mem.h"
#include "sort.h"
#include "sarray_sort.h"

int main()
{
  struct rec { double d; slong l; sint i; float f; };
  buffer buf = {0,0,0};
  struct rec rec[500];
  uint i;
  
  for(i=0;i<500;++i) {
    sint num1 = rand() & 0xff;
    slong num2 = rand();
    num2<<=(CHAR_BIT)*sizeof(int)-1;
    num2|=rand();
    num2<<=(CHAR_BIT)*sizeof(int)-1;
    num2|=rand();
    num2= num2<0?-num2:num2;
    rec[i].d = num2;
    rec[i].f = num2;
    rec[i].l = num2;
    rec[i].i = num1;
  }
  sarray_sort_2(struct rec,rec,500, i,0, l,1, &buf);
  for(i=0;i<500;++i)
    printf("%g\t%g\t%ld\t%d\n",
      rec[i].d,rec[i].f,(long)rec[i].l,(int)rec[i].i);

  printf("\n");
  sarray_sort(struct rec,rec,500, l,1, &buf);
  for(i=0;i<500;++i)
    printf("%g\t%g\t%ld\t%d\n",
      rec[i].d,rec[i].f,(long)rec[i].l,(int)rec[i].i);
  buffer_free(&buf);
  return 0;
}

