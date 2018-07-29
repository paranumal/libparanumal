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
#include <stdlib.h>
#include <stdarg.h>
#include <limits.h>
#include "errmem.h"
#include "types.h"
#include "sort.h"
#include "tuple_list.h"

const unsigned mi=2, ml=2, mr=2;

void test1()
{
  buffer buf;
  uint i,j;
  tuple_list tl;
  buffer_init(&buf,1024);
  tuple_list_init_max(&tl,mi,ml,mr,500);
  tl.n=tl.max;
  for(i=0;i<tl.n;++i) {
    int num1 = rand();
    slong num2 = rand();
    num2<<=(CHAR_BIT)*sizeof(int)-1;
    num2|=rand();
    num2<<=(CHAR_BIT)*sizeof(int)-1;
    num2|=rand();
    for(j=0;j<mi;++j) tl.vi[i*mi+j]=num1;
    for(j=0;j<ml;++j) tl.vl[i*ml+j]=num2;
    for(j=0;j<mr;++j) tl.vr[i*mr+j]=num1;
  }
  tuple_list_sort(&tl,0,&buf);
  for(i=0;i<tl.n;++i) {
    for(j=0;j<mi;++j) printf(" %016llx",(long long)tl.vi[i*mi+j]);
    for(j=0;j<ml;++j) printf(" %016llx",(long long)tl.vl[i*ml+j]);
    for(j=0;j<mr;++j) printf(" %g"     ,(double)   tl.vr[i*mr+j]);
    printf("\n");
  }
  printf("on the long:\n");
  tuple_list_sort(&tl,mi,&buf);
  for(i=0;i<tl.n;++i) {
    for(j=0;j<mi;++j) printf(" %016llx",(long long)tl.vi[i*mi+j]);
    for(j=0;j<ml;++j) printf(" %016llx",(long long)tl.vl[i*ml+j]);
    for(j=0;j<mr;++j) printf(" %g"     ,(double)   tl.vr[i*mr+j]);
    printf("\n");
  }
  buffer_free(&buf);
  tuple_list_free(&tl);
}

int main()
{
  test1();
  return 0;
}

