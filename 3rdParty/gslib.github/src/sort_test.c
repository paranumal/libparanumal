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

#define SMALL 22
#define NUM   500
#define SI 9

ulong A[NUM][SI], Av[NUM];
uint  B[NUM][SI], Bv[NUM];

uint P[NUM], Q[NUM];

int main()
{
  buffer buf = {0,0,0};
  uint i;

  /*buffer_init(&buf, sortp_long_worksize(NUM,0));*/

#if 0
  printf("\nsource:\n");
#endif
  for(i=0;i!=NUM;++i) {
    A[i][0]=rand();
    A[i][0]<<=CHAR_BIT*sizeof(int)-1;
    A[i][0]^=rand();
    A[i][0]<<=CHAR_BIT*sizeof(int)-1;
    A[i][0]^=rand();
    if(0) A[i][0]&=0x000ff00;
    B[i][0]=A[i][0];
#if 0    
    printf("%016lx\t%016lx\n",(unsigned long)A[i][0],(unsigned long)B[i][0]);
#endif
  }
#if 0
  printf("\n");
#endif
  printf("merge sort:\n");
  for(i=0;i!=SMALL;++i) Q[i]=SMALL-1-i;
  sortv_long(Av,  &A[0][0],SMALL,sizeof(ulong[SI]), &buf);
  sortp_long(&buf,0, &A[0][0],SMALL,sizeof(ulong[SI]));
    memcpy(P,buf.ptr,SMALL*sizeof(uint));
  memcpy(buf.ptr,Q,SMALL*sizeof(uint));
  sortp_long(&buf,1, &A[0][0],SMALL,sizeof(ulong[SI]));
    memcpy(Q,buf.ptr,SMALL*sizeof(uint));
  for(i=0;i!=SMALL;++i)
    printf("%u\t%u\t%016lx\t%d\t%d\n",(unsigned)P[i],(unsigned)Q[i],
           (unsigned long)A[P[i]][0],
           A[P[i]][0]==A[Q[i]][0],
           Av[i]==A[P[i]][0]);
  printf("\n");
  printf("radix sort:\n");
  for(i=0;i!=NUM;++i) Q[i]=NUM-1-i;
  sortv_long(Av,  &A[0][0],NUM,sizeof(ulong[SI]), &buf);
  sortp_long(&buf,0, &A[0][0],NUM,sizeof(ulong[SI]));
    memcpy(P,buf.ptr,NUM*sizeof(uint));
  memcpy(buf.ptr,Q,NUM*sizeof(uint));
  sortp_long(&buf,1, &A[0][0],NUM,sizeof(ulong[SI]));
    memcpy(Q,buf.ptr,NUM*sizeof(uint));
  for(i=0;i!=NUM;++i)
    printf("%u\t%u\t%016lx\t%d\t%d\n",(unsigned)P[i],(unsigned)Q[i],
           (unsigned long)A[P[i]][0],
           A[P[i]][0]==A[Q[i]][0],
           Av[i]==A[P[i]][0]);

  printf("\nsmall integers:\n");
  printf("\n");

  printf("heap sort:\n");
  for(i=0;i!=SMALL;++i) Q[i]=SMALL-1-i;
  sortv(Q,  Q,SMALL,sizeof(uint), &buf);
  for(i=0;i!=SMALL;++i) printf("\t%u\n",(unsigned)Q[i]);

  printf("merge sort:\n");
  for(i=0;i!=SMALL;++i) Q[i]=SMALL-1-i;
  sortv(Bv,  &B[0][0],SMALL,sizeof(uint[SI]), &buf);
  sortp(&buf,0, &B[0][0],SMALL,sizeof(uint[SI]));
    memcpy(P,buf.ptr,SMALL*sizeof(uint));
  memcpy(buf.ptr,Q,SMALL*sizeof(uint));
  sortp(&buf,1, &B[0][0],SMALL,sizeof(uint[SI]));
    memcpy(Q,buf.ptr,SMALL*sizeof(uint));
  for(i=0;i!=SMALL;++i)
    printf("%u\t%u\t%016lx\t%d\t%d\n",(unsigned)P[i],(unsigned)Q[i],
           (unsigned long)B[P[i]][0],
           B[P[i]][0]==B[Q[i]][0],
           B[P[i]][0]==Bv[i]);
  printf("\n");
  printf("radix sort:\n");
  for(i=0;i!=NUM;++i) Q[i]=NUM-1-i;
  sortv(Bv,  &B[0][0],NUM,sizeof(uint[SI]), &buf);
  sortp(&buf,0, &B[0][0],NUM,sizeof(uint[SI]));
    memcpy(P,buf.ptr,NUM*sizeof(uint));
  memcpy(buf.ptr,Q,NUM*sizeof(uint));
  sortp(&buf,1, &B[0][0],NUM,sizeof(uint[SI]));
    memcpy(Q,buf.ptr,NUM*sizeof(uint));
  for(i=0;i!=NUM;++i)
    printf("%u\t%u\t%016lx\t%d\t%d\n",(unsigned)P[i],(unsigned)Q[i],
           (unsigned long)B[P[i]][0],
           B[P[i]][0]==B[Q[i]][0],
           B[P[i]][0]==Bv[i]);
  buffer_free(&buf);
  return 0;
}

