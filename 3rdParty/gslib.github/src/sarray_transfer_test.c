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
#include "c99.h"
#include "name.h"
#include "fail.h"
#include "types.h"
#include "comm.h"
#include "mem.h"
#include "sort.h"
#include "sarray_sort.h"
#include "crystal.h"
#include "sarray_transfer.h"

typedef struct {
  double d;
  ulong l,l2;
  uint i;
  uint p;
} r_work;

int main(int narg, char *arg[])
{
  comm_ext world; int np;
  struct comm comm;
  struct crystal crystal;
  struct array A, A0=null_array; r_work *row, *row_0;
  uint i;
#ifdef MPI
  MPI_Init(&narg,&arg);
  world = MPI_COMM_WORLD;
  MPI_Comm_size(world,&np);
#else
  world=0, np=1;
#endif

  comm_init(&comm,world);
  crystal_init(&crystal,&comm);

  array_init(r_work,&A,np*3), A.n=np*3, row=A.ptr;
  for(i=0;i<A.n;++i) {
    row[i].i = rand();
    row[i].l = row[i].l2 = rand();
    row[i].p = rand()%np;
    row[i].d = rand()/(double)rand();
  }
  
  sarray_sort_3(r_work,row,A.n, i,0, l,1, p,0, &crystal.data);
  
  for(i=0;i<A.n;++i)
    printf("%02d send -> %02d: %08x %08x %d %g\n",
      (int)comm.id,(int)row[i].p,(int)row[i].i,
      (int)row[i].l,(int)row[i].p,row[i].d);
  
  array_cat(r_work,&A0, row,A.n);
  
  sarray_transfer(r_work,&A, p,1, &crystal);

  row=A.ptr;
  for(i=0;i<A.n;++i)
    printf("%02d recv <- %02d: %08x %08x %d %g\n",
      (int)comm.id,(int)row[i].p,(int)row[i].i,
      (int)row[i].l,(int)row[i].p,row[i].d);

  sarray_transfer(r_work,&A, p,1, &crystal);
  sarray_sort_3(r_work,row,A.n, i,0, l,1, p,0, &crystal.data);
  if(A.n!=A0.n)
    fail(1,__FILE__,__LINE__,"final array has different length than original");
  row=A.ptr, row_0=A0.ptr;
  for(i=0;i<A.n;++i)
    if(   row[i].d != row_0[i].d
       || row[i].l != row_0[i].l
       || row[i].l2!= row_0[i].l2
       || row[i].i != row_0[i].i
       || row[i].p != row_0[i].p)
      fail(1,__FILE__,__LINE__,"final array differs from original");
      
  array_free(&A0);
  array_free(&A);
  crystal_free(&crystal);

  fflush(stdout); comm_barrier(&comm);
  if(comm.id==0) printf("tests passed\n"), fflush(stdout);
  
  comm_free(&comm);
  
#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
