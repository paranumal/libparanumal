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
#include "crystal.h"

int main(int narg, char *arg[])
{
  comm_ext world; int np;
  struct comm comm;
  struct crystal cr;
  uint i,sum, *data, *end;
#ifdef MPI
  MPI_Init(&narg,&arg);
  world = MPI_COMM_WORLD;
  MPI_Comm_size(world,&np);
#else
  world=0, np=1;
#endif

  comm_init(&comm,world);
  
  crystal_init(&cr,&comm);

  cr.data.n = (4+(comm.id&1))*comm.np;
  buffer_reserve(&cr.data,cr.data.n*sizeof(uint));
  data = cr.data.ptr;
  for(i=0;i<comm.np;++i, data+=3+data[2]) {
    data[0] = i, data[1] = comm.id, data[2] = 1;
    data[3] = 2*comm.id;
    if(comm.id&1) data[2] = 2, data[4] = data[3]+1;
  }

#if 0
  data = cr.data.ptr, end = data + cr.data.n;
  for(;data!=end; data+=3+data[2]) {
    uint i;
    printf("%u -> %u:",data[1],data[0]);
    for(i=0;i<data[2];++i) printf(" %u",data[3+i]);
    printf("\n");
  }
#endif
  
  crystal_router(&cr);

#if 0
  printf("\n");
  data = cr.data.ptr, end = data + cr.data.n;
  for(;data!=end; data+=3+data[2]) {
    uint i;
    printf("%u <- %u:",data[0],data[1]);
    for(i=0;i<data[2];++i) printf(" %u",data[3+i]);
    printf("\n");
  }
#endif
  
  if(cr.data.n != comm.np*4 + (comm.np/2))
    fail(1,__FILE__,__LINE__,"failure on %u",comm.id);
  sum = 0;
  data = cr.data.ptr, end = data + cr.data.n;
  for(;data!=end; data+=3+data[2]) {
    sum+=data[1];
    if(data[3]!=data[1]*2)
      fail(1,__FILE__,__LINE__,"failure on %u",comm.id);
    if(data[1]&1 && (data[2]!=2 || data[4]!=data[3]+1))
      fail(1,__FILE__,__LINE__,"failure on %u",comm.id);
  }
  if(sum != comm.np*(comm.np-1)/2)
    fail(1,__FILE__,__LINE__,"failure on %u",comm.id);

  crystal_free(&cr);
  comm_free(&comm);

  diagnostic("",__FILE__,__LINE__,
    "test successful %u/%u",(unsigned)comm.id,(unsigned)comm.np);
  
#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
