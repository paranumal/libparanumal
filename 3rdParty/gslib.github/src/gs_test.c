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
#include "gs_defs.h"
#include "gs.h"

typedef double T;
const gs_dom dom = gs_double;

static void test(const struct comm *comm)
{
  struct gs_data *gsh;
  const uint np = comm->np;
  slong *id = tmalloc(slong,np+4);
  T *v = tmalloc(T,np+4);
  uint i;
  id[0] = -(slong)(np+10+3*comm->id);
  for(i=0;i<np;++i) id[i+1] = -(sint)(i+1);
  id[np+1] = comm->id+1;
  id[np+2] = comm->id+1;
  id[np+3] = np-comm->id;
  gsh = gs_setup(id,np+4,comm,0,gs_auto,1);
  free(id);
  
  for(i=0;i<np+4;++i) v[i] = 1;
  gs(v,dom,gs_add,0,gsh,0);
  if(comm->id==0) for(i=0;i<np+4;++i) printf("%g\n",v[i]);
  if(comm->id==0) printf("\n");
  for(i=0;i<np+4;++i) v[i] = 1;
  gs(v,dom,gs_add,1,gsh,0);
  if(comm->id==0) for(i=0;i<np+4;++i) printf("%g\n",v[i]);

  gs_free(gsh);
  free(v);
}

int main(int narg, char *arg[])
{
  comm_ext world; int np;
  struct comm comm;
  
#ifdef MPI
  MPI_Init(&narg,&arg);
  world = MPI_COMM_WORLD;
  MPI_Comm_size(world,&np);
#else
  world=0, np=1;
#endif

  comm_init(&comm,world);

  test(&comm);
  
  comm_free(&comm);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
