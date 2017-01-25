#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../../src/c99.h"
#include "../../src/name.h"
#include "../../src/fail.h"
#include "../../src/types.h"
#include "../../src/mem.h"
#include "../../src/comm.h"
#include "../../src/gs_defs.h"
#include "../../src/gs.h"

typedef double T;
const gs_dom dom = gs_double;

int main(int narg, char *arg[])
{
  comm_ext world; int np;
  double t1,t2,gs_time,mpi_time;
  struct gs_data *gsh;
  struct comm comm;
  double *localData,*recvBuf;
  slong *glo_num;
  int i,j,nid,nsamples;
  T *v;

#ifdef MPI
  MPI_Init(&narg,&arg);
  world = MPI_COMM_WORLD;
  MPI_Comm_size(world,&np);

#else
  world=0, np=1;
#endif

  comm_init(&comm,world);
  MPI_Comm_rank(world,&nid);

  glo_num = malloc(sizeof(slong)*1);
  glo_num[0] = 1;

  gsh = gs_setup(glo_num,1,&comm,0,gs_auto,1);

  nsamples  = 500000;
  localData = malloc(sizeof(int)*1);
  recvBuf = malloc(sizeof(int)*1);

  MPI_Barrier(world);
  t1 = MPI_Wtime();
  for(j=0;j<nsamples;j++){
    recvBuf[0] = ((nid+1.0)/np)*2;
    
    /* gs_irecv(recvBuf,gs_int,gs_add,0,gsh,0); */
    /* gs_isend(recvBuf,gs_int,gs_add,0,gsh,0); */
    /* gs_wait(recvBuf,gs_int,gs_add,0,gsh,0); */
    gs(recvBuf,gs_double,gs_add,0,gsh,0);
  }

  MPI_Barrier(world);
  t2 = MPI_Wtime();
  gs_time = t2 - t1;
  MPI_Barrier(world);
  t1 = MPI_Wtime();
  for(j=0;j<nsamples;j++){
    localData[0] = ((nid+1.0)/np)*2;
    MPI_Allreduce(localData,localData,1,MPI_DOUBLE,MPI_SUM,world);
  }

  MPI_Barrier(world);
  t2 = MPI_Wtime();
  mpi_time = t2 - t1;
  if(nid==0) printf("gs_sum: %f mpi_sum %f\n",recvBuf[0],localData[0]);
  if(nid==0)printf("gs_time: %f mpi_time: %f\n",gs_time,mpi_time);
#ifdef MPI
  MPI_Barrier(world);
  MPI_Finalize();
#endif

  return 0;
}
