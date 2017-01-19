#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../../src/c99.h"
#include "../../src/name.h"
#include "../../src/fail.h"
#include "../../src/types.h"
#include "../../src/comm.h"
#include "../../src/mem.h"
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
  int *localData,*recvBuf;
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

  glo_num = malloc(sizeof(slong)*np);

  for(i=1;i<=np;i++){
    j = nid+1;
    if(j>=i){
      glo_num[i-1] = (i-1)*np-i*(i-1)/2 + j-1;
    } else {
      glo_num[i-1] = (j-1)*np - j*(j-1)/2+i-1;
    }
    if(j==i){
      glo_num[i-1] = 0;
    }
  }

  gsh = gs_setup(glo_num,np,&comm,0,gs_auto,1);

  nsamples  = 10000;
  localData = malloc(sizeof(int)*np);
  recvBuf = malloc(sizeof(int)*np);

  MPI_Barrier(world);
  t1 = MPI_Wtime();
  for(j=0;j<nsamples;j++){
    for(i=0;i<np;i++){
      localData[i] = nid+i;
      recvBuf[i] = nid+i;
    }
    
    gs(recvBuf,gs_int,gs_add,0,gsh,0);

    for(i=0;i<np;i++){
      recvBuf[i] = recvBuf[i] - localData[i];
    }
  }
  MPI_Barrier(world);
  t2 = MPI_Wtime();
  gs_time = t2 - t1;
  MPI_Barrier(world);
  t1 = MPI_Wtime();
  for(j=0;j<nsamples;j++){
    for(i=0;i<np;i++){
      localData[i] = nid+i;
    }
    
    MPI_Alltoall(localData,1,MPI_INT,recvBuf,1,MPI_INT,world);
  }

  MPI_Barrier(world);
  t2 = MPI_Wtime();
  mpi_time = t2 - t1;

  if(nid==0)printf("gs_time: %f mpi_time: %f\n",gs_time,mpi_time);
#ifdef MPI
  MPI_Barrier(world);
  //  MPI_Finalize();
#endif

  return 0;
}
