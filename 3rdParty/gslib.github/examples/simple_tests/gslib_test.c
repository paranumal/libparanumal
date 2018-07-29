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
  FILE *inp;
  struct gs_data *gsh;
  struct comm comm;
  int *targetCoreIndexing,index;
  slong *sendbuf,*recvbuf;
  int *sendcounts,*displs,*duplicate_count;
  int ret,totalElements,maxNp,v1,v2,v3,v4;
  int arrayOffset,kOld,localBufSpace,maxArray;
  int elementIndex,targetCore,i,j,nid,k;
  char buffer[1024];
  char inpFile[128];
  int **mat,fail;
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

  if(nid==0){
    inp = fopen(arg[1],"r");
    mat = (int **) malloc(sizeof(int *)*np);
    //Only supports 2d right now
    fgets(buffer,1024,inp);
    sscanf(buffer,"%d%d%d%d%d%d%d",&totalElements,&maxArray,&v2,&maxNp,&i,&j,&k);
    maxArray = maxArray+1; //Index off of 1, not 0
    MPI_Bcast(&maxArray,1,MPI_INT,0,world);

    duplicate_count = malloc(sizeof(int)*maxArray);

    //scanf("%d%d%d%d%d%d%d",&totalElements,&v1,&v2,&maxNp,&i,&j,&k);
    //    fscanf(inp,"%d%d%d%d%d%d%d",&totalElements,&v1,&v2,&maxNp,&i,&j,&k);

    //Allocate memory for each core
    //Allocate memory for the elements
    for(i=0;i<np;i++) mat[i] = (int*) malloc(sizeof(int)*4*totalElements);

    //Set an array for indexing the elements in the matrix
    targetCoreIndexing = (int *) malloc(sizeof(int)*np);
    for(i=0;i<np;i++) {
      targetCoreIndexing[i] = 0;
    }

    //read in the data
    for(i=0;i<totalElements;i++){

      fgets(buffer,1024,inp);
      sscanf(buffer,"%d%d%d%d%d", &elementIndex, &v1, &v2, &v3, &v4);
      //scanf("%d%d%d%d%d", &elementIndex, &v1, &v2, &v3, &v4);
      //      ret = fscanf(inp,"%d%d%d%d%d", &elementIndex, &v1, &v2, &v3, &v4);
      //if(ret == EOF) printf("EOF\n");

      targetCore = elementIndex%np;

      arrayOffset = targetCoreIndexing[targetCore];

      //Store the vertices

      k = 0+4*arrayOffset;
      mat[targetCore][k] = v1;

      k++;
      mat[targetCore][k] = v2;

      k++;
      mat[targetCore][k] = v3;

      k++;
      mat[targetCore][k] = v4;

      targetCoreIndexing[targetCore] = targetCoreIndexing[targetCore]+1;
    }
    
    sendbuf    = malloc(sizeof(long)*4*totalElements);
    sendcounts = malloc(sizeof(int)*np);
    displs     = malloc(sizeof(int)*np);
    k=0;
    kOld=0;

    //Distribute data

    for(i=0;i<np;i++) {
      //Fill the send buffer
      for(j=0;j<targetCoreIndexing[i];j++){
        sendbuf[k] = mat[i][0+4*j];
        k++;
        sendbuf[k] = mat[i][1+4*j];
        k++;
        sendbuf[k] = mat[i][2+4*j];
        k++;
        sendbuf[k] = mat[i][3+4*j];
        k++;
      }
      sendcounts[i] = (k-kOld);
      if(i!=0){
        MPI_Send(&sendcounts[i],1,MPI_INT,i,0,world);
      }
      displs[i]     = kOld;
      kOld = k;
    }

    for(i=0;i<maxArray;i++){
      duplicate_count[i] = 0;
    }

    //Count duplicates
    for(i=0;i<4*totalElements;i++){
      duplicate_count[sendbuf[i]]++;
    }

    //Bcast duplicates
    MPI_Bcast(duplicate_count,maxArray,MPI_INT,0,world);

    localBufSpace = sendcounts[0];
    recvbuf = malloc(sizeof(long)*localBufSpace);
    v = malloc(sizeof(double)*localBufSpace);
    MPI_Scatterv(sendbuf,sendcounts,displs,MPI_LONG,recvbuf,localBufSpace,MPI_LONG,0,world);
    if(np!=1){
      for(i=0;i<sendcounts[0];i++){
        recvbuf[i] = sendbuf[i];
      }
    }
  } else {
    MPI_Recv(&localBufSpace,1,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Bcast(&maxArray,1,MPI_INT,0,world);
    duplicate_count = malloc(sizeof(int)*maxArray);
    MPI_Bcast(duplicate_count,maxArray,MPI_INT,0,world);
    //printf("Nid: %d maxEle: %d\n",nid,localBufSpace);
    recvbuf = malloc(sizeof(long)*localBufSpace);
    v = malloc(sizeof(double)*localBufSpace);
    MPI_Scatterv(sendbuf,sendcounts,displs,MPI_LONG,recvbuf,localBufSpace,MPI_LONG,0,world);
  }

  gsh = gs_setup(recvbuf,localBufSpace,&comm,0,gs_pairwise,1);
#pragma acc enter data create(v[0:localBufSpace])
#pragma acc enter data copyin(recvbuf[0:localBufSpace])

#pragma acc parallel loop present(v[0:localBufSpace],recvbuf[0:localBufSpace])
  for(i=0;i<localBufSpace;i++){
    v[i] = recvbuf[i];
  }

  gs(v,dom,gs_add,0,gsh,0);

#pragma acc update host(v[0:localBufSpace])
  fail = 0;
  //Check v
  for(i=0;i<localBufSpace;i++){
    if(v[i]!=duplicate_count[recvbuf[i]]*recvbuf[i]){
      printf("Add failure on core %d index %d\n",nid,i);
      printf("v[%d] %f recv %d %d\n",i,v[i],duplicate_count[recvbuf[i]],recvbuf[i]);
      fail = 1;
    }
  }
  
  if(fail==0) printf("Add success! on %d\n",nid);
  //Fill v
#pragma acc parallel loop present(v[0:localBufSpace],recvbuf[0:localBufSpace])
  for(i=0;i<localBufSpace;i++){
    v[i] = recvbuf[i];
  }

  gs(v,dom,gs_mul,0,gsh,0);

#pragma acc update host(v[0:localBufSpace])
  fail = 0;
  //Check v
  for(i=0;i<localBufSpace;i++){
    if(v[i]!=pow(recvbuf[i],duplicate_count[recvbuf[i]])){
      printf("Mult failure on core %d index %d\n",nid,i);
      printf("v[%d] %f recv %d %d\n",i,v[i],recvbuf[i],duplicate_count[recvbuf[i]]);
      fail = 1;
    }
  }
  
  if(fail==0) printf("Mult success! on %d\n",nid);

  comm_free(&comm);
  if(nid==0){
    for(i=0;i<np;i++) free(mat[i]);
    free(mat);
    free(targetCoreIndexing);
    free(sendbuf);
    free(sendcounts);
    free(displs);
  }
  free(recvbuf);

#ifdef MPI
  MPI_Barrier(world);
  MPI_Finalize();
#endif

  return 0;
}
