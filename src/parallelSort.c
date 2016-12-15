#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"
/* use this for iint */
#include "mesh2D.h"

void mergeLists(size_t sz,
		iint N1, char *v1,
		iint N2, char *v2,
		char *v3,
		int (*compare)(const void *, const void *),
		void (*match)(void *, void *)){
    
  iint n1 = 0, n2 = 0, n3 = 0;
    
  // merge two lists from v1 and v2
  for(n3=0;n3<N1+N2;++n3){
    if(n1<N1 && n2<N2){
      int c = compare(v1+n1*sz,v2+n2*sz);
      if(c==-1){
	memcpy(v3+n3*sz, v1+n1*sz, sz);
	++n1;
      }
      else{
	memcpy(v3+n3*sz, v2+n2*sz, sz);
	++n2;
      }
    }
    else if(n1<N1){
      memcpy(v3+n3*sz, v1+n1*sz, sz);
      ++n1;
    }
    else if(n2<N2){
      memcpy(v3+n3*sz, v2+n2*sz, sz);
      ++n2;
    }
  }
  
  // scan for matches
  for(n3=0;n3<N1+N2-1;++n3){
    if(!compare(v3+n3*sz,v3+(n3+1)*sz)){
      match(v3+n3*sz, v3+(n3+1)*sz);
    }
  }
    
  /* copy result back to v1, v2 */
  memcpy(v1, v3,       N1*sz);
  memcpy(v2, v3+sz*N1, N2*sz);
}

// assumes N is even and the same on all ranks
void parallelSort(iint N, void *vv, size_t sz,
		  int (*compare)(const void *, const void *),
		  void (*match)(void *, void *)
		  ){
   
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
    
  /* cast void * to char * */
  char *v = (char*) vv;

  /* sort faces by their vertex number pairs */
  qsort(v, N, sz, compare);
    
  /* now do progressive merges */
  iint NA=N/2, NB = N/2, NC = N/2;
    
  MPI_Request recvA, recvC;
  MPI_Request sendA, sendC;
  MPI_Status status;
  int tag = 999;
    
  /* temporary buffer for incoming data */
  void *A = (void*) calloc(NA, sz);
  void *B = v;
  void *C = v+NB*sz;
  
  /* temporary space for merge sort */
  void *tmp = (void*) calloc(N, sz);

  /* max and min elements out of place hop one process at each step */
  for(iint step=0;step<size-1;++step){
      
    /* send C, receive A */
    if(rank<size-1)
      MPI_Isend(C, NC*sz, MPI_CHAR,  rank+1, tag, MPI_COMM_WORLD, &sendC);
    if(rank>0)
      MPI_Irecv(A, NA*sz, MPI_CHAR,  rank-1, tag, MPI_COMM_WORLD, &recvA);
      
    if(rank<size-1)
      MPI_Wait(&sendC, &status);
    if(rank>0)
      MPI_Wait(&recvA, &status);
      
    /* merge sort A & B */
    if(rank>0) 
      mergeLists(sz, NA, (char*)A, NB, (char*)B, (char*)tmp, compare, match);
      
    /* send A, receive C */
    if(rank>0)
      MPI_Isend(A, NA*sz, MPI_CHAR, rank-1, tag, MPI_COMM_WORLD, &sendA);
    if(rank<size-1)
      MPI_Irecv(C, NC*sz, MPI_CHAR, rank+1, tag, MPI_COMM_WORLD, &recvC);
      
    if(rank>0)
      MPI_Wait(&sendA, &status);
    if(rank<size-1)
      MPI_Wait(&recvC, &status);
      
    /* merge sort B & C */
    mergeLists(sz, NB, (char*)B, NC, (char*)C, (char*)tmp, compare, match);
      
  }
    
  free(tmp);
  free(A);
}
