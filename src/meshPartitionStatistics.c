
#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

void meshPartitionStatistics(mesh_t *mesh){

  /* get MPI rank and size */
  int rank, size;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  
  /* now gather statistics on connectivity between processes */
  int *comms = (int*) calloc(size, sizeof(int));
  int Ncomms = 0;

  /* count elements with neighbors on each other rank ranks */
  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      if(mesh->EToP[e*mesh->Nfaces+f]!=-1){
	++comms[mesh->EToP[e*mesh->Nfaces+f]];
	++Ncomms;
      }
    }
  }

  int Nmessages = 0;
  for(int r=0;r<size;++r)
    if(comms[r]>0)
      ++Nmessages;

  for(int r=0;r<size;++r){
    MPI_Barrier(MPI_COMM_WORLD);
    if(r==rank){
      fflush(stdout);
      printf("r: %02d [", rank);
      for(int s=0;s<size;++s){
	printf(" %04d", comms[s]);
      }
      printf("] (Nelements=%d, Nmessages=%d, Ncomms=%d)\n", mesh->Nelements,Nmessages, Ncomms);
      fflush(stdout);
    }
  }
  
  free(comms);
}
