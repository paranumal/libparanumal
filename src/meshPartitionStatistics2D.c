
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"
 

void meshPartitionStatistics2D(mesh2D *mesh){

  /* get MPI rank and size */
  int rank, size;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  /* make sure mesh is connected */
  meshParallelConnect2D(mesh);
  
  /* now gather statistics on connectivity between processes */
  iint *comms = (iint*) calloc(size, sizeof(iint));
  iint Ncomms = 0;

  /* count elements with neighbors on each other rank ranks */
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint f=0;f<mesh->Nfaces;++f){
      if(mesh->EToP[e*mesh->Nfaces+f]!=-1){
	++comms[mesh->EToP[e*mesh->Nfaces+f]];
	++Ncomms;
      }
    }
  }

  iint Nmessages = 0;
  for(iint r=0;r<size;++r)
    if(comms[r]>0)
      ++Nmessages;

  for(iint r=0;r<size;++r){
    MPI_Barrier(MPI_COMM_WORLD);
    if(r==rank){
      fflush(stdout);
      printf("r: %02d [", rank);
      for(iint s=0;s<size;++s){
	printf(" %04d", comms[s]);
      }
      printf("] (Nmessages=%d, Ncomms=%d)\n", Nmessages, Ncomms);
      fflush(stdout);
    }
  }
  
  free(comms);

}
