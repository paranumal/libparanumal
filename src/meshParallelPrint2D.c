#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"
#include "mesh2D.h"

void meshParallelPrint2D(mesh2D *mesh){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  printf("rank %d: Nelements=" intFormat " Nnodes=" intFormat "\n", 
	 rank, mesh->Nelements, mesh->Nnodes);
  
#if 0
  printf("EToV:\n");
  for(int e=0;e<mesh->Nelements;++e){
    printf("%d %d %d\n", 
	   mesh->EToV[e*mesh->Nverts+0],
	   mesh->EToV[e*mesh->Nverts+1],
	   mesh->EToV[e*mesh->Nverts+2]);
  }
#endif

  int *otherNelements = (int*) calloc(size, sizeof(int));
  MPI_Allgather(&(mesh->Nelements), 1, MPI_int,
		otherNelements, 1, MPI_int, 
		MPI_COMM_WORLD);
  
  int *elementStarts = (int*) calloc(size, sizeof(int));
  for(int r=1;r<size;++r){
    elementStarts[r] = elementStarts[r-1]+otherNelements[r-1];
  }

  for(int r1=0;r1<size;++r1){
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==r1){
      fflush(stdout);
      if(r1==0)
	printf("EToE:\n");
      for(int e1=0;e1<mesh->Nelements;++e1){
	int id = e1*mesh->Nfaces;
	for(int f1=0;f1<mesh->Nfaces;++f1){
	  int e2 = mesh->EToE[id+f1];
	  int f2 = mesh->EToF[id+f1];
	  int r2 = mesh->EToP[id+f1];
	  if(e2==-1 || f2==-1) 
	    printf("(" intFormat " " intFormat ")=>X (" intFormat "," intFormat ")\n", 
		   e1+elementStarts[r1], f1, e2, f2);
	  else{
	    
	    if(r2!=-1)
	      e2 += elementStarts[r2];
	    else
	      e2 += elementStarts[r1];
	    
	    
	    printf("(" intFormat " " intFormat ")=>(" intFormat " " intFormat ")\n", 
		   e1+elementStarts[r1], f1, e2, f2);
	  }
	}
	fflush(stdout);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
}
