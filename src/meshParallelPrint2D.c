#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"
#include "mesh2D.h"

void meshParallelPrint2D(mesh2D *mesh){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  printf("rank %d: Nelements=" iintFormat " Nnodes=" iintFormat "\n", 
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

  iint *otherNelements = (iint*) calloc(size, sizeof(iint));
  MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT,
		otherNelements, 1, MPI_IINT, 
		MPI_COMM_WORLD);
  
  iint *elementStarts = (iint*) calloc(size, sizeof(iint));
  for(iint r=1;r<size;++r){
    elementStarts[r] = elementStarts[r-1]+otherNelements[r-1];
  }

  for(int r1=0;r1<size;++r1){
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==r1){
      fflush(stdout);
      if(r1==0)
	printf("EToE:\n");
      for(iint e1=0;e1<mesh->Nelements;++e1){
	iint id = e1*mesh->Nfaces;
	for(iint f1=0;f1<mesh->Nfaces;++f1){
	  iint e2 = mesh->EToE[id+f1];
	  iint f2 = mesh->EToF[id+f1];
	  iint r2 = mesh->EToP[id+f1];
	  if(e2==-1 || f2==-1) 
	    printf("(" iintFormat " " iintFormat ")=>X (" iintFormat "," iintFormat ")\n", 
		   e1+elementStarts[r1], f1, e2, f2);
	  else{
	    
	    if(r2!=-1)
	      e2 += elementStarts[r2];
	    else
	      e2 += elementStarts[r1];
	    
	    
	    printf("(" iintFormat " " iintFormat ")=>(" iintFormat " " iintFormat ")\n", 
		   e1+elementStarts[r1], f1, e2, f2);
	  }
	}
	fflush(stdout);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
}
