#include <stdlib.h>
#include <stdio.h>

#include "mesh3D.h"

void meshParallelPrint3D(mesh3D *mesh){

  int rank, size;
  rank = mesh->rank;
  size = mesh->size;

  printf("rank %d: Nelements=" dlongFormat " Nnodes=" hlongFormat "\n", 
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

  dlong *otherNelements = (dlong*) calloc(size, sizeof(dlong));
  MPI_Allgather(&(mesh->Nelements), 1, MPI_DLONG,
                    otherNelements, 1, MPI_DLONG, 
                    mesh->comm);
  
  hlong *elementStarts = (hlong*) calloc(size, sizeof(hlong));
  for(int r=1;r<size;++r){
    elementStarts[r] = elementStarts[r-1]+otherNelements[r-1];
  }

  for(int r1=0;r1<size;++r1){
    MPI_Barrier(mesh->comm);
    if(rank==r1){
      fflush(stdout);
      if(r1==0)
        printf("EToE:\n");
      for(dlong e1=0;e1<mesh->Nelements;++e1){
        dlong id = e1*mesh->Nfaces;
        for(int f1=0;f1<mesh->Nfaces;++f1){
          hlong e2 = (hlong) mesh->EToE[id+f1];
          int f2 = mesh->EToF[id+f1];
          int r2 = mesh->EToP[id+f1];
          if(e2==-1 || f2==-1) 
            printf("(" hlongFormat " " "%d" ")=>X (" hlongFormat "," "%d" ")\n", 
                   e1+elementStarts[r1], f1, e2, f2);
          else{
            
            if(r2!=-1)
              e2 += elementStarts[r2];
            else
              e2 += elementStarts[r1];
            
            
            printf("(" hlongFormat " " "%d" ")=>(" hlongFormat " " "%d" ")\n", 
                   e1+elementStarts[r1], f1, e2, f2);
          }
        }
        fflush(stdout);
      }
    }
    MPI_Barrier(mesh->comm);
  }
  free(otherNelements);
  free(elementStarts);
}
