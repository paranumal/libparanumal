
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"

void meshMRABSetupP2D(mesh2D *mesh, dfloat *EToDT, int maxLevels) {

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  //find global min and max dt
  dfloat dtmin, dtmax;
  dtmin = EToDT[0];
  dtmax = EToDT[0];
  for (int e=1;e<mesh->Nelements;e++) {
    dtmin = mymin(dtmin,EToDT[e]);
    dtmax = mymax(dtmax,EToDT[e]);
  }
  dfloat dtGmin, dtGmax;
  MPI_Allreduce(&dtmin, &dtGmin, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);    
  MPI_Allreduce(&dtmax, &dtGmax, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);    


  if (rank==0) {
    printf("----------- MRAB Setup ------------------------------\n");
    printf("dtmin = %g, dtmax = %g\n", dtGmin, dtGmax);
  }

  //number of levels
  mesh->MRABNlevels = mymin(floor(log2(dtGmax/dtGmin))+1,maxLevels);

  //shift dtGmin so that we have an integer number of steps
  mesh->NtimeSteps = mesh->finalTime/(pow(2,mesh->MRABNlevels-1)*dtGmin);
  dtGmin = mesh->finalTime/(pow(2,mesh->MRABNlevels-1)*mesh->NtimeSteps);

  mesh->dt = dtGmin; 

  //compute the level of each element
  mesh->MRABlevel = (int *) calloc(mesh->Nelements+mesh->totalHaloPairs,sizeof(int));
  int *MRABsendBuffer;
  if (mesh->totalHaloPairs) 
    MRABsendBuffer = (int *) calloc(mesh->totalHaloPairs,sizeof(int));
  for(int lev=0; lev<mesh->MRABNlevels; lev++){             
    dfloat dtlev = dtGmin*pow(2,lev);   
    for(int e=0;e<mesh->Nelements;++e){
      if(EToDT[e] >=dtlev) 
        mesh->MRABlevel[e] = lev;
    }
  }

  //enforce one level difference between neighbours
  for (int lev=0; lev < mesh->MRABNlevels; lev++){
    if (mesh->totalHaloPairs) 
      meshHaloExchange(mesh, sizeof(int), mesh->MRABlevel, MRABsendBuffer, mesh->MRABlevel+mesh->Nelements);
    for (int e =0; e<mesh->Nelements;e++) {
      if (mesh->MRABlevel[e] > lev+1) { //find elements at least 2 levels higher than lev
        for (int f=0;f<mesh->Nfaces;f++) { //check for a level lev neighbour
          int eP = mesh->EToE[mesh->Nfaces*e+f];
          if (eP > -1) 
            if (mesh->MRABlevel[eP] == lev)
              mesh->MRABlevel[e] = lev + 1;  //if one exists, lower the level of this element to lev-1
        }
      }
    }
  }
  if (mesh->totalHaloPairs) free(MRABsendBuffer);

  //this could change the number of MRAB levels there are, so find the new max level
  mesh->MRABNlevels = 0;
  for (int e=0;e<mesh->Nelements;e++)
    mesh->MRABNlevels = (mesh->MRABlevel[e]>mesh->MRABNlevels) ? mesh->MRABlevel[e] : mesh->MRABNlevels;
  mesh->MRABNlevels++;
  int localNlevels = mesh->MRABNlevels;
  MPI_Allreduce(&localNlevels, &(mesh->MRABNlevels), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);    
  mesh->NtimeSteps = mesh->finalTime/(pow(2,mesh->MRABNlevels-1)*dtGmin);

  //now we need to perform a weighted repartitioning of the mesh to optimize MRAB
  if (size>1) {
    dfloat *weights = (dfloat *) calloc(mesh->Nelements,sizeof(dfloat));
    
    //TODO refine this estimate
    for (int e=0; e<mesh->Nelements;e++) {
      int N = mesh->N[e];
      weights[e] = (mesh->Np[N]*mesh->Np[N])*pow(2,mesh->MRABNlevels-mesh->MRABlevel[e]);
    }
    
    if (rank==0) printf("Repartitioning for MRAB...\n");
    meshMRABWeightedPartitionTriP2D(mesh,weights,mesh->MRABNlevels, mesh->MRABlevel);
  }

  //construct element and halo lists
  mesh->MRABelementIds = (int **) calloc(mesh->MRABNlevels,sizeof(int*));
  mesh->MRABhaloIds = (int **) calloc(mesh->MRABNlevels,sizeof(int*));
  
  mesh->MRABNelements = (int *) calloc(mesh->MRABNlevels,sizeof(int));
  mesh->MRABNhaloElements = (int *) calloc(mesh->MRABNlevels,sizeof(int));

  mesh->MRABNelP = (int **) calloc(mesh->MRABNlevels,sizeof(int*));
  mesh->MRABNhaloEleP = (int **) calloc(mesh->MRABNlevels,sizeof(int*));
  mesh->MRABelIdsP = (int ***) calloc(mesh->MRABNlevels,sizeof(int**));
  mesh->MRABhaloIdsP = (int ***) calloc(mesh->MRABNlevels,sizeof(int**));
  for (int lev=0;lev<mesh->MRABNlevels;lev++) {
    mesh->MRABNelP[lev] = (int *) calloc(mesh->NMax+1,sizeof(int));
    mesh->MRABNhaloEleP[lev] = (int *) calloc(mesh->NMax+1,sizeof(int));
    mesh->MRABelIdsP[lev] = (int **) calloc(mesh->NMax+1,sizeof(int *));
    mesh->MRABhaloIdsP[lev] = (int **) calloc(mesh->NMax+1,sizeof(int *)); 
  }  

  for (int e=0;e<mesh->Nelements;e++) {
    mesh->MRABNelements[mesh->MRABlevel[e]]++;
    mesh->MRABNelP[mesh->MRABlevel[e]][mesh->N[e]]++;
    for (int f=0;f<mesh->Nfaces;f++) { 
      int eP = mesh->EToE[mesh->Nfaces*e+f];
      if (eP > -1) {
        if (mesh->MRABlevel[eP] == mesh->MRABlevel[e]-1) {//check for a level lev-1 neighbour
          mesh->MRABNhaloElements[mesh->MRABlevel[e]]++;
          mesh->MRABNhaloEleP[mesh->MRABlevel[e]][mesh->N[e]]++;
          break;
        }
      }
    }
  }

  for (int lev =0;lev<mesh->MRABNlevels;lev++){
    mesh->MRABelementIds[lev] = (int *) calloc(mesh->MRABNelements[lev],sizeof(int));
    mesh->MRABhaloIds[lev] = (int *) calloc(mesh->MRABNhaloElements[lev],sizeof(int));
    int cnt  =0;
    int cnt2 =0;

    int pcnt[mesh->NMax+1];
    int pcnt2[mesh->NMax+1];

    for (int p=0;p<=mesh->NMax;p++) {
      pcnt[p] = 0;
      pcnt2[p] = 0;
      mesh->MRABelIdsP[lev][p]   = (int *) calloc(mesh->MRABNelP[lev][p],sizeof(int));
      mesh->MRABhaloIdsP[lev][p] = (int *) calloc(mesh->MRABNhaloEleP[lev][p],sizeof(int)); 
    }

    for (int e=0;e<mesh->Nelements;e++){
      if (mesh->MRABlevel[e] == lev) {
        mesh->MRABelementIds[lev][cnt++] = e;
        mesh->MRABelIdsP[lev][mesh->N[e]][pcnt[mesh->N[e]]++] = e;
      
        for (int f=0;f<mesh->Nfaces;f++) { 
          int eP = mesh->EToE[mesh->Nfaces*e+f];
          if (eP > -1) {
            if (mesh->MRABlevel[eP] == lev-1) {//check for a level lev-1 neighbour
              mesh->MRABhaloIds[lev][cnt2++] = e;
              mesh->MRABhaloIdsP[lev][mesh->N[e]][pcnt2[mesh->N[e]]++] = e;
              break;
            }
          }
        }
      }
    }
  }

  //offset index
  mesh->MRABshiftIndex = (int *) calloc(mesh->MRABNlevels,sizeof(int));

  if (rank==0){
    printf("| Rank | Level | Nelements | Level/Level Boundary Elements | \n");
    printf("------------------------------------------------------------\n");
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for (int r =0;r<size;r++) {
    if (r==rank) {
      for (int lev =0; lev<mesh->MRABNlevels; lev++) 
        printf("|  %d,    %d,      %d,        %d     \n", rank, lev, mesh->MRABNelements[lev], mesh->MRABNhaloElements[lev]);
      printf("------------------------------------------------------------\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}