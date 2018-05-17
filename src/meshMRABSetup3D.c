
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mesh3D.h"

dfloat meshMRABSetup3D(mesh3D *mesh, dfloat *EToDT, int maxLevels, dfloat finalTime) {

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  //find global min and max dt
  dfloat dtmin, dtmax;
  dtmin = EToDT[0];
  dtmax = EToDT[0];
  for (dlong e=1;e<mesh->Nelements;e++) {
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
  mesh->NtimeSteps = finalTime/(pow(2,mesh->MRABNlevels-1)*dtGmin);
  dtGmin = finalTime/(pow(2,mesh->MRABNlevels-1)*mesh->NtimeSteps);

 

  //compute the level of each element
  mesh->MRABlevel = (dlong *) calloc(mesh->Nelements+mesh->totalHaloPairs,sizeof(int));
  int *MRABsendBuffer = (int *) calloc(mesh->totalHaloPairs,sizeof(int));
  for(int lev=0; lev<mesh->MRABNlevels; lev++){             
    dfloat dtlev = dtGmin*pow(2,lev);   
    for(dlong e=0;e<mesh->Nelements;++e){
      if(EToDT[e] >=dtlev) 
        mesh->MRABlevel[e] = lev;
    }
  }

  //enforce one level difference between neighbours
  for (int lev=0; lev < mesh->MRABNlevels; lev++){
    if (mesh->totalHaloPairs) 
      meshHaloExchange(mesh, sizeof(int), mesh->MRABlevel, MRABsendBuffer, mesh->MRABlevel+mesh->Nelements);
    for (dlong e =0; e<mesh->Nelements;e++) {
      if (mesh->MRABlevel[e] > lev+1) { //find elements at least 2 levels higher than lev
        for (int f=0;f<mesh->Nfaces;f++) { //check for a level lev neighbour
          dlong eP = mesh->EToE[mesh->Nfaces*e+f];
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
  for (dlong e=0;e<mesh->Nelements;e++)
    mesh->MRABNlevels = (mesh->MRABlevel[e]>mesh->MRABNlevels) ? mesh->MRABlevel[e] : mesh->MRABNlevels;
  mesh->MRABNlevels++;
  int localNlevels = mesh->MRABNlevels;
  MPI_Allreduce(&localNlevels, &(mesh->MRABNlevels), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);    
  mesh->NtimeSteps = mesh->finalTime/(pow(2,mesh->MRABNlevels-1)*dtGmin);

  //now we need to perform a weighted repartitioning of the mesh to optimize MRAB
  if (size>1) {
    //for the moment, just weigth the elements by the number or RHS evals per MRAB step
    // TODO: We should make this an input parameter later to handle other problems. 
    dfloat *weights = (dfloat *) calloc(mesh->Nelements,sizeof(dfloat));
    for (dlong e=0; e<mesh->Nelements;e++) {
      weights[e] = pow(2,mesh->MRABNlevels-mesh->MRABlevel[e]);
    }
    
    if (rank==0) printf("Repartitioning for MRAB...\n");
    meshMRABWeightedPartitionTet3D(mesh,weights,mesh->MRABNlevels, mesh->MRABlevel);
  }

  //construct element and halo lists
  mesh->MRABelementIds = (int **) calloc(mesh->MRABNlevels,sizeof(int*));
  mesh->MRABhaloIds = (int **) calloc(mesh->MRABNlevels,sizeof(int*));
  
  mesh->MRABNelements = (int *) calloc(mesh->MRABNlevels,sizeof(int));
  mesh->MRABNhaloElements = (int *) calloc(mesh->MRABNlevels,sizeof(int));

  for (dlong e=0;e<mesh->Nelements;e++) {
    mesh->MRABNelements[mesh->MRABlevel[e]]++;
    for (int f=0;f<mesh->Nfaces;f++) { 
      dlong eP = mesh->EToE[mesh->Nfaces*e+f];
      if (eP > -1) {
        if (mesh->MRABlevel[eP] == mesh->MRABlevel[e]-1) {//check for a level lev-1 neighbour
          mesh->MRABNhaloElements[mesh->MRABlevel[e]]++;
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
    for (dlong e=0;e<mesh->Nelements;e++){
      if (mesh->MRABlevel[e] == lev) {
        mesh->MRABelementIds[lev][cnt++] = e;
      
        for (int f=0;f<mesh->Nfaces;f++) { 
          dlong eP = mesh->EToE[mesh->Nfaces*e+f];
          if (eP > -1) {
            if (mesh->MRABlevel[eP] == lev-1) {//check for a level lev-1 neighbour
              mesh->MRABhaloIds[lev][cnt2++] = e;
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

  return dtGmin;
}