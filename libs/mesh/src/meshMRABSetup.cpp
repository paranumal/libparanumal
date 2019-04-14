/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "mesh.hpp"


dfloat mesh_t::MRABSetup(dfloat *EToDT, int maxLevels, dfloat finalTime) {

  //find global min and max dt
  dfloat dtmin, dtmax;
  dtmin = EToDT[0];
  dtmax = EToDT[0];
  for (dlong e=1;e<Nelements;e++) {
    dtmin = mymin(dtmin,EToDT[e]);
    dtmax = mymax(dtmax,EToDT[e]);
  }
  dfloat dtGmin, dtGmax;
  MPI_Allreduce(&dtmin, &dtGmin, 1, MPI_DFLOAT, MPI_MIN, comm);
  MPI_Allreduce(&dtmax, &dtGmax, 1, MPI_DFLOAT, MPI_MIN, comm);


  if (rank==0) {
    printf("----------- MRAB Setup ------------------------------\n");
    printf("dtmin = %g, dtmax = %g\n", dtGmin, dtGmax);
  }

  //number of levels
  MRABNlevels = mymin(floor(log2(dtGmax/dtGmin))+1,maxLevels);

  //shift dtGmin so that we have an integer number of steps
  int NtimeSteps = finalTime/(pow(2,MRABNlevels-1)*dtGmin);
  dtGmin = finalTime/(pow(2,MRABNlevels-1)*NtimeSteps);



  //compute the level of each element
  MRABlevel = (dlong *) calloc(Nelements+totalHaloPairs,sizeof(int));
  int *MRABsendBuffer=NULL;
  for(int lev=0; lev<MRABNlevels; lev++){
    dfloat dtlev = dtGmin*pow(2,lev);
    for(dlong e=0;e<Nelements;++e){
      if(EToDT[e] >=dtlev)
        MRABlevel[e] = lev;
    }
  }

  //enforce one level difference between neighbours
  if (totalHaloPairs)
    MRABsendBuffer = (int *) calloc(totalHaloPairs,sizeof(int));

  for (int lev=0; lev < MRABNlevels; lev++){
    if (totalHaloPairs)
      this->HaloExchange(sizeof(int), MRABlevel, MRABsendBuffer, MRABlevel+Nelements);
    for (dlong e =0; e<Nelements;e++) {
      if (MRABlevel[e] > lev+1) { //find elements at least 2 levels higher than lev
        for (int f=0;f<Nfaces;f++) { //check for a level lev neighbour
          int eP = EToE[Nfaces*e+f];
          if (eP > -1)
            if (MRABlevel[eP] == lev)
              MRABlevel[e] = lev + 1;  //if one exists, lower the level of this element to lev-1
        }
      }
    }
  }



  if (totalHaloPairs) free(MRABsendBuffer);

  //this could change the number of MRAB levels there are, so find the new max level
  MRABNlevels = 0;
  for (dlong e=0;e<Nelements;e++)
    MRABNlevels = (MRABlevel[e]>MRABNlevels) ? MRABlevel[e] : MRABNlevels;
  MRABNlevels++;
  int localNlevels = MRABNlevels;
  MPI_Allreduce(&localNlevels, &(MRABNlevels), 1, MPI_INT, MPI_MAX, comm);
  // NtimeSteps = finalTime/(pow(2,MRABNlevels-1)*dtGmin);

  printf("MRABNlevels %d \n", MRABNlevels);

  //now we need to perform a weighted repartitioning of the mesh to optimize MRAB
  if (size>1) {
    //for the moment, just weigth the elements by the number or RHS evals per MRAB step
    // TODO: We should make this an input parameter later to handle other problems.
    dfloat *weights = (dfloat *) calloc(Nelements,sizeof(dfloat));
    for (dlong e=0; e<Nelements;e++) {
      weights[e] = pow(2,MRABNlevels-MRABlevel[e]);
    }

    if (rank==0) printf("Repartitioning for MRAB...\n");
    this->MRABWeightedPartition(weights,MRABNlevels, MRABlevel);
  }

  //construct element and halo lists
  MRABelementIds = (int **) calloc(MRABNlevels,sizeof(int*));
  MRABhaloIds = (int **) calloc(MRABNlevels,sizeof(int*));

  MRABNelements = (int *) calloc(MRABNlevels,sizeof(int));
  MRABNhaloElements = (int *) calloc(MRABNlevels,sizeof(int));

  for (dlong e=0;e<Nelements;e++) {
    MRABNelements[MRABlevel[e]]++;
    for (int f=0;f<Nfaces;f++) {
      int eP = EToE[Nfaces*e+f];
      if (eP > -1) {
        if (MRABlevel[eP] == MRABlevel[e]-1) {//check for a level lev-1 neighbour
          MRABNhaloElements[MRABlevel[e]]++;
          break;
        }
      }
    }
  }

  for (int lev =0;lev<MRABNlevels;lev++){
    MRABelementIds[lev] = (int *) calloc(MRABNelements[lev],sizeof(int));
    MRABhaloIds[lev] = (int *) calloc(MRABNhaloElements[lev],sizeof(int));
    int cnt  =0;
    int cnt2 =0;
    for (dlong e=0;e<Nelements;e++){
      if (MRABlevel[e] == lev) {
        MRABelementIds[lev][cnt++] = e;

        for (int f=0;f<Nfaces;f++) {
          dlong eP = EToE[Nfaces*e+f];
          if (eP > -1) {
            if (MRABlevel[eP] == lev-1) {//check for a level lev-1 neighbour
              MRABhaloIds[lev][cnt2++] = e;
              break;
            }
          }
        }
      }
    }
  }

  //offset index
  MRABshiftIndex = (int *) calloc(MRABNlevels,sizeof(int));

  if (rank==0){
    printf("| Rank | Level | Nelements | Level/Level Boundary Elements | \n");
    printf("------------------------------------------------------------\n");
  }
  MPI_Barrier(comm);
  for (int rr =0;rr<size;rr++) {
    if (rr==rank) {
      for (int lev =0; lev<MRABNlevels; lev++)
        printf("|  %d,    %d,      %d,        %d     \n", rank, lev, MRABNelements[lev], MRABNhaloElements[lev]);
      printf("------------------------------------------------------------\n");
    }
    MPI_Barrier(comm);
  }
  MPI_Barrier(comm);


  return dtGmin;
}



#if 0

void meshMRABSetup2D(mesh2D *mesh, dfloat *EToDT, int maxLevels) {

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  //find global min and max dt
  dfloat dtmin, dtmax;
  dtmin = EToDT[0];
  dtmax = EToDT[0];
  for (int e=1;e<Nelements;e++) {
    dtmin = mymin(dtmin,EToDT[e]);
    dtmax = mymax(dtmax,EToDT[e]);
  }
  dfloat dtGmin, dtGmax;
  MPI_Allreduce(&dtmin, &dtGmin, 1, MPI_DFLOAT, MPI_MIN, comm);
  MPI_Allreduce(&dtmax, &dtGmax, 1, MPI_DFLOAT, MPI_MIN, comm);


  if (rank==0) {
    printf("----------- MRAB Setup ------------------------------\n");
    printf("dtmin = %g, dtmax = %g\n", dtGmin, dtGmax);
  }

  //number of levels
  MRABNlevels = mymin(floor(log2(dtGmax/dtGmin))+1,maxLevels);

  //shift dtGmin so that we have an integer number of steps
  NtimeSteps = finalTime/(pow(2,MRABNlevels-1)*dtGmin);
  dtGmin = finalTime/(pow(2,MRABNlevels-1)*NtimeSteps);

  dt = dtGmin;

  //compute the level of each element
  MRABlevel = (int *) calloc(Nelements+totalHaloPairs,sizeof(int));
  int *MRABsendBuffer;
  for(int lev=0; lev<MRABNlevels; lev++){
    dfloat dtlev = dtGmin*pow(2,lev);
    for(int e=0;e<Nelements;++e){
      if(EToDT[e] >=dtlev)
        MRABlevel[e] = lev;
    }
  }

  //enforce one level difference between neighbours
  if (totalHaloPairs)
    MRABsendBuffer = (int *) calloc(totalHaloPairs,sizeof(int));

  for (int lev=0; lev < MRABNlevels; lev++){
    if (totalHaloPairs)
      meshHaloExchange(mesh, sizeof(int), MRABlevel, MRABsendBuffer, MRABlevel+Nelements);
    for (int e =0; e<Nelements;e++) {
      if (MRABlevel[e] > lev+1) { //find elements at least 2 levels higher than lev
        for (int f=0;f<Nfaces;f++) { //check for a level lev neighbour
          int eP = EToE[Nfaces*e+f];
          if (eP > -1)
            if (MRABlevel[eP] == lev)
              MRABlevel[e] = lev + 1;  //if one exists, lower the level of this element to lev-1
        }
      }
    }
  }



  if (totalHaloPairs) free(MRABsendBuffer);

  //this could change the number of MRAB levels there are, so find the new max level
  MRABNlevels = 0;
  for (int e=0;e<Nelements;e++)
    MRABNlevels = (MRABlevel[e]>MRABNlevels) ? MRABlevel[e] : MRABNlevels;
  MRABNlevels++;
  int localNlevels = MRABNlevels;
  MPI_Allreduce(&localNlevels, &(MRABNlevels), 1, MPI_INT, MPI_MAX, comm);
  NtimeSteps = finalTime/(pow(2,MRABNlevels-1)*dtGmin);

  printf("MRABNlevels %d \n", MRABNlevels);

  //now we need to perform a weighted repartitioning of the mesh to optimize MRAB
  if (size>1) {
    //for the moment, just weigth the elements by the number or RHS evals per MRAB step
    // TODO: We should make this an input parameter later to handle other problems.
    dfloat *weights = (dfloat *) calloc(Nelements,sizeof(dfloat));
    for (int e=0; e<Nelements;e++) {
      weights[e] = pow(2,MRABNlevels-MRABlevel[e]);
    }

    if (rank==0) printf("Repartitioning for MRAB...\n");
    meshMRABWeightedPartitionTri2D(mesh,weights,MRABNlevels, MRABlevel);
  }

  //construct element and halo lists
  MRABelementIds = (int **) calloc(MRABNlevels,sizeof(int*));
  MRABhaloIds = (int **) calloc(MRABNlevels,sizeof(int*));

  MRABNelements = (int *) calloc(MRABNlevels,sizeof(int));
  MRABNhaloElements = (int *) calloc(MRABNlevels,sizeof(int));

  for (int e=0;e<Nelements;e++) {
    MRABNelements[MRABlevel[e]]++;
    for (int f=0;f<Nfaces;f++) {
      int eP = EToE[Nfaces*e+f];
      if (eP > -1) {
        if (MRABlevel[eP] == MRABlevel[e]-1) {//check for a level lev-1 neighbour
          MRABNhaloElements[MRABlevel[e]]++;
          break;
        }
      }
    }
  }

  for (int lev =0;lev<MRABNlevels;lev++){
    MRABelementIds[lev] = (int *) calloc(MRABNelements[lev],sizeof(int));
    MRABhaloIds[lev] = (int *) calloc(MRABNhaloElements[lev],sizeof(int));
    int cnt  =0;
    int cnt2 =0;
    for (int e=0;e<Nelements;e++){
      if (MRABlevel[e] == lev) {
        MRABelementIds[lev][cnt++] = e;

        for (int f=0;f<Nfaces;f++) {
          int eP = EToE[Nfaces*e+f];
          if (eP > -1) {
            if (MRABlevel[eP] == lev-1) {//check for a level lev-1 neighbour
              MRABhaloIds[lev][cnt2++] = e;
              break;
            }
          }
        }
      }
    }
  }



  //offset index
  MRABshiftIndex = (int *) calloc(MRABNlevels,sizeof(int));

  if (rank==0){
    printf("| Rank | Level | Nelements | Level/Level Boundary Elements | \n");
    printf("------------------------------------------------------------\n");
  }
  MPI_Barrier(comm);
  for (int rr =0;rr<size;rr++) {
    if (rr==rank) {
      for (int lev =0; lev<MRABNlevels; lev++)
        printf("|  %d,    %d,      %d,        %d     \n", rank, lev, MRABNelements[lev], MRABNhaloElements[lev]);
      printf("------------------------------------------------------------\n");
    }
    MPI_Barrier(comm);
  }
  MPI_Barrier(comm);
}
#endif
