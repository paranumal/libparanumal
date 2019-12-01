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

void mesh_t::MultiRateSetup(dfloat *EToDT) {

  const int maxLevels = 100;

  //find global min and max dt
  dfloat dtmin=1.e9, dtmax=0.0;
  if (Nelements) {
    dtmin = EToDT[0];
    dtmax = EToDT[0];
  }
  for (dlong e=1;e<Nelements;e++) {
    dtmin = mymin(dtmin,EToDT[e]);
    dtmax = mymax(dtmax,EToDT[e]);
  }
  dfloat dtGmin, dtGmax;
  MPI_Allreduce(&dtmin, &dtGmin, 1, MPI_DFLOAT, MPI_MIN, comm);
  MPI_Allreduce(&dtmax, &dtGmax, 1, MPI_DFLOAT, MPI_MIN, comm);

  if (rank==0) {
    printf("--------------- MultiRate Timestepping Setup ----------------\n");
    printf("-------------------------------------------------------------\n");
  }

  //number of levels
  mrNlevels = mymin(floor(log2(dtGmax/dtGmin))+1,maxLevels);

  //compute the level of each element
  mrLevel = (int *) calloc(Nelements+totalHaloPairs,sizeof(int));
  for(int lev=0; lev<mrNlevels; lev++){
    dfloat dtlev = dtGmin*(2<<lev);
    for(dlong e=0;e<Nelements;++e){
      if(EToDT[e] >=dtlev)
        mrLevel[e] = lev;
    }
  }

  //enforce one level difference between neighbours
  for (int lev=0; lev < mrNlevels; lev++){

    halo->Exchange(mrLevel, 1, ogs_int);

    for (dlong e=0; e<Nelements;e++) {
      if (mrLevel[e] > lev+1) { //find elements at least 2 levels higher than lev
        for (int f=0;f<Nfaces;f++) { //check for a level lev neighbour
          dlong eP = EToE[Nfaces*e+f];
          if (eP > -1)
            if (mrLevel[eP] == lev)
              mrLevel[e] = lev + 1;  //if one exists, lower the level of this element to lev-1
        }
      }
    }
  }

  //this could change the number of levels there are, so find the new max level
  mrNlevels = 0;
  for (dlong e=0;e<Nelements;e++)
    mrNlevels = (mrLevel[e]>mrNlevels) ? mrLevel[e] : mrNlevels;
  mrNlevels++;

  int localNlevels = mrNlevels;
  MPI_Allreduce(&localNlevels, &mrNlevels, 1, MPI_INT, MPI_MAX, comm);

  //construct element and halo lists
  // mrElements[lev] - list of all elements with multirate level <= lev
  // mrInterfaceElements[lev] - list of all elements with multirate level = lev,
  //                                with a neighbor of level lev-1
  mrNelements          = (dlong *) calloc(mrNlevels,sizeof(dlong));
  mrInterfaceNelements = (dlong *) calloc(mrNlevels,sizeof(dlong));

  mrElements          = (dlong **) calloc(mrNlevels,sizeof(dlong*));
  mrInterfaceElements = (dlong **) calloc(mrNlevels,sizeof(dlong*));

  for (dlong e=0;e<Nelements;e++) {
    int lev = mrLevel[e];
    for (int l=lev;l<mrNlevels;l++) mrNelements[l]++;

    //check neighbors
    for (int f=0;f<Nfaces;f++) {
      dlong eP = EToE[Nfaces*e+f];
      if (eP > -1) {
        if (mrLevel[eP] == lev-1) {//check for a level lev-1 neighbour
          mrInterfaceNelements[lev]++;
          break;
        }
      }
    }
  }

  //allocate space
  for (int lev =0;lev<mrNlevels;lev++){
    mrElements[lev]          = (dlong *) calloc(mrNelements[lev],sizeof(dlong));
    mrInterfaceElements[lev] = (dlong *) calloc(mrInterfaceNelements[lev],sizeof(dlong));
  }

  int *cnt  = (int *) calloc(mrNlevels,sizeof(int));
  int *cnt2 = (int *) calloc(mrNlevels,sizeof(int));

  //fill element lists
  for (dlong e=0;e<Nelements;e++){
    int lev = mrLevel[e];
    for (int l=lev;l<mrNlevels;l++)
      mrElements[l][cnt[l]++] = e;

    for (int f=0;f<Nfaces;f++) {
      dlong eP = EToE[Nfaces*e+f];
      if (eP > -1) {
        if (mrLevel[eP] == lev-1) {//check for a level lev-1 neighbour
          mrInterfaceElements[lev][cnt2[lev]++] = e;
          break;
        }
      }
    }
  }
  free(cnt); free(cnt2);

  o_mrLevel = device.malloc(Nelements*sizeof(int), mrLevel);
  o_mrNelements = device.malloc(mrNlevels*sizeof(dlong), mrNelements);
  o_mrInterfaceNelements = device.malloc(mrNlevels*sizeof(dlong), mrInterfaceNelements);

  o_mrElements          = new occa::memory[mrNlevels];
  o_mrInterfaceElements = new occa::memory[mrNlevels];

  for (int lev =0;lev<mrNlevels;lev++){
    if (mrNelements[lev])
      o_mrElements[lev]          = device.malloc(mrNelements[lev]*sizeof(dlong), mrElements[lev]);
    if (mrInterfaceNelements[lev])
      o_mrInterfaceElements[lev] = device.malloc(mrInterfaceNelements[lev]*sizeof(dlong), mrInterfaceElements[lev]);
  }

  if (rank==0){
    printf("| Level |     Nelements     | Level/Level Boundary Elements | \n");
    printf("-------------------------------------------------------------\n");
  }

  hlong Ntotal=0;
  for (int lev=0; lev<mrNlevels; lev++) {

    hlong levNelementsLocal = mrNelements[lev];
    hlong levNelements=0;
    MPI_Allreduce(&levNelementsLocal, &levNelements, 1, MPI_HLONG, MPI_SUM, comm);
    levNelements -= Ntotal;
    Ntotal += levNelements;

    dlong levInterfaceNelementsLocal = mrInterfaceNelements[lev];
    dlong levInterfaceNelements=0;
    MPI_Allreduce(&levInterfaceNelementsLocal, &levInterfaceNelements, 1, MPI_DLONG, MPI_SUM, comm);

    if (rank==0)
      printf("|   %3d |      %12lu |                 %12lu  |\n", lev, (size_t)levNelements, (size_t)levInterfaceNelements);
  }
  if (rank==0)
    printf("-------------------------------------------------------------\n");
}
