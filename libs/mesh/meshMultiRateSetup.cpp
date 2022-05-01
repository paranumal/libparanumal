/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

namespace libp {

void mesh_t::MultiRateSetup(memory<dfloat> EToDT) {

  const int maxLevels = 100;

  //find global min and max dt
  dfloat dtmin = std::numeric_limits<dfloat>::max();
  dfloat dtmax = std::numeric_limits<dfloat>::min();
  for (dlong e=0;e<Nelements;e++) {
    dtmin = std::min(dtmin,EToDT[e]);
    dtmax = std::max(dtmax,EToDT[e]);
  }
  comm.Allreduce(dtmin, Comm::Min);
  comm.Allreduce(dtmax, Comm::Max);

  if (rank==0) {
    printf("--------------- MultiRate Timestepping Setup ----------------\n");
    printf("-------------------------------------------------------------\n");
  }

  //number of levels
  mrNlevels = std::min(static_cast <int>(std::floor(std::log2(dtmax/dtmin)))+1,
                                         maxLevels);

  //compute the level of each element
  mrLevel.malloc(Nelements+totalHaloPairs);
  for(int lev=0; lev<mrNlevels; lev++){
    dfloat dtlev = dtmin*(1<<lev);
    for(dlong e=0;e<Nelements;++e){
      if(EToDT[e] >=dtlev) {
        mrLevel[e] = lev;
      }
    }
  }

  //enforce one level difference between neighbours
  for (int lev=0; lev < mrNlevels; lev++){

    halo.Exchange(mrLevel, 1);

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
    mrNlevels = std::max(mrLevel[e],mrNlevels);
  mrNlevels++;

  comm.Allreduce(mrNlevels, Comm::Max);

  //construct element and halo lists
  // mrElements[lev] - list of all elements with multirate level <= lev
  // mrInterfaceElements[lev] - list of all elements with multirate level = lev,
  //                                with a neighbor of level lev-1
  mrNelements.malloc(mrNlevels, 0);
  mrInterfaceNelements.malloc(mrNlevels, 0);

  mrElements.malloc(mrNlevels);
  mrInterfaceElements.malloc(mrNlevels);

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
    mrElements[lev].malloc(mrNelements[lev]);
    mrInterfaceElements[lev].malloc(mrInterfaceNelements[lev]);
  }

  memory<int> cnt(mrNlevels, 0);
  memory<int> cnt2(mrNlevels, 0);

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

  o_mrLevel = platform.malloc<int>(Nelements, mrLevel);
  o_mrNelements = platform.malloc<dlong>(mrNlevels, mrNelements);
  o_mrInterfaceNelements = platform.malloc<dlong>(mrNlevels, mrInterfaceNelements);

  o_mrElements.malloc(mrNlevels);
  o_mrInterfaceElements.malloc(mrNlevels);

  for (int lev =0;lev<mrNlevels;lev++){
    if (mrNelements[lev])
      o_mrElements[lev]          = platform.malloc<dlong>(mrNelements[lev], mrElements[lev]);
    if (mrInterfaceNelements[lev])
      o_mrInterfaceElements[lev] = platform.malloc<dlong>(mrInterfaceNelements[lev], mrInterfaceElements[lev]);
  }

  if (rank==0){
    printf("| Level |     Nelements     | Level/Level Boundary Elements | \n");
    printf("-------------------------------------------------------------\n");
  }

  hlong Ntotal=0;
  for (int lev=0; lev<mrNlevels; lev++) {

    hlong levNelements = mrNelements[lev];
    comm.Allreduce(levNelements);
    levNelements -= Ntotal;
    Ntotal += levNelements;

    dlong levInterfaceNelements = mrInterfaceNelements[lev];
    comm.Allreduce(levInterfaceNelements);

    if (rank==0)
      printf("|   %3d |      %12lu |                 %12lu  |\n", lev, (size_t)levNelements, (size_t)levInterfaceNelements);
  }
  if (rank==0)
    printf("-------------------------------------------------------------\n");
}

} //namespace libp
