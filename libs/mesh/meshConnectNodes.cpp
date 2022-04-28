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

// uniquely label each node with a global index, used for gatherScatter
void mesh_t::ConnectNodes(){

  hlong localNnodes = Np*Nelements;
  hlong gatherNodeStart = localNnodes;
  comm.Scan(localNnodes, gatherNodeStart);
  gatherNodeStart -= localNnodes;

  // form global node numbering
  globalIds.malloc((totalHaloPairs+Nelements)*Np);

  // initialize with local numbering
  #pragma omp parallel for
  for(dlong n=0;n<Nelements*Np;++n){
    globalIds[n] = 1 + n + gatherNodeStart;
  }

  //make a node-wise bc flag by looking at all neighbors
  mapB.malloc((Nelements+totalHaloPairs)*Np, -1);

  #pragma omp parallel for
  for (dlong e=0;e<Nelements;e++) {
    for (int f=0;f<Nfaces;f++) {
      int bc = EToB[f+e*Nfaces];
      if (bc>0) {
        for (int n=0;n<Nfp;n++) {
          const int fid = faceNodes[n+f*Nfp];
          int bcn = mapB[fid+e*Np];
          if (bcn == -1) { //if theres no bc here yet, write it
            mapB[fid+e*Np] = bc;
          } else { //if theres a bc, take the min
            mapB[fid+e*Np] = std::min(bc,bcn);
          }
        }
      }
    }
  }

  hlong gatherChange = 1;

  // keep comparing numbers on positive and negative traces until convergence
  while(gatherChange>0){

    // reset change counter
    gatherChange = 0;

    // send halo data and recv into extension of buffer
    halo.Exchange(globalIds, Np);
    halo.Exchange(mapB, Np);

    // compare trace nodes
    // #pragma omp parallel for
    for(dlong e=0;e<Nelements;++e){

      for(int n=0;n<Nfp*Nfaces;++n){
        dlong id  = e*Nfp*Nfaces + n;
        dlong idM = vmapM[id];
        dlong idP = vmapP[id];
        hlong gidM = globalIds[idM];
        hlong gidP = globalIds[idP];
        int bcM = mapB[idM];
        int bcP = mapB[idP];

        if(gidP<gidM){
          ++gatherChange;
          globalIds[idM] = gidP;
        }

        if (bcP > 0) {
          if (bcM == -1) {
            //if theres no bc here yet, write it
            mapB[idM] = bcP;
            ++gatherChange;
          } else if (bcP<bcM) {
            mapB[idM] = bcP;
            ++gatherChange;
          }
        }
      }
    }

    // sum up changes
    comm.Allreduce(gatherChange);
  }

  o_mapB = platform.malloc<int>(mapB);
}

} //namespace libp
