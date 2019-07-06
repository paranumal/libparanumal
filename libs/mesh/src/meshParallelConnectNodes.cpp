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


// uniquely label each node with a global index, used for gatherScatter
void mesh_t::ParallelConnectNodes(){

  dlong localNodeCount = Np*Nelements;
  dlong *allLocalNodeCounts = (dlong*) calloc(size, sizeof(dlong));

  MPI_Allgather(&localNodeCount,    1, MPI_DLONG,
                allLocalNodeCounts, 1, MPI_DLONG,
                comm);

  hlong gatherNodeStart = 0;
  for(int rr=0;rr<rank;++rr)
    gatherNodeStart += allLocalNodeCounts[rr];

  free(allLocalNodeCounts);

  // form continuous node numbering (local=>virtual gather)
  int *baseRank = (int *) malloc((totalHaloPairs+Nelements)*Np*sizeof(int));
  globalIds = (hlong *) malloc((totalHaloPairs+Nelements)*Np*sizeof(hlong));

  // use local numbering
  for(dlong e=0;e<Nelements;++e){
    for(int n=0;n<Np;++n){
      dlong id = e*Np+n;

      baseRank[id] = rank;
      globalIds[id] = 1 + id + Nnodes + gatherNodeStart;
    }

    // use vertex ids for vertex nodes to reduce iterations
    for(int v=0;v<Nverts;++v){
      dlong id = e*Np + vertexNodes[v];
      hlong gid = EToV[e*Nverts+v] + 1;
      globalIds[id] = gid;
    }
  }

  dlong localChange = 0, gatherChange = 1;

  // keep comparing numbers on positive and negative traces until convergence
  while(gatherChange>0){

    // reset change counter
    localChange = 0;

    // send halo data and recv into extension of buffer
    halo->Exchange(baseRank, Np, ogs_int);
    halo->Exchange(globalIds, Np, ogs_hlong);

    // compare trace nodes
    for(dlong e=0;e<Nelements;++e){
      for(int n=0;n<Nfp*Nfaces;++n){
        dlong id  = e*Nfp*Nfaces + n;
        dlong idM = vmapM[id];
        dlong idP = vmapP[id];
        hlong gidM = globalIds[idM];
        hlong gidP = globalIds[idP];

        int baseRankM = baseRank[idM];
        int baseRankP = baseRank[idP];

        if(gidM<gidP || (gidP==gidM && baseRankM<baseRankP)){
          ++localChange;
          baseRank[idP]  = baseRank[idM];
          globalIds[idP] = globalIds[idM];
        }

        if(gidP<gidM || (gidP==gidM && baseRankP<baseRankM)){
          ++localChange;
          baseRank[idM]  = baseRank[idP];
          globalIds[idM] = globalIds[idP];
        }
      }
    }

    // sum up changes
    MPI_Allreduce(&localChange, &gatherChange, 1, MPI_DLONG, MPI_SUM, comm);
  }

  free(baseRank);
}
