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

#ifdef GLIBCXX_PARALLEL
#include <parallel/algorithm>
using __gnu_parallel::sort;
#else
using std::sort;
#endif

namespace libp {

typedef struct{

  hlong gid;
  dlong element;
  dlong rank;
  dlong dest;

}vertex_t;

// set up halo exchange for entire boundary ring elements for inter-processor MPI
// exchange of trace nodes
void mesh_t::HaloRingSetup(){

  memory<hlong> globalOffset(size+1, 0);

  //gather number of elements on each rank
  hlong localNelements = Nelements;
  comm.Allgather(localNelements, globalOffset+1);

  for(int rr=0;rr<size;++rr)
    globalOffset[rr+1] = globalOffset[rr]+globalOffset[rr+1];


  dlong Ntotal = Nverts*(Nelements+totalHaloPairs);

  memory<int> minRank(Ntotal);
  memory<int> maxRank(Ntotal);

  for (dlong i=0;i<Ntotal;i++) {
    minRank[i] = rank;
    maxRank[i] = rank;
  }

  hlong gatherChange = 1;

  // keep comparing numbers on positive and negative traces until convergence
  while(gatherChange>0){

    // reset change counter
    gatherChange = 0;

    // send halo data and recv into extension of buffer
    halo.Exchange(minRank, Nverts);
    halo.Exchange(maxRank, Nverts);

    // compare trace vertices
    #pragma omp parallel for collapse(2)
    for(dlong e=0;e<Nelements;++e){
      for(int n=0;n<Nfaces*NfaceVertices;++n){
        dlong id  = e*Nfaces*NfaceVertices + n;
        dlong idM = VmapM[id];
        dlong idP = VmapP[id];

        int minRankM = minRank[idM];
        int minRankP = minRank[idP];

        int maxRankM = maxRank[idM];
        int maxRankP = maxRank[idP];

        if(minRankP<minRankM){
          gatherChange=1;
          minRank[idM] = minRankP;
        }

        if(maxRankP>maxRankM){
          gatherChange=1;
          maxRank[idM] = maxRankP;
        }
      }
    }

    // sum up changes
    comm.Allreduce(gatherChange);
  }

  //Make a list of the elements participating in the ring exchange
  //Count the number of shared vertices in the local mesh
  dlong NsendVerts=0;
  for (int e=0;e<Nelements;e++) { //for all global elements
    for (int v=0;v<Nverts;v++) {
      dlong id = e*Nverts + v; //Id of a vertex in a global element
      if ((minRank[id]!=rank)||(maxRank[id]!=rank)) { //vertex is shared
        NsendVerts++;
      }
    }
  }

  memory<vertex_t> vertexSendList(NsendVerts);

  memory<int> vertexSendCounts(size, 0);
  memory<int> vertexRecvCounts(size);
  memory<int> vertexSendOffsets(size+1);
  memory<int> vertexRecvOffsets(size+1);

  NsendVerts=0;
  for (int e=0;e<Nelements;e++) { //for all elements
    for (int v=0;v<Nverts;v++) {
      dlong id = e*Nverts + v;
      if ((minRank[id]!=rank)||(maxRank[id]!=rank)) { //vertex is shared
        vertexSendList[NsendVerts].gid = EToV[id]; //global vertex index
        vertexSendList[NsendVerts].element = e; //local element index
        vertexSendList[NsendVerts].rank = rank;
        vertexSendList[NsendVerts].dest = EToV[id]%size; //destination rank for sorting

        vertexSendCounts[vertexSendList[NsendVerts].dest]++; //count outgoing

        NsendVerts++;
      }
    }
  }

  // sort based on destination (=gid%size)
  sort(vertexSendList.ptr(), vertexSendList.ptr()+NsendVerts,
            [](const vertex_t& a, const vertex_t& b)
              {return a.dest < b.dest;});

  // share counts
  comm.Alltoall(vertexSendCounts, vertexRecvCounts);

  dlong NrecvVerts = 0;
  vertexSendOffsets[0] = 0;
  vertexRecvOffsets[0] = 0;
  for(int rr=0;rr<size;++rr){
    NrecvVerts += vertexRecvCounts[rr];

    vertexSendOffsets[rr+1] = vertexSendOffsets[rr] + vertexSendCounts[rr];
    vertexRecvOffsets[rr+1] = vertexRecvOffsets[rr] + vertexRecvCounts[rr];
  }

  memory<vertex_t> vertexRecvList(NrecvVerts);

  // exchange shared vertices
  comm.Alltoallv(vertexSendList, vertexSendCounts, vertexSendOffsets,
                 vertexRecvList, vertexRecvCounts, vertexRecvOffsets);

  // sort based on globalId to find matches
  sort(vertexRecvList.ptr(), vertexRecvList.ptr()+NrecvVerts,
            [](const vertex_t& a, const vertex_t& b)
              {return a.gid < b.gid;});

  //count the number of unique vertices we have
  dlong Nunique=(NrecvVerts) ? 1:0; //at least one if there's any
  for(dlong n=1;n<NrecvVerts;++n){
    if (vertexRecvList[n].gid != vertexRecvList[n-1].gid) { // new vertex
      Nunique++;
    }
  }

  //Build offsets to unique vertice starts
  memory<dlong> vertexOffsets(Nunique+1);

  vertexOffsets[0] = 0;
  Nunique=(NrecvVerts) ? 1:0;
  for(dlong n=1;n<NrecvVerts;++n){
    if (vertexRecvList[n].gid != vertexRecvList[n-1].gid) { // new vertex
      vertexOffsets[Nunique++] = n;
    }
  }
  vertexOffsets[Nunique] = NrecvVerts;

  //reset counts
  NsendVerts = 0;
  for(int rr=0;rr<size;++rr){
    vertexSendCounts[rr] = 0;
    vertexSendOffsets[rr+1] = 0;
    vertexRecvOffsets[rr+1] = 0;
  }

  // now count how many vertices to send to each rank
  for(dlong n=0;n<Nunique;++n){
    const dlong start = vertexOffsets[n];
    const dlong end   = vertexOffsets[n+1];

    const int cnt = end - start; //number of elements contributing
    for(dlong m=start;m<end;++m){
      vertexSendCounts[vertexRecvList[m].rank] += cnt; // this could get large, but hopefully just O(NsharedVerts)
      NsendVerts += cnt; //total
    }
  }

  //resize send storage
  vertexSendList.malloc(NsendVerts);

  //build list of vertices to send out
  NsendVerts=0;
  for(dlong n=0;n<Nunique;++n){
    const dlong start = vertexOffsets[n];
    const dlong end   = vertexOffsets[n+1];

    for(dlong v1=start;v1<end;++v1){ // for each vertex of this globalId
      for(dlong v2=start;v2<end;++v2){//send info for all vertices
        vertexSendList[NsendVerts] = vertexRecvList[v2];
        vertexSendList[NsendVerts].dest = vertexRecvList[v1].rank; //send to v1's rank
        NsendVerts++;
      }
    }
  }

  // sort based on destination
  sort(vertexSendList.ptr(), vertexSendList.ptr()+NsendVerts,
            [](const vertex_t& a, const vertex_t& b)
              {return a.dest < b.dest;});

  // share counts
  comm.Alltoall(vertexSendCounts, vertexRecvCounts);

  NrecvVerts = 0;
  for(int rr=0;rr<size;++rr){
    NrecvVerts += vertexRecvCounts[rr];

    vertexSendOffsets[rr+1] = vertexSendOffsets[rr] + vertexSendCounts[rr];
    vertexRecvOffsets[rr+1] = vertexRecvOffsets[rr] + vertexRecvCounts[rr];
  }

  //resize recv storage
  vertexRecvList.malloc(NrecvVerts);

  // exchange shared vertices
  comm.Alltoallv(vertexSendList, vertexSendCounts, vertexSendOffsets,
                 vertexRecvList, vertexRecvCounts, vertexRecvOffsets);

  // sort based on rank then element id to find matches
  sort(vertexRecvList.ptr(), vertexRecvList.ptr()+NrecvVerts,
            [](const vertex_t& a, const vertex_t& b) {
              if(a.rank < b.rank) return true;
              if(a.rank > b.rank) return false;

              return a.element < b.element;
            });

  //count the number of elements in the 1-ring
  totalRingElements=(NrecvVerts) ? 1:0;
  for (dlong n=1;n<NrecvVerts;n++) {
    if (vertexRecvList[n].rank==rank) continue; //don't count the local elements

    if (  (vertexRecvList[n].rank!=vertexRecvList[n-1].rank)
        ||(vertexRecvList[n].element!=vertexRecvList[n-1].element)) {
      totalRingElements++;
    }
  }

  //make a list of global element ids taking part in the halo exchange
  memory<hlong> globalElementId(Nelements+totalRingElements);

  //outgoing elements
  for(int e=0;e<Nelements;++e)
    globalElementId[e] = e + globalOffset[rank] + 1;

  //incoming elements
  totalRingElements=0;

  //first one
  if (NrecvVerts>0) {
    if (vertexRecvList[0].rank!=rank) {
      globalElementId[Nelements]
                       = -(vertexRecvList[0].element
                          + globalOffset[vertexRecvList[0].rank] + 1); //negative so doesnt contribute to sum in ogs
      totalRingElements++;
    }
  }

  for (dlong n=1;n<NrecvVerts;n++) {
    if (vertexRecvList[n].rank==rank) continue; //don't count the local elements

    if (  (vertexRecvList[n].rank!=vertexRecvList[n-1].rank)
        ||(vertexRecvList[n].element!=vertexRecvList[n-1].element)) {
      globalElementId[Nelements+totalRingElements++]
                       = -(vertexRecvList[n].element
                          + globalOffset[vertexRecvList[n].rank] + 1); //negative so doesnt contribute to sum in ogs
    }
  }

  // finally we now have a list of all elements that we need to send to form
  //  the 1-ring (to rule them all)

  //make the halo exchange op
  int verbose = 0;
  ringHalo.Setup(Nelements+totalRingElements,
                 globalElementId, comm,
                 ogs::Auto, verbose, platform);
}

} //namespace libp
