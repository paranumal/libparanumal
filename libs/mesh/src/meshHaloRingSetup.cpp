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

typedef struct{

  hlong gid;
  dlong element;
  dlong rank;
  dlong dest;

}vertex_t;

// comparator on destination rank
static int compareDest(const void *a,
                       const void *b){

  vertex_t *va = (vertex_t*) a;
  vertex_t *vb = (vertex_t*) b;

  if(va->dest < vb->dest) return -1;
  if(va->dest > vb->dest) return +1;

  return 0;
}

// comparator on global id
static int compareGlobalId(const void *a,
                           const void *b){

  vertex_t *va = (vertex_t*) a;
  vertex_t *vb = (vertex_t*) b;

  if(va->gid < vb->gid) return -1;
  if(va->gid > vb->gid) return +1;

  return 0;
}

// comparator on rank and element
static int compareRank(const void *a,
                       const void *b){

  vertex_t *va = (vertex_t*) a;
  vertex_t *vb = (vertex_t*) b;

  if(va->rank < vb->rank) return -1;
  if(va->rank > vb->rank) return +1;

  if(va->element < vb->element) return -1;
  if(va->element > vb->element) return +1;

  return 0;
}

// set up halo exchange for entire boundary ring elements for inter-processor MPI
// exchange of trace nodes
void mesh_t::HaloRingSetup(){

  //make a global indexing of element Ids
  hlong *globalOffsets = (hlong *) calloc(size+1,sizeof(hlong));
  hlong localNelements = (hlong) Nelements;

  //gather number of elements on each rank
  MPI_Allgather(&localNelements, 1, MPI_HLONG, globalOffsets+1, 1, MPI_HLONG, comm);

  for(int rr=0;rr<size;++rr)
    globalOffsets[rr+1] = globalOffsets[rr]+globalOffsets[rr+1];

  //use the gs to find what nodes are local to this rank
  dlong Ntotal = Np*Nelements;
  int *minRank = (int *) calloc(Ntotal,sizeof(int));
  int *maxRank = (int *) calloc(Ntotal,sizeof(int));
  for (dlong i=0;i<Ntotal;i++) {
    minRank[i] = rank;
    maxRank[i] = rank;
  }

  ogs->GatherScatter(minRank, ogs_int, ogs_min, ogs_sym); //minRank[n] contains the smallest rank taking part in the gather of node n
  ogs->GatherScatter(maxRank, ogs_int, ogs_max, ogs_sym); //maxRank[n] contains the largest rank taking part in the gather of node n

  //We already made a list of the globally connected element in ParallelGatherScatterSetup
  // NglobalGatherElements and globalGatherElementList contain the count and list

  //Make a list of the elements participating in the ring exchange
  //Count the number of shared vertices in the local mesh
  dlong NsendVerts=0;
  for (int e=0;e<NglobalGatherElements;e++) { //for all global elements
    for (int v=0;v<Nverts;v++) {
      dlong n = vertexNodes[v] + globalGatherElementList[e]*Np; //Id of a vertex in a global element
      if ((minRank[n]!=rank)||(maxRank[n]!=rank)) { //vertex is shared
        NsendVerts++;
      }
    }
  }

  vertex_t *vertexSendList = (vertex_t*) malloc(NsendVerts*sizeof(vertex_t));

  int *vertexSendCounts = (int*) calloc(size, sizeof(int));
  int *vertexRecvCounts = (int*) calloc(size, sizeof(int));
  int *vertexSendOffsets = (int*) calloc(size+1, sizeof(int));
  int *vertexRecvOffsets = (int*) calloc(size+1, sizeof(int));

  // Make the MPI_VERTEX_T data type
  MPI_Datatype MPI_VERTEX_T;
  MPI_Datatype dtype[4] = {MPI_HLONG, MPI_DLONG, MPI_DLONG, MPI_DLONG};
  int blength[4] = {1, 1, 1, 1};
  MPI_Aint addr[4], displ[4];
  MPI_Get_address ( &(vertexSendList[0]        ), addr+0);
  MPI_Get_address ( &(vertexSendList[0].element), addr+1);
  MPI_Get_address ( &(vertexSendList[0].rank   ), addr+2);
  MPI_Get_address ( &(vertexSendList[0].dest   ), addr+3);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  displ[3] = addr[3] - addr[0];
  MPI_Type_create_struct (4, blength, displ, dtype, &MPI_VERTEX_T);
  MPI_Type_commit (&MPI_VERTEX_T);

  NsendVerts=0;
  for (int e=0;e<NglobalGatherElements;e++) { //for all global elements
    for (int v=0;v<Nverts;v++) {
      dlong n = vertexNodes[v] + globalGatherElementList[e]*Np; //Id of a vertex in a global element
      if ((minRank[n]!=rank)||(maxRank[n]!=rank)) { //vertex is shared
        vertexSendList[NsendVerts].gid = globalIds[n]; //global node index
        vertexSendList[NsendVerts].element = globalGatherElementList[e]; //local element index
        vertexSendList[NsendVerts].rank = rank;
        vertexSendList[NsendVerts].dest = globalIds[n]%size; //destination rank for sorting

        vertexSendCounts[vertexSendList[NsendVerts].dest]++; //count outgoing

        NsendVerts++;
      }
    }
  }

  free(minRank); free(maxRank);

  // sort based on destination (=gid%size)
  qsort(vertexSendList, NsendVerts, sizeof(vertex_t), compareDest);

  // share counts
  MPI_Alltoall(vertexSendCounts, 1, MPI_INT,
               vertexRecvCounts, 1, MPI_INT,
               comm);

  dlong NrecvVerts = 0;
  for(int rr=0;rr<size;++rr){
    NrecvVerts += vertexRecvCounts[rr];

    vertexSendOffsets[rr+1] = vertexSendOffsets[rr] + vertexSendCounts[rr];
    vertexRecvOffsets[rr+1] = vertexRecvOffsets[rr] + vertexRecvCounts[rr];
  }

  vertex_t *vertexRecvList = (vertex_t*) malloc(NrecvVerts*sizeof(vertex_t));

  // exchange shared vertices
  MPI_Alltoallv(vertexSendList, vertexSendCounts, vertexSendOffsets, MPI_VERTEX_T,
                vertexRecvList, vertexRecvCounts, vertexRecvOffsets, MPI_VERTEX_T,
                comm);

  // sort based on globalId to find matches
  qsort(vertexRecvList, NrecvVerts, sizeof(vertex_t), compareGlobalId);

  //count the number of unique vertices we have
  dlong Nunique=(NrecvVerts) ? 1:0; //at least one if there's any
  for(dlong n=1;n<NrecvVerts;++n){
    if (vertexRecvList[n].gid != vertexRecvList[n-1].gid) { // new vertex
      Nunique++;
    }
  }

  //Build offsets to unique vertice starts
  dlong *vertexOffsets = (dlong*) calloc(Nunique+1, sizeof(dlong));

  Nunique=(NrecvVerts) ? 1:0;
  for(dlong n=1;n<NrecvVerts;++n){
    if (vertexRecvList[n].gid != vertexRecvList[n-1].gid) { // new vertex
      vertexOffsets[Nunique++] = n;
    }
  }
  vertexOffsets[Nunique] = NrecvVerts;

  //make sure the AlltoAll is done everywhere so we can reuse vertexSend arrays
  MPI_Barrier(comm);

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
  vertexSendList = (vertex_t*) realloc(vertexSendList, NsendVerts*sizeof(vertex_t));

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
  qsort(vertexSendList, NsendVerts, sizeof(vertex_t), compareDest);

  // share counts
  MPI_Alltoall(vertexSendCounts, 1, MPI_INT,
               vertexRecvCounts, 1, MPI_INT,
               comm);

  NrecvVerts = 0;
  for(int rr=0;rr<size;++rr){
    NrecvVerts += vertexRecvCounts[rr];

    vertexSendOffsets[rr+1] = vertexSendOffsets[rr] + vertexSendCounts[rr];
    vertexRecvOffsets[rr+1] = vertexRecvOffsets[rr] + vertexRecvCounts[rr];
  }

  //resize recv storage
  vertexRecvList = (vertex_t*) realloc(vertexRecvList,NrecvVerts*sizeof(vertex_t));

  // exchange shared vertices
  MPI_Alltoallv(vertexSendList, vertexSendCounts, vertexSendOffsets, MPI_VERTEX_T,
                vertexRecvList, vertexRecvCounts, vertexRecvOffsets, MPI_VERTEX_T,
                comm);

  // sort based on rank then element id to find matches
  qsort(vertexRecvList, NrecvVerts, sizeof(vertex_t), compareRank);

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
  hlong *globalElementId = (hlong *) malloc((Nelements+totalRingElements)*sizeof(hlong));

  //outgoing elements
  for(int e=0;e<Nelements;++e)
    globalElementId[e] = e + globalOffsets[rank] + 1;

  //incoming elements
  totalRingElements=0;

  //first one
  if (NrecvVerts>0) {
    if (vertexRecvList[0].rank!=rank) {
      globalElementId[Nelements]
                       = -(vertexRecvList[0].element
                          + globalOffsets[vertexRecvList[0].rank] + 1); //negative so doesnt contribute to sum in ogs
      totalRingElements++;
    }
  }

  for (dlong n=1;n<NrecvVerts;n++) {
    if (vertexRecvList[n].rank==rank) continue; //don't count the local elements

    if (  (vertexRecvList[n].rank!=vertexRecvList[n-1].rank)
        ||(vertexRecvList[n].element!=vertexRecvList[n-1].element)) {
      globalElementId[Nelements+totalRingElements++]
                       = -(vertexRecvList[n].element
                          + globalOffsets[vertexRecvList[n].rank] + 1); //negative so doesnt contribute to sum in ogs
    }
  }

  // finally we now have a list of all elements that we need to send to form
  //  the 1-ring (to rule them all)

  //make the halo exchange op
  int verbose = 0;
  ringHalo = halo_t::Setup(Nelements+totalRingElements, globalElementId, comm,
                           verbose, device);

  //clean up
  free(globalElementId);
  free(globalOffsets);

  MPI_Barrier(comm);
  MPI_Type_free(&MPI_VERTEX_T);
  free(vertexSendList);
  free(vertexRecvList);
  free(vertexSendCounts);
  free(vertexRecvCounts);
  free(vertexSendOffsets);
  free(vertexRecvOffsets);
}
