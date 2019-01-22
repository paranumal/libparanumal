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

#include "elliptic.h"

// 1. create: vertex, element, rank triples
// 2. send to process vertex%size

// 3. on each process sort by vertex, then rank, then element

// 4. for each vertex, send the list back to ranks in the list
// 5. then each rank sorts the list by source-rank and eliminates duplicates and self connections

typedef struct{

  hlong vertex;
  hlong element;
  hlong sourceRank;
  hlong sortRank;
  hlong otherRank;
}vertex_t;


int compareSortRank(const void *a, 
		    const void *b){

  vertex_t *va = (vertex_t*) a;
  vertex_t *vb = (vertex_t*) b;

  if(va->sortRank < vb->sortRank) return -1;
  if(va->sortRank > vb->sortRank) return +1;

  return 0;
}

int compareVertex(const void *a, 
		  const void *b){
  
  vertex_t *va = (vertex_t*) a;
  vertex_t *vb = (vertex_t*) b;
  
  if(va->vertex < vb->vertex) return -1;
  if(va->vertex > vb->vertex) return +1;

  return 0;
}

int compareSourceRank(const void *a, 
		    const void *b){

  vertex_t *va = (vertex_t*) a;
  vertex_t *vb = (vertex_t*) b;

  if(va->sourceRank < vb->sourceRank) return -1;
  if(va->sourceRank > vb->sourceRank) return +1;

  return 0;
}

int compareOtherRankElement(const void *a, 
			    const void *b){
  
  vertex_t *va = (vertex_t*) a;
  vertex_t *vb = (vertex_t*) b;

  if(va->otherRank < vb->otherRank) return -1;
  if(va->otherRank > vb->otherRank) return +1;

  if(va->element < vb->element) return -1;
  if(va->element > vb->element) return +1;

  
  return 0;
}

void ellipticBuildOneRing(elliptic_t *elliptic){

  mesh_t *mesh = elliptic->mesh;

  vertex_t *vertexSendList = (vertex_t*) calloc(mesh->Nelements*mesh->Nverts, sizeof(vertex_t));

  hlong *vertexSendCounts = (hlong*) calloc(mesh->size, sizeof(hlong));
  hlong *vertexRecvCounts = (hlong*) calloc(mesh->size, sizeof(hlong));

  hlong cnt = 0;
  for(hlong e=0;e<mesh->Nelements;++e){
    for(int v=0;v<mesh->Nverts;++v){
      vertexSendList[cnt].vertex = mesh->EToV[e*mesh->Nverts+v];
      vertexSendList[cnt].element = e;
      vertexSendList[cnt].sourceRank = mesh->rank;
      vertexSendList[cnt].sortRank = vertexSendList[cnt].vertex%mesh->size;
      ++vertexSendCounts[vertexSendList[cnt].sortRank];
      ++cnt;
    }
  }
  
  // sort based on sortRank
  qsort(vertexSendList, cnt, sizeof(vertex_t), compareSortRank);

#if 0
  for(int v=0;v<cnt;++v){
    printf("rank: %d, vertex: %d, element: %d, sourceRank: %d, sortRank: %d\n",
	   mesh->rank,
	   vertexSendList[v].vertex,
	   vertexSendList[v].element,
	   vertexSendList[v].sourceRank,
	   vertexSendList[v].sortRank);
  }
#endif
  
  // send sortRankCounts (hackety)
  MPI_Alltoall(vertexSendCounts, 1, MPI_HLONG,
	       vertexRecvCounts, 1, MPI_HLONG,
	       mesh->comm);

  // exchange vertices
  hlong *vertexSendDispls = (hlong*) calloc(mesh->size+1, sizeof(hlong));
  hlong *vertexRecvDispls = (hlong*) calloc(mesh->size+1, sizeof(hlong));
  hlong NvertexSend = 0;
  hlong NvertexRecv = 0;
  for(int r=0;r<mesh->size;++r){

    NvertexSend += vertexSendCounts[r];
    NvertexRecv += vertexRecvCounts[r];
    
    vertexSendCounts[r] *= sizeof(vertex_t); // hack-hack-hack
    vertexRecvCounts[r] *= sizeof(vertex_t); // hack-hack-hack
    
    vertexSendDispls[r+1] = vertexSendDispls[r] + vertexSendCounts[r];
    vertexRecvDispls[r+1] = vertexRecvDispls[r] + vertexRecvCounts[r];

    printf("rank %d: send %d, %d  recv %d, %d\n",
	   mesh->rank,
	   vertexSendCounts[r], vertexSendDispls[r],
	   vertexRecvCounts[r], vertexRecvDispls[r]);
  }

  vertex_t *vertexRecvList = (vertex_t*) calloc(NvertexRecv, sizeof(vertex_t)); // hack-hack-hack
  
  MPI_Alltoallv(vertexSendList, vertexSendCounts, vertexSendDispls, MPI_CHAR,
		vertexRecvList, vertexRecvCounts, vertexRecvDispls, MPI_CHAR,
		mesh->comm);

  // sort received vertex based on vertex number
  qsort(vertexRecvList, NvertexRecv, sizeof(vertex_t), compareVertex);  

  // count number of unique received vertices
  hlong NvertexUniqueRecv = (NvertexRecv>0) ? 1:0;
  for(hlong n=1;n<NvertexRecv;++n){
    if(compareVertex(vertexRecvList+n, vertexRecvList+n-1)!=0) { // new vertex
      ++NvertexUniqueRecv;
    }
  }

  // find offset of the start of each new unique vertex  in sorted list
  hlong *vertexUniqueRecvDispls = (hlong*) calloc(NvertexUniqueRecv+1, sizeof(hlong));

  cnt = 1;
  vertexUniqueRecvDispls[0] = 0;
  for(hlong n=1;n<NvertexRecv;++n){
    if(compareVertex(vertexRecvList+n, vertexRecvList+n-1)!=0) { // new vertex
      vertexUniqueRecvDispls[cnt] = n;
      ++cnt;
    }
  }
  vertexUniqueRecvDispls[cnt] = NvertexRecv; // cap at end
  
  // now count how many vertices to send to each rank
  hlong *vertexHaloSendCounts = (hlong*) calloc(mesh->size, sizeof(hlong));
  hlong Ntotal = 0;
  for(hlong n=0;n<NvertexUniqueRecv;++n){
    int NuniqueRecvMultiplicity = vertexUniqueRecvDispls[n+1] - vertexUniqueRecvDispls[n];
    for(hlong m=vertexUniqueRecvDispls[n];m<vertexUniqueRecvDispls[n+1];++m){
      vertexHaloSendCounts[vertexRecvList[m].sourceRank] += NuniqueRecvMultiplicity; // watch out for this
      Ntotal += NuniqueRecvMultiplicity;
    }
  }

  vertex_t *vertexHaloSendList = (vertex_t*) calloc(Ntotal, sizeof(vertex_t));
  cnt = 0;
  for(hlong n=0;n<NvertexUniqueRecv;++n){
    int NuniqueRecvMultiplicity = vertexUniqueRecvDispls[n+1] - vertexUniqueRecvDispls[n];
    hlong start = vertexUniqueRecvDispls[n];
    hlong end   = vertexUniqueRecvDispls[n+1];
    for(hlong v1=start;v1<end;++v1){ // vertex v1 to be sent back with list of conns
      for(hlong v2=start;v2<end;++v2){
	vertexHaloSendList[cnt] = vertexRecvList[v2];
	vertexHaloSendList[cnt].otherRank = vertexRecvList[v1].sourceRank;
	++cnt;
      }
    }
  }

  hlong NvertexHaloSend = cnt;
  
  // sort Halo send list based on source rank
  qsort(vertexRecvList, NvertexHaloSend, sizeof(vertex_t), compareSourceRank);   // check qsort counts

  // now figure out how many halo vertices to expect
  hlong *vertexHaloRecvCounts = (hlong*) calloc(mesh->size, sizeof(hlong));
  MPI_Alltoall(vertexHaloSendCounts, 1, MPI_HLONG,
	       vertexHaloRecvCounts, 1, MPI_HLONG,
	       mesh->comm);

  // find displacements for
  hlong *vertexHaloSendDispls = (hlong*) calloc(mesh->size+1, sizeof(hlong));
  hlong *vertexHaloRecvDispls = (hlong*) calloc(mesh->size+1, sizeof(hlong));
  for(int r=0;r<mesh->size;++r){
    vertexHaloSendCounts[r] *= sizeof(vertex_t);
    vertexHaloRecvCounts[r] *= sizeof(vertex_t);
    vertexHaloSendDispls[r+1] = vertexHaloSendDispls[r] + vertexHaloSendCounts[r];
    vertexHaloRecvDispls[r+1] = vertexHaloRecvDispls[r] + vertexHaloRecvCounts[r];
  }
  hlong NvertexHaloRecv = vertexHaloRecvDispls[mesh->size];
  
  vertex_t *vertexHaloRecvList = (vertex_t*) calloc(NvertexHaloRecv, sizeof(vertex_t)); // hack-hack-hack
  
  MPI_Alltoallv(vertexHaloSendList, vertexHaloSendCounts, vertexHaloSendDispls, MPI_CHAR,
		vertexHaloRecvList, vertexHaloRecvCounts, vertexHaloRecvDispls, MPI_CHAR,
		mesh->comm);

  // finally we now have a list of all elements that are needed to form the 1-ring (to rule them all)

  // sort the list by "other rank then element"
  qsort(vertexHaloRecvList, NvertexHaloRecv, sizeof(vertex_t), compareOtherRankElement);   // check qsort counts

  // remove duplicates and local elements from halo list in place
  hlong NnonLocalHaloElements = 0;
  for(hlong v=0;v<NvertexHaloRecv;++v){
    if(vertexHaloRecvList[v].sourceRank != mesh->rank){
      printf("%d %d %d ([rank] receives [element] from [other rank])\n",
	     mesh->rank, vertexHaloRecvList[v].element, vertexHaloRecvList[v].otherRank);
      ++NnonLocalHaloElements;
    }
  }

  MPI_Finalize();
  exit(0);
  
  vertex_t *nonLocalHaloElements = (vertex_t*) calloc(NnonLocalHaloElements, sizeof(vertex_t));
  cnt = 0;

  for(hlong v=0;v<NvertexHaloRecv;++v){
    if(vertexHaloRecvList[v].otherRank != mesh->rank){ // rule out local elements
      if(vertexHaloRecvList[v].otherRank == nonLocalHaloElements[cnt].otherRank){ // in same rank as last added
	if(vertexHaloRecvList[v].element != nonLocalHaloElements[cnt].element){ // but not the same element (not dupe)
	  nonLocalHaloElements[cnt] = vertexHaloRecvList[v];
	  ++cnt;
	}
      }else{
	nonLocalHaloElements[cnt] = vertexHaloRecvList[v]; 
	++cnt;
      }
    }
  }

  NnonLocalHaloElements = cnt;

#if 0
  for(hlong v=0;v<NnonLocalHaloElements;++v){
    printf("%d %d %d ([rank] receives [element] from [other rank])\n",
	   mesh->rank, nonLocalHaloElements[v].element, nonLocalHaloElements[v].otherRank);
  }
#endif
}
