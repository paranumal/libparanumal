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

typedef struct {

  dlong element, elementN;
  int face, faceN, rankN;

}facePair_t;

/* comparison function that orders halo element/face
   based on their indexes */
static int compareHaloFaces(const void *a,
                     const void *b){

  facePair_t *fa = (facePair_t*) a;
  facePair_t *fb = (facePair_t*) b;

  if(fa->rankN < fb->rankN) return -1;
  if(fa->rankN > fb->rankN) return +1;

  if(fa->elementN < fb->elementN) return -1;
  if(fa->elementN > fb->elementN) return +1;

  if(fa->faceN < fb->faceN) return -1;
  if(fa->faceN > fb->faceN) return +1;

  return 0;
}


// set up halo infomation for inter-processor MPI
// exchange of trace nodes
void mesh_t::HaloSetup(){

  // non-blocking MPI isend/irecv requests (used in meshHaloExchange)
  haloSendRequests = calloc(size, sizeof(MPI_Request));
  haloRecvRequests = calloc(size, sizeof(MPI_Request));

  // count number of halo element nodes to swap
  totalHaloPairs = 0;
  NhaloPairs = (int*) calloc(size, sizeof(int));
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      int r = EToP[e*Nfaces+f]; // rank of neighbor
      if(r!=-1){
        totalHaloPairs += 1;
        NhaloPairs[r] += 1;
      }
    }
  }

  // count number of MPI messages in halo exchange
  NhaloMessages = 0;
  for(int r=0;r<size;++r)
    if(NhaloPairs[r])
      ++NhaloMessages;

  // create a list of element/faces with halo neighbor
  facePair_t *haloElements =
    (facePair_t*) calloc(totalHaloPairs, sizeof(facePair_t));

  dlong cnt = 0;
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      dlong ef = e*Nfaces+f;
      if(EToP[ef]!=-1){
        haloElements[cnt].element  = e;
        haloElements[cnt].face     = f;
        haloElements[cnt].elementN = EToE[ef];
        haloElements[cnt].faceN    = EToF[ef];
        haloElements[cnt].rankN    = EToP[ef];
        ++cnt;
      }
    }
  }

  // sort the face pairs in order the destination requires
  qsort(haloElements, totalHaloPairs, sizeof(facePair_t), compareHaloFaces);

  // record the outgoing order for elements
  haloElementList = (dlong*) calloc(totalHaloPairs, sizeof(dlong));
  for(dlong i=0;i<totalHaloPairs;++i){
    dlong e = haloElements[i].element;
    haloElementList[i] = e;
  }

  // record the outgoing node ids for trace nodes
  haloGetNodeIds = (dlong*) calloc(totalHaloPairs*Nfp, sizeof(dlong));
  haloPutNodeIds = (dlong*) calloc(totalHaloPairs*Nfp, sizeof(dlong));

  cnt = 0;
  for(dlong i=0;i<totalHaloPairs;++i){
    dlong eM = haloElements[i].element;
    int fM = haloElements[i].face;
    int fP = haloElements[i].faceN;
    for(int n=0;n<Nfp;++n){
      haloGetNodeIds[cnt] = eM*Np + faceNodes[fM*Nfp+n];
      ++cnt;
    }
  }

  // now arrange for incoming nodes
  cnt = Nelements;
  dlong ncnt = 0;
  for(int r=0;r<size;++r){
    for(dlong e=0;e<Nelements;++e){
      for(int f=0;f<Nfaces;++f){
	dlong ef = e*Nfaces+f;
	if(EToP[ef]==r){
	  EToE[ef] = cnt;
	  int fP = EToF[ef];
	  for(int n=0;n<Nfp;++n){
	    haloPutNodeIds[ncnt] = cnt*Np + faceNodes[fP*Nfp+n];
	    ++ncnt;
	  }
	  ++cnt; // next halo element
	}
      }
    }
  }

  // create halo extension for x,y arrays
  dlong totalHaloNodes = totalHaloPairs*Np;
  dlong localNodes     = Nelements*Np;

  // temporary send buffer
  dfloat *sendBuffer = (dfloat*) calloc(totalHaloNodes, sizeof(dfloat));

  // extend x,y arrays to hold coordinates of node coordinates of elements in halo
  x = (dfloat*) realloc(x, (localNodes+totalHaloNodes)*sizeof(dfloat));
  y = (dfloat*) realloc(y, (localNodes+totalHaloNodes)*sizeof(dfloat));
  if(dim==3)
    z = (dfloat*) realloc(z, (localNodes+totalHaloNodes)*sizeof(dfloat));

  // send halo data and recv into extended part of arrays
  this->HaloExchange(Np*sizeof(dfloat), x, sendBuffer, x + localNodes);
  this->HaloExchange(Np*sizeof(dfloat), y, sendBuffer, y + localNodes);
  if(dim==3)
    this->HaloExchange(Np*sizeof(dfloat), z, sendBuffer, z + localNodes);

  // grab EX,EY,EZ from halo
  EX = (dfloat*) realloc(EX, (Nelements+totalHaloPairs)*Nverts*sizeof(dfloat));
  EY = (dfloat*) realloc(EY, (Nelements+totalHaloPairs)*Nverts*sizeof(dfloat));
  if(dim==3)
    EZ = (dfloat*) realloc(EZ, (Nelements+totalHaloPairs)*Nverts*sizeof(dfloat));

  // send halo data and recv into extended part of arrays
  this->HaloExchange(Nverts*sizeof(dfloat), EX, sendBuffer, EX + Nverts*Nelements);
  this->HaloExchange(Nverts*sizeof(dfloat), EY, sendBuffer, EY + Nverts*Nelements);
  if(dim==3)
    this->HaloExchange(Nverts*sizeof(dfloat), EZ, sendBuffer, EZ + Nverts*Nelements);

  free(haloElements);
  free(sendBuffer);
}
