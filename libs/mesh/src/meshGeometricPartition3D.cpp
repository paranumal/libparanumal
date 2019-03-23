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

// 20 bits per coordinate
#define bitRange 20

// spread bits of i by introducing zeros between binary bits
static unsigned long long int bitSplitter3D(unsigned int i){

  unsigned long long int mask = 1;
  unsigned long long int li = i;
  unsigned long long int lj = 0;

  for(int b=0;b<bitRange;++b){
    lj |= ((li & mask) << 2*b); // bit b moves to bit 3b
    mask <<= 1;
  }

  return lj;
}

// compute Morton index of (ix,iy) relative to a bitRange x bitRange  Morton lattice
static unsigned long long int mortonIndex3D(unsigned int ix, unsigned int iy, unsigned int iz){

  // spread bits of ix apart (introduce zeros)
  unsigned long long int sx = bitSplitter3D(ix);
  unsigned long long int sy = bitSplitter3D(iy);
  unsigned long long int sz = bitSplitter3D(iz);

  // interleave bits of ix and iy
  unsigned long long int mi = sx | (sy<<1) | (sz<<2);

  return mi;
}

// capsule for element vertices + Morton index
typedef struct {

  unsigned long long int index;

  dlong element;

  int type;

  // use 8 for maximum vertices per element
  hlong v[8];

  dfloat EX[8], EY[8], EZ[8];

}element_t;

// compare the Morton indices for two element capsules
static int compareElements(const void *a, const void *b){

  element_t *ea = (element_t*) a;
  element_t *eb = (element_t*) b;

  if(ea->index < eb->index) return -1;
  if(ea->index > eb->index) return  1;

  return 0;

}

// stub for the match function needed by parallelSort
static void bogusMatch3D(void *a, void *b){ }

// geometric partition of elements in 3D mesh using Morton ordering + parallelSort
void mesh3D::GeometricPartition(){

  dlong maxNelements;
  MPI_Allreduce(&(Nelements), &maxNelements, 1, MPI_DLONG, MPI_MAX,
		comm);
  maxNelements = 2*((maxNelements+1)/2);

  // fix maxNelements
  element_t *elements
    = (element_t*) calloc(maxNelements, sizeof(element_t));

  // local bounding box of element centers
  dfloat minvx = 1e9, maxvx = -1e9;
  dfloat minvy = 1e9, maxvy = -1e9;
  dfloat minvz = 1e9, maxvz = -1e9;

  // compute element centers on this process
  for(dlong n=0;n<Nverts*Nelements;++n){
    minvx = mymin(minvx, EX[n]);
    maxvx = mymax(maxvx, EX[n]);
    minvy = mymin(minvy, EY[n]);
    maxvy = mymax(maxvy, EY[n]);
    minvz = mymin(minvz, EZ[n]);
    maxvz = mymax(maxvz, EZ[n]);
  }

  // find global bounding box of element centers
  dfloat gminvx, gminvy, gminvz, gmaxvx, gmaxvy, gmaxvz;
  MPI_Allreduce(&minvx, &gminvx, 1, MPI_DFLOAT, MPI_MIN, comm);
  MPI_Allreduce(&minvy, &gminvy, 1, MPI_DFLOAT, MPI_MIN, comm);
  MPI_Allreduce(&minvz, &gminvz, 1, MPI_DFLOAT, MPI_MIN, comm);
  MPI_Allreduce(&maxvx, &gmaxvx, 1, MPI_DFLOAT, MPI_MAX, comm);
  MPI_Allreduce(&maxvy, &gmaxvy, 1, MPI_DFLOAT, MPI_MAX, comm);
  MPI_Allreduce(&maxvz, &gmaxvz, 1, MPI_DFLOAT, MPI_MAX, comm);

  // choose sub-range of Morton lattice coordinates to embed element centers in
  unsigned long long int Nboxes = (((unsigned long long int)1)<<(bitRange-1));

  // compute Morton index for each element
  for(dlong e=0;e<Nelements;++e){

    // element center coordinates
    dfloat cx = 0, cy = 0, cz = 0;
    for(int n=0;n<Nverts;++n){
      cx += EX[e*Nverts+n];
      cy += EY[e*Nverts+n];
      cz += EZ[e*Nverts+n];
    }
    cx /= Nverts;
    cy /= Nverts;
    cz /= Nverts;

    // encapsulate element, vertices, Morton index, vertex coordinates
    elements[e].element = e;
    for(int n=0;n<Nverts;++n){
      elements[e].v[n] = EToV[e*Nverts+n];
      elements[e].EX[n] = EX[e*Nverts+n];
      elements[e].EY[n] = EY[e*Nverts+n];
      elements[e].EZ[n] = EZ[e*Nverts+n];
    }

    elements[e].type = elementInfo[e];

    dfloat maxlength = mymax(gmaxvx-gminvx, mymax(gmaxvy-gminvy, gmaxvz-gminvz));

    // avoid stretching axes
    unsigned long long int ix = (cx-gminvx)*Nboxes/maxlength;
    unsigned long long int iy = (cy-gminvy)*Nboxes/maxlength;
    unsigned long long int iz = (cz-gminvz)*Nboxes/maxlength;

    elements[e].index = mortonIndex3D(ix, iy, iz);
  }

  // pad element array with dummy elements
  for(dlong e=Nelements;e<maxNelements;++e){
    elements[e].element = -1;
    elements[e].index = mortonIndex3D(Nboxes+1, Nboxes+1, Nboxes+1);
  }

  // odd-even parallel sort of element capsules based on their Morton index
  parallelSort(size, rank, comm,
	       maxNelements, elements, sizeof(element_t),
	       compareElements,
	       bogusMatch3D);

#if 0
  // count number of elements that end up on this process
  int cnt = 0;
  for(int e=0;e<maxNelements;++e)
    cnt += (elements[e].element != -1);

  // reset number of elements and element-to-vertex connectivity from returned capsules
  free(EToV);
  free(EX);
  free(EY);
  free(EZ);

  Nelements = cnt;
  EToV = (int*) calloc(cnt*Nverts, sizeof(int));
  EX = (dfloat*) calloc(cnt*Nverts, sizeof(dfloat));
  EY = (dfloat*) calloc(cnt*Nverts, sizeof(dfloat));
  EZ = (dfloat*) calloc(cnt*Nverts, sizeof(dfloat));

  cnt = 0;
  for(int e=0;e<maxNelements;++e){
    if(elements[e].element != -1){
      for(int n=0;n<Nverts;++n){
	EToV[cnt*Nverts + n] = elements[e].v[n];
	EX[cnt*Nverts + n]   = elements[e].EX[n];
	EY[cnt*Nverts + n]   = elements[e].EY[n];
	EZ[cnt*Nverts + n]   = elements[e].EZ[n];
      }
      ++cnt;
    }
  }
#else
  // compress and renumber elements
  dlong sk  = 0;
  for(dlong e=0;e<maxNelements;++e){
    if(elements[e].element != -1){
      elements[sk] = elements[e];
      ++sk;
    }
  }

  dlong localNelements = sk;

  /// redistribute elements to improve balancing
  dlong *globalNelements = (dlong *) calloc(size,sizeof(dlong));
  hlong *starts = (hlong *) calloc(size+1,sizeof(hlong));

  MPI_Allgather(&localNelements, 1, MPI_DLONG, globalNelements, 1,  MPI_DLONG, comm);

  for(int r=0;r<size;++r)
    starts[r+1] = starts[r]+globalNelements[r];

  hlong allNelements = starts[size];

  // decide how many to keep on each process
  hlong chunk = allNelements/size;
  int remainder = (int) (allNelements - chunk*size);

  int *Nsend = (int *) calloc(size, sizeof(int));
  int *Nrecv = (int *) calloc(size, sizeof(int));
  // int *Ncount = (int *) calloc(size, sizeof(int));
  int *sendOffsets = (int*) calloc(size, sizeof(int));
  int *recvOffsets = (int*) calloc(size, sizeof(int));


  // Make the MPI_ELEMENT_T data type
  MPI_Datatype MPI_ELEMENT_T;
  MPI_Datatype dtype[7] = {MPI_LONG_LONG_INT, MPI_DLONG, MPI_INT,
                            MPI_HLONG, MPI_DFLOAT, MPI_DFLOAT, MPI_DFLOAT};
  int blength[7] = {1, 1, 1, 8, 8, 8, 8};
  MPI_Aint addr[7], displ[7];
  MPI_Get_address ( &(elements[0]        ), addr+0);
  MPI_Get_address ( &(elements[0].element), addr+1);
  MPI_Get_address ( &(elements[0].type   ), addr+2);
  MPI_Get_address ( &(elements[0].v[0]   ), addr+3);
  MPI_Get_address ( &(elements[0].EX[0]  ), addr+4);
  MPI_Get_address ( &(elements[0].EY[0]  ), addr+5);
  MPI_Get_address ( &(elements[0].EZ[0]  ), addr+6);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  displ[3] = addr[3] - addr[0];
  displ[4] = addr[4] - addr[0];
  displ[5] = addr[5] - addr[0];
  displ[6] = addr[6] - addr[0];
  MPI_Type_create_struct (7, blength, displ, dtype, &MPI_ELEMENT_T);
  MPI_Type_commit (&MPI_ELEMENT_T);


  for(dlong e=0;e<localNelements;++e){

    // global element index
    elements[e].element = starts[rank]+e;

    // 0, chunk+1, 2*(chunk+1) ..., remainder*(chunk+1), remainder*(chunk+1) + chunk
    int r;
    if(elements[e].element<remainder*(chunk+1))
      r = elements[e].element/(chunk+1);
    else
      r = remainder + ((elements[e].element-remainder*(chunk+1))/chunk);

    ++Nsend[r];
  }

  // find send offsets
  for(int r=1;r<size;++r)
    sendOffsets[r] = sendOffsets[r-1] + Nsend[r-1];

  // exchange byte counts
  MPI_Alltoall(Nsend, 1, MPI_INT, Nrecv, 1, MPI_INT, comm);

  // count incoming clusters
  dlong newNelements = 0;
  for(int r=0;r<size;++r)
    newNelements += Nrecv[r];

  for(int r=1;r<size;++r)
    recvOffsets[r] = recvOffsets[r-1] + Nrecv[r-1];

  element_t *tmpElements = (element_t *) calloc(newNelements, sizeof(element_t));

  // exchange parallel clusters
  MPI_Alltoallv(elements, Nsend, sendOffsets, MPI_ELEMENT_T,
                tmpElements, Nrecv, recvOffsets, MPI_ELEMENT_T, comm);

  MPI_Barrier(comm);
  MPI_Type_free(&MPI_ELEMENT_T);

  // replace elements with inbound elements
  if (elements) free(elements);
  elements = tmpElements;

  // reset number of elements and element-to-vertex connectivity from returned capsules
  free(EToV);
  free(EX);
  free(EY);
  free(EZ);
  free(elementInfo);

  Nelements = newNelements;
  EToV = (hlong*) calloc(newNelements*Nverts, sizeof(hlong));
  EX = (dfloat*) calloc(newNelements*Nverts, sizeof(dfloat));
  EY = (dfloat*) calloc(newNelements*Nverts, sizeof(dfloat));
  EZ = (dfloat*) calloc(newNelements*Nverts, sizeof(dfloat));
  elementInfo = (hlong*) calloc(newNelements, sizeof(hlong));

  for(dlong e=0;e<newNelements;++e){
    for(int n=0;n<Nverts;++n){
      EToV[e*Nverts + n] = elements[e].v[n];
      EX[e*Nverts + n]   = elements[e].EX[n];
      EY[e*Nverts + n]   = elements[e].EY[n];
      EZ[e*Nverts + n]   = elements[e].EZ[n];
    }
    elementInfo[e] = elements[e].type;
  }
  if (elements) free(elements);
#endif
}
