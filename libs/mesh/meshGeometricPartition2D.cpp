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
#include "mesh2D.hpp"

#define bitRange 15

#if 0

/// THIS SECTION ------------------------------------------------------------------------------------>
// taken from: http://and-what-happened.blogspot.com/2011/08/fast-2d-and-3d-hilbert-curves-and.html

unsigned int Morton_2D_Encode_16bit( unsigned int index1, unsigned int index2 )
{ // pack 2 16-bit indices into a 32-bit Morton code
  index1 &= 0x0000ffff;
  index2 &= 0x0000ffff;
  index1 |= ( index1 << 8 );
  index2 |= ( index2 << 8 );
  index1 &= 0x00ff00ff;
  index2 &= 0x00ff00ff;
  index1 |= ( index1 << 4 );
  index2 |= ( index2 << 4 );
  index1 &= 0x0f0f0f0f;
  index2 &= 0x0f0f0f0f;
  index1 |= ( index1 << 2 );
  index2 |= ( index2 << 2 );
  index1 &= 0x33333333;
  index2 &= 0x33333333;
  index1 |= ( index1 << 1 );
  index2 |= ( index2 << 1 );
  index1 &= 0x55555555;
  index2 &= 0x55555555;
  return( index1 | ( index2 << 1 ) );
}

unsigned int MortonToHilbert2D( const unsigned int morton, const unsigned int bits )
{
  unsigned int hilbert = 0;
  unsigned int remap = 0xb4;
  unsigned int block = ( bits << 1 );
  while( block )
    {
      block -= 2;
      unsigned int mcode = ( ( morton >> block ) & 3 );
      unsigned int hcode = ( ( remap >> ( mcode << 1 ) ) & 3 );
      remap ^= ( 0x82000028 >> ( hcode << 3 ) );
      hilbert = ( ( hilbert << 2 ) + hcode );
    }
  return( hilbert );
}


unsigned int hilbert2D(unsigned int index1, unsigned int index2){

  unsigned int morton = Morton_2D_Encode_16bit(index1,index2);

  return MortonToHilbert2D(morton, 16);
}

/// THIS SECTION TO HERE <--------------------------------------------------------------------------------

// spread bits of i by introducing zeros between binary bits
unsigned long long int bitSplitter(unsigned int i){

  unsigned long long int mask = 1;
  unsigned long long int li = i;
  unsigned long long int lj = 0;

  for(int b=0;b<bitRange;++b){
    lj |=  (li & mask) << b;
    mask <<= 1;
  }

  return lj;

}

// compute Morton index of (ix,iy) relative to a bitRange x bitRange  Morton lattice
unsigned long long int mortonIndex2D(unsigned int ix, unsigned int iy){

  // spread bits of ix apart (introduce zeros)
  unsigned long long int sx = bitSplitter(ix);
  unsigned long long int sy = bitSplitter(iy);

  // interleave bits of ix and iy
  unsigned long long int mi = sx | (sy<<1);

  return mi;
}

#else

// from: https://en.wikipedia.org/wiki/Hilbert_curve

//rotate/flip a quadrant appropriately
static void rot(unsigned int n, unsigned int *x, unsigned int *y, unsigned int rx, unsigned int ry) {
  if (ry == 0) {
    if (rx == 1) {
      *x = n-1 - *x;
      *y = n-1 - *y;
    }

    //Swap x and y
    int t  = *x;
    *x = *y;
    *y = t;
  }
}



//convert (x,y) to d
static unsigned int hilbert2D (unsigned int n, unsigned int x, unsigned int y) {
  unsigned int rx, ry, s, d=0;
  for (s=n/2; s>0; s/=2) {
    rx = (x & s) > 0;
    ry = (y & s) > 0;
    d += s * s * ((3 * rx) ^ ry);
    rot(s, &x, &y, rx, ry);
  }
  return d;
}

#endif

// capsule for element vertices + Morton index
typedef struct {

  unsigned long long int index;

  dlong element;

  int type;

  // 4 for maximum number of vertices per element in 2D
  hlong v[4];

  dfloat EX[4], EY[4];

}element_t;

// compare the Morton indices for two element capsules
static int compareElements2D(const void *a, const void *b){

  element_t *ea = (element_t*) a;
  element_t *eb = (element_t*) b;

  if(ea->index < eb->index) return -1;
  if(ea->index > eb->index) return  1;

  return 0;

}

// stub for the match function needed by parallelSort
static void bogusMatch(void *a, void *b){ }

// geometric partition of elements in 2D mesh using Morton ordering + parallelSort
void mesh2D::GeometricPartition(){

  dlong maxNelements;
  MPI_Allreduce(&(Nelements), &maxNelements, 1, MPI_DLONG, MPI_MAX, comm);
  maxNelements = 2*((maxNelements+1)/2);

  // fix maxNelements
  element_t *elements
    = (element_t*) calloc(maxNelements, sizeof(element_t));

  // local bounding box of element centers
  dfloat mincx = 1e9, maxcx = -1e9;
  dfloat mincy = 1e9, maxcy = -1e9;

  // compute element centers on this process
  for(dlong e=0;e<Nelements;++e){
    dfloat cx = 0, cy = 0;
    for(int n=0;n<Nverts;++n){
      cx += EX[e*Nverts+n];
      cy += EY[e*Nverts+n];
    }
    cx /= Nverts;
    cy /= Nverts;

    mincx = mymin(mincx, cx);
    maxcx = mymax(maxcx, cx);
    mincy = mymin(mincy, cy);
    maxcy = mymax(maxcy, cy);
  }

  dfloat delta = 1e-1;
  mincx -= delta;
  mincy -= delta;
  maxcx += delta;
  maxcy += delta;

  // find global bounding box of element centers
  dfloat gmincx, gmincy, gmaxcx, gmaxcy;
  MPI_Allreduce(&mincx, &gmincx, 1, MPI_DFLOAT, MPI_MIN, comm);
  MPI_Allreduce(&mincy, &gmincy, 1, MPI_DFLOAT, MPI_MIN, comm);
  MPI_Allreduce(&maxcx, &gmaxcx, 1, MPI_DFLOAT, MPI_MAX, comm);
  MPI_Allreduce(&maxcy, &gmaxcy, 1, MPI_DFLOAT, MPI_MAX, comm);

  dfloat maxlength = mymax(gmaxcx-gmincx, gmaxcy-gmincy);

  // choose sub-range of Morton lattice coordinates to embed element centers in
  unsigned int Nboxes = (((unsigned int)1)<<(bitRange));

  // compute Morton index for each element
  for(dlong e=0;e<Nelements;++e){

    // element center coordinates
    dfloat cx = 0, cy = 0;
    for(int n=0;n<Nverts;++n){
      cx += EX[e*Nverts+n];
      cy += EY[e*Nverts+n];
    }
    cx /= Nverts;
    cy /= Nverts;

    // encapsulate element, vertices, Morton index, vertex coordinates
    elements[e].element = e;
    for(int n=0;n<Nverts;++n){
      elements[e].v[n] = EToV[e*Nverts+n];
      elements[e].EX[n] = EX[e*Nverts+n];
      elements[e].EY[n] = EY[e*Nverts+n];
    }

    elements[e].type = elementInfo[e];

    unsigned int ix = (cx-gmincx)*Nboxes/maxlength;
    unsigned int iy = (cy-gmincy)*Nboxes/maxlength;

    //elements[e].index = mortonIndex2D(ix, iy);
    elements[e].index = hilbert2D(Nboxes, ix, iy);
  }

  // pad element array with dummy elements
  for(dlong e=Nelements;e<maxNelements;++e){
    elements[e].element = -1;

    elements[e].index = hilbert2D(Nboxes, Nboxes-1, Nboxes-1);

    //    elements[e].index = hilbert2D(Nboxes+1, Nboxes+1);
    //    elements[e].index = mortonIndex2D(Nboxes+1, Nboxes+1);
  }

  // odd-even parallel sort of element capsules based on their Morton index
  parallelSort(size, rank, comm,
	       maxNelements, elements, sizeof(element_t),
	       compareElements2D,
	       bogusMatch);


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
  // TODO: We need a safer version of this for very large meshes.
  // if dlong is a long long int Nsend and/or sendOffsets may overflow int
  dlong *globalNelements = (dlong *) calloc(size,sizeof(dlong));
  hlong *starts = (hlong *) calloc(size+1,sizeof(hlong));

  MPI_Allgather(&localNelements, 1, MPI_DLONG, globalNelements, 1,  MPI_DLONG, comm);

  for(int rr=0;rr<size;++rr)
    starts[rr+1] = starts[rr]+globalNelements[rr];

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
  MPI_Datatype dtype[6] = {MPI_LONG_LONG_INT, MPI_DLONG, MPI_INT,
                            MPI_HLONG, MPI_DFLOAT, MPI_DFLOAT};
  int blength[6] = {1, 1, 1, 4, 4, 4};
  MPI_Aint addr[6], displ[6];
  MPI_Get_address ( &(elements[0]        ), addr+0);
  MPI_Get_address ( &(elements[0].element), addr+1);
  MPI_Get_address ( &(elements[0].type   ), addr+2);
  MPI_Get_address ( &(elements[0].v[0]   ), addr+3);
  MPI_Get_address ( &(elements[0].EX[0]  ), addr+4);
  MPI_Get_address ( &(elements[0].EY[0]  ), addr+5);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  displ[3] = addr[3] - addr[0];
  displ[4] = addr[4] - addr[0];
  displ[5] = addr[5] - addr[0];
  MPI_Type_create_struct (6, blength, displ, dtype, &MPI_ELEMENT_T);
  MPI_Type_commit (&MPI_ELEMENT_T);

  for(dlong e=0;e<localNelements;++e){

    // global element index
    elements[e].element = starts[rank]+e;

    // 0, chunk+1, 2*(chunk+1) ..., remainder*(chunk+1), remainder*(chunk+1) + chunk
    int rr;
    if(elements[e].element<remainder*(chunk+1))
      rr = elements[e].element/(chunk+1);
    else
      rr = remainder + ((elements[e].element-remainder*(chunk+1))/chunk);

    ++Nsend[rr];
  }

  // find send offsets
  for(int rr=1;rr<size;++rr)
    sendOffsets[rr] = sendOffsets[rr-1] + Nsend[rr-1];

  // exchange byte counts
  MPI_Alltoall(Nsend, 1, MPI_INT, Nrecv, 1, MPI_INT, comm);

  // count incoming clusters
  dlong newNelements = 0;
  for(int rr=0;rr<size;++rr)
    newNelements += Nrecv[rr];

  for(int rr=1;rr<size;++rr)
    recvOffsets[rr] = recvOffsets[rr-1] + Nrecv[rr-1];

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
  free(elementInfo);

  Nelements = newNelements;
  EToV = (hlong*) calloc(newNelements*Nverts, sizeof(hlong));
  EX = (dfloat*) calloc(newNelements*Nverts, sizeof(dfloat));
  EY = (dfloat*) calloc(newNelements*Nverts, sizeof(dfloat));
  elementInfo = (hlong*) calloc(newNelements, sizeof(hlong));

  for(dlong e=0;e<newNelements;++e){
    for(int n=0;n<Nverts;++n){
      EToV[e*Nverts + n] = elements[e].v[n];
      EX[e*Nverts + n]   = elements[e].EX[n];
      EY[e*Nverts + n]   = elements[e].EY[n];
    }
    elementInfo[e] = elements[e].type;
  }
  if (elements) free(elements);
}
