#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"

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

#define bitRange 10

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

// capsule for element vertices + Morton index
typedef struct {

  unsigned long long int index;

  iint element;

  int type;

  // 4 for maximum number of vertices per element in 2D
  iint v[4];

  dfloat EX[4], EY[4];

}element_t;

// compare the Morton indices for two element capsules
int compareElements(const void *a, const void *b){

  element_t *ea = (element_t*) a;
  element_t *eb = (element_t*) b;

  if(ea->index < eb->index) return -1;
  if(ea->index > eb->index) return  1;

  return 0;

}

// stub for the match function needed by parallelSort
void bogusMatch(void *a, void *b){ }

// geometric partition of elements in 2D mesh using Morton ordering + parallelSort
void meshGeometricPartition2D(mesh2D *mesh){

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  iint maxNelements;
  MPI_Allreduce(&(mesh->Nelements), &maxNelements, 1, MPI_IINT, MPI_MAX, MPI_COMM_WORLD);
  maxNelements = 2*((maxNelements+1)/2);

  // fix maxNelements
  element_t *elements
    = (element_t*) calloc(maxNelements, sizeof(element_t));

  // local bounding box of element centers
  dfloat mincx = 1e9, maxcx = -1e9;
  dfloat mincy = 1e9, maxcy = -1e9;

  // compute element centers on this process
  for(iint e=0;e<mesh->Nelements;++e){
    dfloat cx = 0, cy = 0;
    for(iint n=0;n<mesh->Nverts;++n){
      cx += mesh->EX[e*mesh->Nverts+n];
      cy += mesh->EY[e*mesh->Nverts+n];
    }
    cx /= mesh->Nverts;
    cy /= mesh->Nverts;

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
  MPI_Allreduce(&mincx, &gmincx, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&mincy, &gmincy, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&maxcx, &gmaxcx, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&maxcy, &gmaxcy, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);

  // choose sub-range of Morton lattice coordinates to embed element centers in
  unsigned int Nboxes = (((unsigned int)1)<<(bitRange-1));

  // compute Morton index for each element
  for(iint e=0;e<mesh->Nelements;++e){

    // element center coordinates
    dfloat cx = 0, cy = 0;
    for(iint n=0;n<mesh->Nverts;++n){
      cx += mesh->EX[e*mesh->Nverts+n];
      cy += mesh->EY[e*mesh->Nverts+n];
    }
    cx /= mesh->Nverts;
    cy /= mesh->Nverts;

    // encapsulate element, vertices, Morton index, vertex coordinates
    elements[e].element = e;
    for(iint n=0;n<mesh->Nverts;++n){
      elements[e].v[n] = mesh->EToV[e*mesh->Nverts+n];
      elements[e].EX[n] = mesh->EX[e*mesh->Nverts+n];
      elements[e].EY[n] = mesh->EY[e*mesh->Nverts+n];
    }

    elements[e].type = mesh->elementInfo[e];

    unsigned int ix = (cx-gmincx)*Nboxes/(gmaxcx-gmincx);
    unsigned int iy = (cy-gmincy)*Nboxes/(gmaxcy-gmincy);

    //elements[e].index = mortonIndex2D(ix, iy);
    elements[e].index = hilbert2D(ix, iy);
  }

  // pad element array with dummy elements
  for(iint e=mesh->Nelements;e<maxNelements;++e){
    elements[e].element = -1;
    elements[e].index = hilbert2D(Nboxes+1, Nboxes+1);
    //    elements[e].index = mortonIndex2D(Nboxes+1, Nboxes+1);
  }

  // odd-even parallel sort of element capsules based on their Morton index
  parallelSort(maxNelements, elements, sizeof(element_t),
	       compareElements,
	       bogusMatch);


  // compress and renumber elements
  iint sk  = 0;
  for(iint e=0;e<maxNelements;++e){
    if(elements[e].element != -1){
      elements[sk] = elements[e];
      ++sk;
    }
  }

  iint localNelements = sk;

  /// redistribute elements to improve balancing
  iint *globalNelements = (iint *) calloc(size,sizeof(iint));
  iint *starts = (iint *) calloc(size+1,sizeof(iint));

  MPI_Allgather(&localNelements, 1, MPI_IINT, globalNelements, 1,  MPI_IINT, MPI_COMM_WORLD);

  for(iint r=0;r<size;++r)
    starts[r+1] = starts[r]+globalNelements[r];

  iint allNelements = starts[size];

  // decide how many to keep on each process
  int chunk = allNelements/size;
  int remainder = allNelements - chunk*size;

  iint *Nsend = (iint *) calloc(size, sizeof(iint));
  iint *Nrecv = (iint *) calloc(size, sizeof(iint));
  iint *Ncount = (iint *) calloc(size, sizeof(iint));
  iint *sendOffsets = (iint*) calloc(size, sizeof(iint));
  iint *recvOffsets = (iint*) calloc(size, sizeof(iint));


  for(iint e=0;e<localNelements;++e){

    // global element index
    elements[e].element = starts[rank]+e;

    // 0, chunk+1, 2*(chunk+1) ..., remainder*(chunk+1), remainder*(chunk+1) + chunk
    iint r;
    if(elements[e].element<remainder*(chunk+1))
      r = elements[e].element/(chunk+1);
    else
      r = remainder + ((elements[e].element-remainder*(chunk+1))/chunk);

    ++Nsend[r];
  }

  // find send offsets
  for(iint r=1;r<size;++r)
    sendOffsets[r] = sendOffsets[r-1] + Nsend[r-1];

  // exchange byte counts
  MPI_Alltoall(Nsend, 1, MPI_IINT, Nrecv, 1, MPI_IINT, MPI_COMM_WORLD);

  // count incoming clusters
  iint newNelements = 0;
  for(iint r=0;r<size;++r){
    newNelements += Nrecv[r];
    Nrecv[r] *= sizeof(element_t);
    Nsend[r] *= sizeof(element_t);
    sendOffsets[r] *= sizeof(element_t);
  }
  for(iint r=1;r<size;++r)
    recvOffsets[r] = recvOffsets[r-1] + Nrecv[r-1];

  element_t *tmpElements = (element_t *) calloc(newNelements, sizeof(element_t));

  // exchange parallel clusters
  MPI_Alltoallv(elements, Nsend, sendOffsets, MPI_CHAR,
                tmpElements, Nrecv, recvOffsets, MPI_CHAR, MPI_COMM_WORLD);

  // replace elements with inbound elements
  if (elements) free(elements);
  elements = tmpElements;

  // reset number of elements and element-to-vertex connectivity from returned capsules
  free(mesh->EToV);
  free(mesh->EX);
  free(mesh->EY);
  free(mesh->elementInfo);

  mesh->Nelements = newNelements;
  mesh->EToV = (iint*) calloc(newNelements*mesh->Nverts, sizeof(iint));
  mesh->EX = (dfloat*) calloc(newNelements*mesh->Nverts, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(newNelements*mesh->Nverts, sizeof(dfloat));
  mesh->elementInfo = (int*) calloc(newNelements, sizeof(int));

  for(iint e=0;e<newNelements;++e){
    for(iint n=0;n<mesh->Nverts;++n){
      mesh->EToV[e*mesh->Nverts + n] = elements[e].v[n];
      mesh->EX[e*mesh->Nverts + n]   = elements[e].EX[n];
      mesh->EY[e*mesh->Nverts + n]   = elements[e].EY[n];
    }
    mesh->elementInfo[e] = elements[e].type;
  }
}
