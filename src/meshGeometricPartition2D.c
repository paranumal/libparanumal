#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"

#define bitRange 30

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

  iint maxNelements;
  MPI_Allreduce(&(mesh->Nelements), &maxNelements, 1, MPI_IINT, MPI_MAX,
		MPI_COMM_WORLD);
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

  dfloat delta = 1e-4;
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
  unsigned long long int Nboxes = (((unsigned long long int)1)<<(bitRange-1));
  
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

    unsigned long long int ix = (cx-gmincx)*Nboxes/(gmaxcx-gmincx);
    unsigned long long int iy = (cy-gmincy)*Nboxes/(gmaxcy-gmincy);
			
    elements[e].index = mortonIndex2D(ix, iy);
  }

  // pad element array with dummy elements
  for(iint e=mesh->Nelements;e<maxNelements;++e){
    elements[e].element = -1;
    elements[e].index = mortonIndex2D(Nboxes+1, Nboxes+1);
  }

  // odd-even parallel sort of element capsules based on their Morton index
  parallelSort(maxNelements, elements, sizeof(element_t),
	       compareElements, 
	       bogusMatch);

  // count number of elements that end up on this process
  iint cnt = 0;
  for(iint e=0;e<maxNelements;++e)
    cnt += (elements[e].element != -1);

  // reset number of elements and element-to-vertex connectivity from returned capsules
  free(mesh->EToV);
  free(mesh->EX);
  free(mesh->EY);

  mesh->Nelements = cnt;
  mesh->EToV = (iint*) calloc(cnt*mesh->Nverts, sizeof(iint));
  mesh->EX = (dfloat*) calloc(cnt*mesh->Nverts, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(cnt*mesh->Nverts, sizeof(dfloat));

  cnt = 0;
  for(iint e=0;e<maxNelements;++e){
    if(elements[e].element != -1){
      for(iint n=0;n<mesh->Nverts;++n){
	mesh->EToV[cnt*mesh->Nverts + n] = elements[e].v[n];
	mesh->EX[cnt*mesh->Nverts + n]   = elements[e].EX[n];
	mesh->EY[cnt*mesh->Nverts + n]   = elements[e].EY[n];
      }
      ++cnt;
    }
  }

}
