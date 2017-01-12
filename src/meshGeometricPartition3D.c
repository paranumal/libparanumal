#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mesh3D.h"

// 20 bits per coordinate
#define bitRange 20

// spread bits of i by introducing zeros between binary bits
unsigned long long int bitSplitter(unsigned int i){
  
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
unsigned long long int mortonIndex3D(unsigned int ix, unsigned int iy, unsigned int iz){
  
  // spread bits of ix apart (introduce zeros)
  unsigned long long int sx = bitSplitter(ix);
  unsigned long long int sy = bitSplitter(iy);
  unsigned long long int sz = bitSplitter(iz);
  
  // interleave bits of ix and iy
  unsigned long long int mi = sx | (sy<<1) | (sz<<2); 
  
  return mi;
}

// capsule for element vertices + Morton index
typedef struct {
  
  unsigned long long int index;
  
  iint element;

  // use 8 for maximum vertices per element
  iint v[8];

  dfloat EX[8], EY[8], EZ[8];

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

// geometric partition of elements in 3D mesh using Morton ordering + parallelSort
void meshGeometricPartition3D(mesh3D *mesh){

  iint maxNelements;
  MPI_Allreduce(&(mesh->Nelements), &maxNelements, 1, MPI_IINT, MPI_MAX,
		MPI_COMM_WORLD);
  maxNelements = 2*((maxNelements+1)/2);
  
  // fix maxNelements
  element_t *elements 
    = (element_t*) calloc(maxNelements, sizeof(element_t));

  // local bounding box of element centers
  dfloat minvx = 1e9, maxvx = -1e9;
  dfloat minvy = 1e9, maxvy = -1e9;
  dfloat minvz = 1e9, maxvz = -1e9;

  // compute element centers on this process
  for(iint n=0;n<mesh->Nverts*mesh->Nelements;++n){
    minvx = mymin(minvx, mesh->EX[n]);
    maxvx = mymax(minvx, mesh->EX[n]);
    minvy = mymin(minvy, mesh->EY[n]);
    maxvy = mymax(minvy, mesh->EY[n]);
    minvz = mymin(minvz, mesh->EZ[n]);
    maxvz = mymax(minvz, mesh->EZ[n]);
  }
  
  // find global bounding box of element centers
  dfloat gminvx, gminvy, gminvz, gmaxvx, gmaxvy, gmaxvz;
  MPI_Allreduce(&minvx, &gminvx, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&minvy, &gminvy, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&minvz, &gminvz, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&maxvx, &gmaxvx, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&maxvy, &gmaxvy, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&maxvz, &gmaxvz, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);

  // choose sub-range of Morton lattice coordinates to embed element centers in
  unsigned long long int Nboxes = (((unsigned long long int)1)<<(bitRange-1));
  
  // compute Morton index for each element
  for(iint e=0;e<mesh->Nelements;++e){

    // element center coordinates
    dfloat cx = 0, cy = 0, cz = 0;
    for(iint n=0;n<mesh->Nverts;++n){
      cx += mesh->EX[e*mesh->Nverts+n];
      cy += mesh->EY[e*mesh->Nverts+n];
      cz += mesh->EZ[e*mesh->Nverts+n];
    }
    cx /= mesh->Nverts;
    cy /= mesh->Nverts;
    cz /= mesh->Nverts;

    // encapsulate element, vertices, Morton index, vertex coordinates
    elements[e].element = e;
    for(iint n=0;n<mesh->Nverts;++n){
      elements[e].v[n] = mesh->EToV[e*mesh->Nverts+n];
      elements[e].EX[n] = mesh->EX[e*mesh->Nverts+n];
      elements[e].EY[n] = mesh->EY[e*mesh->Nverts+n];
      elements[e].EZ[n] = mesh->EZ[e*mesh->Nverts+n];
    }

    unsigned long long int ix = (cx-gminvx)*Nboxes/(gmaxvx-gminvx);
    unsigned long long int iy = (cy-gminvy)*Nboxes/(gmaxvy-gminvy);
    unsigned long long int iz = (cz-gminvz)*Nboxes/(gmaxvz-gminvz);
			
    elements[e].index = mortonIndex3D(ix, iy, iz);
  }

  // pad element array with dummy elements
  for(iint e=mesh->Nelements;e<maxNelements;++e){
    elements[e].element = -1;
    elements[e].index = mortonIndex3D(Nboxes+1, Nboxes+1, Nboxes+1);
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
  free(mesh->EZ);

  mesh->Nelements = cnt;
  mesh->EToV = (iint*) calloc(cnt*mesh->Nverts, sizeof(iint));
  mesh->EX = (dfloat*) calloc(cnt*mesh->Nverts, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(cnt*mesh->Nverts, sizeof(dfloat));
  mesh->EZ = (dfloat*) calloc(cnt*mesh->Nverts, sizeof(dfloat));

  cnt = 0;
  for(iint e=0;e<maxNelements;++e){
    if(elements[e].element != -1){
      for(iint n=0;n<mesh->Nverts;++n){
	mesh->EToV[cnt*mesh->Nverts + n] = elements[e].v[n];
	mesh->EX[cnt*mesh->Nverts + n]   = elements[e].EX[n];
	mesh->EY[cnt*mesh->Nverts + n]   = elements[e].EY[n];
	mesh->EZ[cnt*mesh->Nverts + n]   = elements[e].EZ[n];
      }
      ++cnt;
    }
  }
}
