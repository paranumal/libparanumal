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

#include "elliptic.hpp"
#include "mesh/meshDefines2D.h"
#include "mesh/meshDefines3D.h"

#define nonZero_t parAlmond::parCOO::nonZero_t

typedef struct{
  hlong gnum;
  int   lnum;
}globalNode_t;

int compareGlobalNodes2(const void *a, const void *b){
  
  globalNode_t *ea = (globalNode_t*) a;
  globalNode_t *eb = (globalNode_t*) b;
  
  if(ea->gnum < eb->gnum) return -1;
  if(ea->gnum > eb->gnum) return +1;
  
  return 0;
}


void elliptic_t::BuildOperatorMatrixContinuousDevice(occa::memory &o_A,
						     dlong &Annz) {

  // build scanner and sorter
  occa::properties kernelInfo = mesh.props;
  
  if(sizeof(hlong)==8)
    kernelInfo["defines/hlong"]= "long long int";
  if(sizeof(hlong)==4)
    kernelInfo["defines/hlong"]= "int";
  
  hlong BIG_NUM = ((hlong)1) << (8*sizeof(hlong)-2);
  kernelInfo["defines/" "BIG_NUM"] = (const int64_t)BIG_NUM;
  
  if(sizeof(hlong)==8)
    kernelInfo["defines/hlong"]= "long long int";
  if(sizeof(hlong)==4)
    kernelInfo["defines/hlong"]= "int";

  kernelInfo["defines/" "p_G00ID"]= G00ID;
  kernelInfo["defines/" "p_G01ID"]= G01ID;
  kernelInfo["defines/" "p_G02ID"]= G02ID;
  kernelInfo["defines/" "p_G11ID"]= G11ID;
  kernelInfo["defines/" "p_G12ID"]= G12ID;
  kernelInfo["defines/" "p_G22ID"]= G22ID;
  kernelInfo["defines/" "p_GWJID"]= GWJID;
  
  char kernelName[BUFSIZ];
  switch(mesh.elementType){
  case TRIANGLES:
    sprintf(kernelName, "ellipticBuildOperatorMatrixContinuousTri2D");
    break;
  case QUADRILATERALS:
    sprintf(kernelName, "ellipticBuildOperatorMatrixContinuousQuad2D");
    break;
  case TETRAHEDRA:
    sprintf(kernelName, "ellipticBuildOperatorMatrixContinuousTet3D");
    break;
  case HEXAHEDRA:
    sprintf(kernelName, "ellipticBuildOperatorMatrixContinuousHex3D");
    break;
  }
  
  occa::kernel buildMatrixKernel =
    platform.buildKernel(DELLIPTIC "/okl/ellipticBuildOperatorMatrixContinuous.okl",
			 kernelName,
			 kernelInfo);

  // ---------------------------------------------------------------------
  // #1. create one-ring mesh 


  // A. build one ring mesh
  mesh.HaloRingSetup();

  // B. transfer into halo
  dlong Nel = mesh.Nelements;
  dlong ringNel = mesh.totalRingElements;
  dlong allNel = Nel+ringNel;

  // i. import ggeo
  size_t NggeoBlk;
  if(mesh.elementType == TRIANGLES ||
     mesh.elementType == TETRAHEDRA)
    NggeoBlk = mesh.Nggeo;
  else
    NggeoBlk = mesh.Nggeo*mesh.Np;
  
  occa::memory o_ringGgeo =
    platform.device.malloc(allNel*NggeoBlk*sizeof(dfloat));

  o_ringGgeo.copyFrom(mesh.o_ggeo, Nel*NggeoBlk*sizeof(dfloat), 0, 0);
					   
  mesh.ringHalo->Exchange(o_ringGgeo, NggeoBlk, ogs_dfloat);

  // ii. import globalNumbers
  size_t NnumBlk = mesh.Np;
  occa::memory o_ringMaskedGlobalNumbering =
    platform.device.malloc(allNel*NnumBlk*sizeof(hlong));

  occa::memory o_maskedGlobalNumbering = platform.device.malloc(Nel*NnumBlk*sizeof(hlong), maskedGlobalNumbering);
  
  o_ringMaskedGlobalNumbering.copyFrom(o_maskedGlobalNumbering, Nel*NnumBlk*sizeof(hlong), 0, 0);

  mesh.ringHalo->Exchange(o_ringMaskedGlobalNumbering, NnumBlk, ogs_hlong);

  hlong *ringMaskedGlobalNumbering = (hlong*) calloc(allNel*NnumBlk,sizeof(hlong));
  o_ringMaskedGlobalNumbering.copyTo(ringMaskedGlobalNumbering);

  // iii. build output map
  hlong Ngather = ogsMasked->Ngather;
  hlong *globalRowStarts = (hlong*) calloc(mesh.size+1,sizeof(hlong));
  MPI_Allgather(&Ngather, 1, MPI_HLONG, globalRowStarts+1, 1, MPI_HLONG, mesh.comm);
  for(int r=0;r<mesh.size;++r) {
    globalRowStarts[r+1] = globalRowStarts[r]+globalRowStarts[r+1];
  }

  hlong gatherStart = globalRowStarts[mesh.rank];
  hlong gatherEnd   = globalRowStarts[mesh.rank+1];

  // iv. for each element build map to shuffle output
  globalNode_t *gnodes = (globalNode_t*) calloc(mesh.Np, sizeof(globalNode_t));
  dlong    *colMap = (dlong*) calloc(mesh.Np*allNel, sizeof(dlong));
  for(dlong e=0;e<allNel;++e){
    for(int n=0;n<mesh.Np;++n){
      hlong id = ringMaskedGlobalNumbering[e*mesh.Np+n];
      gnodes[n].gnum = (id<0) ? BIG_NUM:id;
      gnodes[n].lnum = n;
    }
    qsort(gnodes, mesh.Np, sizeof(globalNode_t), compareGlobalNodes2);
    for(int n=0;n<mesh.Np;++n){
      colMap[e*mesh.Np+gnodes[n].lnum] = n; // order to map write index into
    }
  }

  occa::memory o_colMap = platform.device.malloc(allNel*mesh.Np*sizeof(dlong), colMap);

  // find where to put the row segments
  dlong *rowCount = (dlong*) calloc(mesh.Np*allNel, sizeof(dlong));
  for(dlong e=0;e<allNel;++e){
    for(int n=0;n<mesh.Np;++n){
      hlong id = ringMaskedGlobalNumbering[e*mesh.Np+n];
      if(gatherStart<=id && id<gatherEnd){
	++(rowCount[id-gatherStart]);
      }
    }
  }
  dlong *rowStarts = (dlong*) calloc(Ngather+1, sizeof(dlong));
  for(dlong n=0;n<Ngather;++n){
    rowStarts[n+1] = rowStarts[n] + rowCount[n];
  }
  
  dlong *rowMap = (dlong*) calloc(mesh.Np*allNel, sizeof(dlong));
  memset(rowMap, -1, mesh.Np*allNel*sizeof(dlong));
  for(dlong e=0;e<allNel;++e){
    for(int n=0;n<mesh.Np;++n){
      
      hlong id = ringMaskedGlobalNumbering[e*mesh.Np+n];
      if(gatherStart<=id && id<gatherEnd){
	rowMap[e*mesh.Np+n] = (rowStarts[id-gatherStart]++);
      }
    }
  }

  dlong maxDegree = 0;
  for(dlong g=0;g<Ngather;++g){
    dlong deg = (g==0) ? rowStarts[g] : rowStarts[g]-rowStarts[g-1];
    ///    printf("rowCount[%d] = %d\n", g,deg);
    maxDegree = mymax(maxDegree, deg);
  }
  printf("maxDegree=%d\n", maxDegree);
  
  kernelInfo["defines/" "MAX_DEGREE"]= (int) maxDegree;

  occa::kernel squeezeGapsKernel = 
    platform.buildKernel(DELLIPTIC "/okl/ellipticBuildOperatorMatrixContinuous.okl",
			 "ellipticSqueezeGaps",
			 kernelInfo);
  
  occa::kernel assembleMatrixKernel = 
    platform.buildKernel(DELLIPTIC "/okl/ellipticBuildOperatorMatrixContinuous.okl",
			 "ellipticAssembleEntries",
			 kernelInfo);
  
  occa::memory o_rowStarts = platform.device.malloc(Ngather*sizeof(dlong), rowStarts);

  dlong allNgather = rowStarts[Ngather-1];
  printf("gatherStart = %lld, gatherEnd = %lld, Ngather=%d, allNgather=%d\n",
	 gatherStart, gatherEnd, Ngather, allNgather*mesh.Np);
  
  occa::memory o_rowMap = platform.device.malloc(allNel*mesh.Np*sizeof(dlong), rowMap);
  
  o_A = platform.device.malloc(allNgather*mesh.Np*sizeof(nonZero_t));
  occa::memory o_rowCounts = platform.malloc(Ngather*sizeof(dlong));
  occa::memory o_AL2 = platform.malloc(allNgather*mesh.Np*sizeof(nonZero_t)); // wasteful until fixed

  platform.device.finish();
  double ticA = MPI_Wtime();

  

  
  switch(mesh.elementType){
  case TRIANGLES:
    buildMatrixKernel(allNel, o_ringMaskedGlobalNumbering,  o_rowMap, o_colMap,
		      mesh.o_S, mesh.o_MM, o_ringGgeo, 
		      lambda, o_A);

    break;
  case QUADRILATERALS:
    buildMatrixKernel(allNel, o_ringMaskedGlobalNumbering, o_rowMap, o_colMap,
		      mesh.o_D, o_ringGgeo,
		      lambda, o_A);

    break;
  case TETRAHEDRA:
    buildMatrixKernel(allNel, o_ringMaskedGlobalNumbering, o_rowMap, o_colMap,
		      mesh.o_S, mesh.o_MM, o_ringGgeo, 
		      lambda, o_A);
    break;
  case HEXAHEDRA:
    buildMatrixKernel(allNel, o_ringMaskedGlobalNumbering, o_rowMap, o_colMap,
		      mesh.o_D,  o_ringGgeo,
		      lambda, o_A);
    break;
  }

  occa::memory o_newCounts = platform.device.malloc((Ngather+1)*sizeof(dlong));
  assembleMatrixKernel(Ngather, o_rowStarts, o_A, o_AL2, o_newCounts);
  o_A.free();
  
  dlong *newCounts = (dlong*) calloc(Ngather, sizeof(dlong));
  o_newCounts.copyTo(newCounts);
  dlong *newStarts = (dlong*) calloc(Ngather+1, sizeof(dlong));
  for(int n=0;n<Ngather;++n){
    newStarts[n+1] = newStarts[n] + newCounts[n];
  }
  Annz = newStarts[Ngather];
  printf("Annz = %d\n", Annz);
  occa::memory o_newStarts = platform.device.malloc((Ngather+1)*sizeof(dlong), newStarts);
  
  o_A = platform.device.malloc(Annz*sizeof(nonZero_t));
  squeezeGapsKernel(Ngather, o_rowStarts, o_newStarts, o_AL2, o_A);
  platform.device.finish();
  double ticB = MPI_Wtime();

  printf("matrix build on device took %g\n", ticB-ticA);

  free(ringMaskedGlobalNumbering);
  free(globalRowStarts);
  free(gnodes);
  free(rowCount);
  free(rowStarts);
  free(rowMap);
  free(newCounts);
  free(newStarts);
}

  

