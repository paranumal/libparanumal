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

static void compressMatrixMultiDevice(platform_t &platform, mesh_t &mesh, ogs_t *ogsMasked, int includeLast,
				      occa::memory &o_AL,
				      dlong nnzLocal,
				      deviceSort_t &sorter,
				      deviceScan_t &scanner,
				      occa::memory &o_A,
				      dlong &Annz,
				      parAlmond::parCOO& A){


  // 1. sort based on row (fastest) then column in each row
  sorter.sort(nnzLocal, o_AL);

  // 2. compactify
  occa::memory o_AL2;
  nnzLocal = scanner.trashCompactor(platform, nnzLocal, sizeof(nonZero_t), includeLast, o_AL, o_AL2);
  o_AL.free();
  
  // 3. copy to host
  nonZero_t *h_AL = (nonZero_t*) calloc(nnzLocal, sizeof(nonZero_t));
  o_AL2.copyTo(h_AL);
  o_AL2.free();
  
  // number of degrees of freedom on this rank (after gathering)
  hlong Ngather = ogsMasked->Ngather;

  // every gathered degree of freedom has its own global id
  A.globalRowStarts = (hlong*) calloc(mesh.size+1,sizeof(hlong));
  A.globalColStarts = (hlong*) calloc(mesh.size+1,sizeof(hlong));
  MPI_Allgather(&Ngather, 1, MPI_HLONG, A.globalRowStarts+1, 1, MPI_HLONG, mesh.comm);
  for(int r=0;r<mesh.size;++r) {
    A.globalRowStarts[r+1] = A.globalRowStarts[r]+A.globalRowStarts[r+1];
    A.globalColStarts[r+1] = A.globalRowStarts[r+1];
  }
  
  int *AsendCounts  = (int*) calloc(mesh.size, sizeof(int));
  int *ArecvCounts  = (int*) calloc(mesh.size, sizeof(int));
  int *AsendOffsets = (int*) calloc(mesh.size+1, sizeof(int));
  int *ArecvOffsets = (int*) calloc(mesh.size+1, sizeof(int));
  
  // count how many non-zeros to send to each process
  int rr=0;
  for(dlong n=0;n<nnzLocal;++n) {
    const hlong id = h_AL[n].row;
    while(id>=A.globalRowStarts[rr+1]) rr++;
    AsendCounts[rr]++;
  }

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(AsendCounts, 1, MPI_INT, ArecvCounts, 1, MPI_INT, mesh.comm);

  // find send and recv offsets for gather
  A.nnz = 0;
  for(int r=0;r<mesh.size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    A.nnz += ArecvCounts[r];
  }

  nonZero_t *entriesIn = (nonZero_t*) calloc(A.nnz, sizeof(nonZero_t));

  // determine number to receive
  MPI_Alltoallv(h_AL,      AsendCounts, AsendOffsets, parAlmond::MPI_NONZERO_T,
		entriesIn, ArecvCounts, ArecvOffsets, parAlmond::MPI_NONZERO_T,
		mesh.comm);

  // 1. load onto device and sort
  occa::memory o_AL3 = platform.device.malloc(A.nnz*sizeof(nonZero_t), entriesIn);
  sorter.sort(A.nnz, o_AL3);

  // 2. free up host storage
  free(entriesIn);
  
  // 3. assemble matrix (gather duplicates)
  occa::memory o_AL4;
  includeLast = 1; // the bc nodes are already removed
  Annz = scanner.trashCompactor(platform, A.nnz, sizeof(nonZero_t), includeLast, o_AL3, o_A);
  o_AL3.free();

  // 4. build final host storage of compressed matrix
  //  A.entries = (nonZero_t*) calloc(A.nnz, sizeof(nonZero_t));
  //  o_AL4.copyTo(A.entries);

  // release buffers
  free(AsendCounts);
  free(ArecvCounts);
  free(AsendOffsets);
  free(ArecvOffsets);
  free(h_AL);
}


void elliptic_t::BuildOperatorMatrixContinuousDevice(occa::kernel &buildMatrixKernel, occa::memory &o_maskedGlobalNumbering, hlong BIG_NUM,
						     deviceSort_t &sorter, deviceScan_t &scanner, parAlmond::parCOO& A, occa::memory &o_A, dlong &Annz) {

  dlong nnzLocal = mesh.Np*mesh.Np*mesh.Nelements+1;


  occa::memory o_AL =  platform.malloc(nnzLocal*sizeof(nonZero_t));

  // make sure at least one fake entry is in local matrix
  nonZero_t dummyEnt;
  dummyEnt.row = BIG_NUM;
  dummyEnt.col = BIG_NUM;
  dummyEnt.val = 0;
  
  o_AL.copyFrom(&dummyEnt, sizeof(nonZero_t), (nnzLocal-1)*sizeof(nonZero_t));
  
  switch(mesh.elementType){
  case TRIANGLES:
    buildMatrixKernel(mesh.Nelements, o_maskedGlobalNumbering,
		      mesh.o_S, mesh.o_MM, mesh.o_ggeo,
		      lambda, o_AL);

    break;
  case QUADRILATERALS:
    buildMatrixKernel(mesh.Nelements, o_maskedGlobalNumbering,
		      mesh.o_D, mesh.o_ggeo,
		      lambda, o_AL);

    break;
  case TETRAHEDRA:
    buildMatrixKernel(mesh.Nelements, o_maskedGlobalNumbering,
		      mesh.o_S, mesh.o_MM, mesh.o_ggeo, lambda, o_AL);
    break;
  case HEXAHEDRA:
    buildMatrixKernel(mesh.Nelements, o_maskedGlobalNumbering,
		      mesh.o_D,  mesh.o_ggeo,
		      lambda, o_AL);
    break;
  }

  int includeLast = 0; 
  
  // assemble on device + MPI
  // [ warning - destroys o_AL ]
  // output is to host A
  compressMatrixMultiDevice(platform, mesh, ogsMasked, includeLast, o_AL, nnzLocal, sorter, scanner, o_A, Annz, A);
  
}
