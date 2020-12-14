/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"

namespace parAlmond {

int Nrefs = 0;

//NC: Hard code these for now. Should be sufficient for GPU devices, but needs attention for CPU
const int blockSize = 256;
int NonzerosPerBlock = 1024; //should be a multiple of blockSize for good unrolling

occa::kernel SpMVcsrKernel1;
occa::kernel SpMVcsrKernel2;
occa::kernel SpMVmcsrKernel;

occa::kernel kcycleCombinedOp1Kernel;
occa::kernel kcycleCombinedOp2Kernel;
occa::kernel kcycleWeightedCombinedOp1Kernel;
occa::kernel kcycleWeightedCombinedOp2Kernel;
occa::kernel vectorAddInnerProdKernel;
occa::kernel vectorAddWeightedInnerProdKernel;

MPI_Datatype MPI_NONZERO_T;

void buildParAlmondKernels(platform_t& platform){

  // Make the MPI_NONZERO_T data type
  parCOO::nonZero_t NZ;
  MPI_Datatype dtype[3] = {MPI_HLONG, MPI_HLONG, MPI_DFLOAT};
  int blength[3] = {1, 1, 1};
  MPI_Aint addr[3], displ[3];
  MPI_Get_address ( &(NZ.row), addr+0);
  MPI_Get_address ( &(NZ.col), addr+1);
  MPI_Get_address ( &(NZ.val), addr+2);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  MPI_Type_create_struct (3, blength, displ, dtype, &MPI_NONZERO_T);
  MPI_Type_commit (&MPI_NONZERO_T);

  //seed rng
  int rank=platform.rank;
  double seed = (double) rank;
  srand48(seed);

  //build kernels
  occa::properties kernelInfo = platform.props;

  kernelInfo["defines/" "p_BLOCKSIZE"]= blockSize;
  kernelInfo["defines/" "p_NonzerosPerBlock"]= NonzerosPerBlock;

  if (rank==0) {printf("Compiling parALMOND Kernels...");fflush(stdout);}

  SpMVcsrKernel1  = platform.buildKernel(PARALMOND_DIR"/okl/SpMVcsr.okl",  "SpMVcsr1",  kernelInfo);
  SpMVcsrKernel2  = platform.buildKernel(PARALMOND_DIR"/okl/SpMVcsr.okl",  "SpMVcsr2",  kernelInfo);
  SpMVmcsrKernel  = platform.buildKernel(PARALMOND_DIR"/okl/SpMVmcsr.okl", "SpMVmcsr1", kernelInfo);

  vectorAddInnerProdKernel = platform.buildKernel(PARALMOND_DIR"/okl/vectorAddInnerProd.okl", "vectorAddInnerProd", kernelInfo);
  vectorAddWeightedInnerProdKernel = platform.buildKernel(PARALMOND_DIR"/okl/vectorAddInnerProd.okl", "vectorAddWeightedInnerProd", kernelInfo);

  kcycleCombinedOp1Kernel = platform.buildKernel(PARALMOND_DIR"/okl/kcycleCombinedOp.okl", "kcycleCombinedOp1", kernelInfo);
  kcycleCombinedOp2Kernel = platform.buildKernel(PARALMOND_DIR"/okl/kcycleCombinedOp.okl", "kcycleCombinedOp2", kernelInfo);
  kcycleWeightedCombinedOp1Kernel = platform.buildKernel(PARALMOND_DIR"/okl/kcycleCombinedOp.okl", "kcycleWeightedCombinedOp1", kernelInfo);
  kcycleWeightedCombinedOp2Kernel = platform.buildKernel(PARALMOND_DIR"/okl/kcycleCombinedOp.okl", "kcycleWeightedCombinedOp2", kernelInfo);

  if(rank==0) printf("done.\n");
}

void freeParAlmondKernels() {

  SpMVcsrKernel1.free();
  SpMVcsrKernel2.free();
  SpMVmcsrKernel.free();

  kcycleCombinedOp1Kernel.free();
  kcycleCombinedOp2Kernel.free();
  kcycleWeightedCombinedOp1Kernel.free();
  kcycleWeightedCombinedOp2Kernel.free();
  vectorAddInnerProdKernel.free();
  vectorAddWeightedInnerProdKernel.free();
}


} //namespace parAlmond
