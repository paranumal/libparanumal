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


#include "core.hpp"
#include "ogsKernels.hpp"

//convert a macro command into a string
#define _STR(x) #x
#define STR(x) _STR(x)

namespace ogs {

  int Nrefs = 0;

  occa::stream dataStream;

#define DEFINE_GATHERSCATTER_KERNEL(T,OP) \
  occa::kernel gatherScatterKernel_##T##_##OP;

#define DEFINE_GATHER_KERNEL(T,OP) \
  occa::kernel gatherKernel_##T##_##OP;

#define DEFINE_SCATTER_KERNEL(T) \
  occa::kernel scatterKernel_##T;

#define DEFINE_HALO_KERNEL(T) \
  occa::kernel haloExtractKernel_##T;

#define DEFINE_KERNELS(T)                        \
  OGS_FOR_EACH_OP(T,DEFINE_GATHERSCATTER_KERNEL) \
  OGS_FOR_EACH_OP(T,DEFINE_GATHER_KERNEL)        \
  DEFINE_SCATTER_KERNEL(T)                       \
  DEFINE_HALO_KERNEL(T)

  OGS_FOR_EACH_TYPE(DEFINE_KERNELS)


void initKernels(MPI_Comm& comm, occa::device& device) {

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  dataStream = device.createStream();

  occa::properties kernelInfo;
  occaDeviceProperties(device, kernelInfo);

#define DEFINE_OCCA_ADD_INIT(T) \
  kernelInfo["defines/init_" STR(T) "_add"] = (T)  0;                              \
  kernelInfo["defines/init_" STR(T) "_mul"] = (T)  1;                              \
  kernelInfo["defines/init_" STR(T) "_min"] = (T)  std::numeric_limits<T>::max(); \
  kernelInfo["defines/init_" STR(T) "_max"] = (T) -std::numeric_limits<T>::max();

//OCCA properties don't have an operator+ for long long int, so alias it to int64_t
typedef int64_t long_long;
  OGS_FOR_EACH_TYPE(DEFINE_OCCA_ADD_INIT)

  kernelInfo["includes"] += DOGS "/include/ogsDefs.h";

  if (rank==0) {printf("Compiling GatherScatter Kernels...");fflush(stdout);}

#define DEFINE_GATHERSCATTER_BUILD(T,OP)                                           \
  gatherScatterKernel_##T##_##OP = buildKernel(device,                             \
                                             DOGS "/okl/gatherScatter.okl",        \
                                             "gatherScatter_" STR(T) "_" STR(OP),  \
                                             kernelInfo, comm);                    \

#define DEFINE_GATHER_BUILD(T,OP)                                                  \
  gatherKernel_##T##_##OP = buildKernel(device,                                    \
                                             DOGS "/okl/gatherScatter.okl",        \
                                             "gather_" STR(T) "_" STR(OP),         \
                                             kernelInfo, comm);                    \

#define DEFINE_SCATTER_BUILD(T)                                                    \
  scatterKernel_##T      = buildKernel(device,                                     \
                                             DOGS "/okl/gatherScatter.okl",        \
                                             "scatter_" STR(T),                    \
                                             kernelInfo, comm);                    \

#define DEFINE_HALO_BUILD(T)                                                       \
  haloExtractKernel_##T  = buildKernel(device,                                     \
                                             DOGS "/okl/haloExtract.okl",          \
                                             "haloExtract_" STR(T),                \
                                             kernelInfo, comm);                    \

#define DEFINE_BUILD(T)                         \
  OGS_FOR_EACH_OP(T,DEFINE_GATHERSCATTER_BUILD) \
  OGS_FOR_EACH_OP(T,DEFINE_GATHER_BUILD)        \
  DEFINE_SCATTER_BUILD(T)                       \
  DEFINE_HALO_BUILD(T)

  OGS_FOR_EACH_TYPE(DEFINE_BUILD)

  if(rank==0) printf("done.\n");
}

void freeKernels() {

#define DEFINE_GATHERSCATTER_FREE(T,OP)      \
  gatherScatterKernel_##T##_##OP.free();

#define DEFINE_GATHER_FREE(T,OP)      \
  gatherKernel_##T##_##OP.free();

#define DEFINE_SCATTER_FREE(T)       \
  scatterKernel_##T.free();

#define DEFINE_HALO_FREE(T)          \
  haloExtractKernel_##T.free();

#define DEFINE_FREE(T)                         \
  OGS_FOR_EACH_OP(T,DEFINE_GATHERSCATTER_FREE) \
  OGS_FOR_EACH_OP(T,DEFINE_GATHER_FREE)        \
  DEFINE_SCATTER_FREE(T)                       \
  DEFINE_HALO_FREE(T)

  OGS_FOR_EACH_TYPE(DEFINE_FREE)
}

} //namespace ogs

