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

#ifndef OGS_KERNELS_HPP
#define OGS_KERNELS_HPP

#include <limits>
#include "ogs.hpp"
#include "ogsDefs.h"

#define DEFINE_ADD_OGS_INIT(T)                                  \
  static T init_##T##_add = (T)  0;                             \
  static T init_##T##_mul = (T)  1;                             \
  static T init_##T##_min = (T)  std::numeric_limits<T>::max(); \
  static T init_##T##_max = (T) -std::numeric_limits<T>::max();

namespace ogs {

extern int Nrefs;

extern occa::stream dataStream;

void initKernels(MPI_Comm& comm, occa::device& device);

void freeKernels();

//Setup a gslib struct
void *gsSetup(MPI_Comm meshComm,
              dlong NuniqueBases,
              hlong *gatherGlobalNodes,
              int nonsymm, int verbose);

void gsUnique(hlong *gatherGlobalNodes,
              dlong NuniqueBases,
              MPI_Comm meshComm);

void gsFree(void* gs);

#define DEFINE_GATHERSCATTER_KERNEL(T,OP) \
  extern occa::kernel gatherScatterKernel_##T##_##OP;

#define DEFINE_GATHER_KERNEL(T,OP) \
  extern occa::kernel gatherKernel_##T##_##OP;

#define DEFINE_SCATTER_KERNEL(T) \
  extern occa::kernel scatterKernel_##T;

#define DEFINE_HALO_KERNEL(T) \
  extern occa::kernel haloExtractKernel_##T;

#define DEFINE_KERNELS(T)                        \
  OGS_FOR_EACH_OP(T,DEFINE_GATHERSCATTER_KERNEL) \
  OGS_FOR_EACH_OP(T,DEFINE_GATHER_KERNEL)        \
  DEFINE_SCATTER_KERNEL(T)                       \
  DEFINE_HALO_KERNEL(T)

OGS_FOR_EACH_TYPE(DEFINE_KERNELS)

#undef DEFINE_GATHERSCATTER_KERNEL
#undef DEFINE_GATHER_KERNEL
#undef DEFINE_SCATTER_KERNEL
#undef DEFINE_HALO_KERNEL
#undef DEFINE_KERNELS

void occaGatherScatterStart(occa::memory& o_v,
                            const int Nentries, const int Nvectors, const dlong stride,
                            const ogs_type type, const ogs_op op, ogs_t &ogs);
void occaGatherScatterFinish(occa::memory& o_v,
                            const int Nentries, const int Nvectors, const dlong stride,
                            const ogs_type type, const ogs_op op, ogs_t &ogs);

void occaGatherStart(occa::memory& o_gv, occa::memory& o_v,
                     const int Nentries, const int Nvectors,
                     const dlong gstride, const dlong stride,
                     const ogs_type type, const ogs_op op, ogs_t &ogs);
void occaGatherFinish(occa::memory& o_gv, occa::memory& o_v,
                     const int Nentries, const int Nvectors,
                     const dlong gstride, const dlong stride,
                     const ogs_type type, const ogs_op op, ogs_t &ogs);

void occaScatterStart(occa::memory& o_v, occa::memory& o_gv,
                      const int Nentries, const int Nvectors,
                      const dlong stride, const dlong gstride,
                      const ogs_type type, const ogs_op op, ogs_t &ogs);
void occaScatterFinish(occa::memory& o_v, occa::memory& o_gv,
                      const int Nentries, const int Nvectors,
                      const dlong stride, const dlong gstride,
                      const ogs_type type, const ogs_op op, ogs_t &ogs);

void hostGatherScatter(void* v, const int Nentries, const int Nvectors,
                       const dlong stride, const ogs_type type,
                       const ogs_op op, ogs_t &ogs);

void hostGather(void* gv, void* v, const int Nentries, const int Nvectors,
                const dlong gstride, const dlong stride,
                const ogs_type type, const ogs_op op, ogs_t &ogs);

void hostScatter(void* v, void* gv, const int Nentries, const int Nvectors,
                 const dlong stride, const dlong gstride,
                 const ogs_type type, const ogs_op op, ogs_t &ogs);

void occaGatherScatterKernel(const dlong Ngather,
                             const int Nentries,
                             const int Nvectors,
                             const dlong stride,
                             occa::memory& o_gatherStarts,
                             occa::memory& o_gatherIds,
                             const ogs_type type,
                             const ogs_op op,
                             occa::memory&  o_v);

void occaGatherKernel(const dlong Ngather,
                      const int Nentries,
                      const int Nvectors,
                      const dlong stride,
                      const dlong gtride,
                      occa::memory& o_gatherStarts,
                      occa::memory& o_gatherIds,
                      const ogs_type type,
                      const ogs_op op,
                      occa::memory& o_v,
                      occa::memory& o_gv);

void occaScatterKernel(const dlong Ngather,
                       const int Nentries,
                       const int Nvectors,
                       const dlong gtride,
                       const dlong stride,
                       occa::memory& o_gatherStarts,
                       occa::memory& o_gatherIds,
                       const ogs_type type,
                       const ogs_op op,
                       occa::memory& o_gv,
                       occa::memory& o_v);

void hostGatherScatterKernel(const dlong Ngather,
                             const int Nentries,
                             const int Nvectors,
                             const dlong stride,
                             dlong* gatherStarts,
                             dlong* gatherIds,
                             const ogs_type type,
                             const ogs_op op,
                             void* v);

void hostGatherKernel(const dlong Ngather,
                      const int Nentries,
                      const int Nvectors,
                      const dlong stride,
                      const dlong gstride,
                      const dlong *gatherStarts,
                      const dlong *gatherIds,
                      const ogs_type type,
                      const ogs_op op,
                      const void *v,
                      void *gv);

void hostScatterKernel(const dlong Ngather,
                       const int Nentries,
                       const int Nvectors,
                       const dlong gstride,
                       const dlong stride,
                       const dlong *gatherStarts,
                       const dlong *gatherIds,
                       const ogs_type type,
                       const ogs_op op,
                       const void *gv,
                       void *v);

void gsGatherScatter(void* v, const int Nentries, const int Nvectors,
                     const dlong stride, const ogs_type type, const ogs_op op,
                     void *gsh);
void gsGather(void* v, const int Nentries, const int Nvectors,
              const dlong stride, const ogs_type type, const ogs_op op,
              void *gsh);
void gsScatter(void* v, const int Nentries, const int Nvectors,
              const dlong stride, const ogs_type type, const ogs_op op,
              void *gsh);

void occaExtractKernel(const dlong Nsend,
                       const int Nentries,
                       occa::memory& o_haloSendIds,
                       const ogs_type type,
                       occa::memory& o_v,
                       occa::memory& o_haloBuf);

void hostExtractKernel(const dlong Nsend,
                       const int Nentries,
                       dlong *haloSendIds,
                       const ogs_type type,
                       void *v,
                       void *hostBuf);
}

#endif
