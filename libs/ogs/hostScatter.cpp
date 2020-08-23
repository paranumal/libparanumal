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

#include "ogs.hpp"
#include "ogs/ogsKernels.hpp"

namespace ogs {

OGS_DEFINE_TYPE_SIZES()

void hostScatter(void* v,
                 void* gv,
                 const int Nentries,
                 const int Nvectors,
                 const dlong stride,
                 const dlong gstride,
                 const ogs_type type,
                 const ogs_op op,
                 const ogs_transpose trans,
                 ogs_t &ogs){

  const size_t Nbytes = ogs_type_size[type];

  ogs.reallocHostBuffer(Nbytes*Nentries*Nvectors);

  if (trans == ogs_sym)
    LIBP_ABORT(string("Calling ogs::Scatter in ogs_sym mode not supported."))

  if (ogs.haloGather.Nrows)
    for (int i=0;i<Nvectors;i++)
      memcpy((char*)ogs.hostBuf+ogs.Nhalo*Nbytes*Nentries*i,
             (char*)gv+ogs.localGather.Nrows*Nbytes + gstride*Nbytes*i,
             ogs.haloGather.Nrows*Nbytes*Nentries);

  // MPI based scatter using gslib
  // (must use ogs_notrans so the negative ids don't contribute to op)
  gsGatherScatter(ogs.hostBuf, Nentries, Nvectors, ogs.Nhalo,
                  type, op, ogs_notrans, ogs.gsh);

  dlong NhaloScatter = (trans == ogs_trans) ? ogs.haloGather.Nrows : ogs.haloScatter.Nrows;

  if (NhaloScatter) {
    if (trans == ogs_trans)
      hostScatterKernel(ogs.haloGather.Nrows, Nentries, Nvectors, ogs.Nhalo, stride,
                        ogs.haloGather.rowStarts, ogs.haloGather.colIds,
                        type, op, ogs.hostBuf, v);
    else
      hostScatterKernel(ogs.haloScatter.Nrows, Nentries, Nvectors, ogs.Nhalo, stride,
                        ogs.haloScatter.rowStarts, ogs.haloScatter.colIds,
                        type, op, ogs.hostBuf, v);
  }

  // scatter interior nodes
  if (ogs.Nlocal) {
    if (trans == ogs_trans)
      hostScatterKernel(ogs.localGather.Nrows, Nentries, Nvectors, gstride, stride,
                        ogs.localGather.rowStarts, ogs.localGather.colIds,
                        type, op, gv, v);
    else
      hostScatterKernel(ogs.localScatter.Nrows, Nentries, Nvectors, gstride, stride,
                        ogs.localScatter.rowStarts, ogs.localScatter.colIds,
                        type, op, gv, v);
  }
}

/*------------------------------------------------------------------------------
  The basic gather kernel
------------------------------------------------------------------------------*/
#define DEFINE_SCATTER(T)                                                       \
static void hostScatterKernel_##T(const dlong N,                                \
                                  const int   Nentries,                         \
                                  const int   Nvectors,                         \
                                  const dlong gstride,                          \
                                  const dlong stride,                           \
                                  const dlong *scatterStarts,                   \
                                  const dlong *scatterIds,                      \
                                  const     T *gatherq,                         \
                                            T *q)                               \
{                                                                               \
  for(dlong n=0;n<N*Nentries*Nvectors;++n){                                     \
    const int m     = n/(N*Nentries);                                           \
    const dlong vid = n%(N*Nentries);                                           \
    const dlong gid = vid/Nentries;                                             \
    const int k     = vid%Nentries;                                             \
    const T gq = gatherq[k+gid*Nentries+m*gstride];                             \
    const dlong start = scatterStarts[gid];                                     \
    const dlong end = scatterStarts[gid+1];                                     \
    for(dlong g=start;g<end;++g){                                               \
      const dlong id = scatterIds[g];                                           \
      q[k+id*Nentries+m*stride] = gq;                                           \
    }                                                                           \
  }                                                                             \
}

#define DEFINE_PROCS(T) \
  DEFINE_SCATTER(T)

OGS_FOR_EACH_TYPE(DEFINE_PROCS)

#define SWITCH_TYPE_CASE(T) case ogs_##T: { WITH_TYPE(T); break; }
#define SWITCH_TYPE(type) switch(type) { \
    OGS_FOR_EACH_TYPE(SWITCH_TYPE_CASE) case ogs_type_n: break; }

void hostScatterKernel(const dlong N,
                       const int Nentries,
                       const int Nvectors,
                       const dlong gstride,
                       const dlong stride,
                       const dlong *scatterStarts,
                       const dlong *scatterIds,
                       const ogs_type type,
                       const ogs_op op,
                       const void *gv,
                       void *v) {

#define WITH_TYPE(T)                          \
  hostScatterKernel_##T(N,                    \
                        Nentries,             \
                        Nvectors,             \
                        gstride,              \
                        stride,               \
                        scatterStarts,        \
                        scatterIds,           \
                        (T*) gv,              \
                        (T*) v);

  SWITCH_TYPE(type)

#undef  WITH_TYPE
}

} //namespace ogs