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
#include "ogsKernels.hpp"

namespace ogs {

OGS_DEFINE_TYPE_SIZES()
OGS_FOR_EACH_TYPE(DEFINE_ADD_OGS_INIT)

void hostGather(void* gv,
                void* v,
                const int Nentries,
                const int Nvectors,
                const dlong gstride,
                const dlong stride,
                const ogs_type type,
                const ogs_op op,
                const ogs_transpose trans,
                ogs_t &ogs){

  const size_t Nbytes = ogs_type_size[type];

  if (trans == ogs_sym)
    LIBP_ABORT(string("Calling ogs::Gather in ogs_sym mode not supported."))

  ogs.reallocHostBuffer(Nbytes*Nentries*Nvectors);

  dlong NhaloGather = (trans == ogs_notrans) ? ogs.NhaloGather : ogs.NhaloScatter;

  // gather halo nodes
  if (NhaloGather) {
    if (trans == ogs_notrans)
      hostGatherKernel(ogs.NhaloGather, Nentries, Nvectors, stride, ogs.Nhalo,
                       ogs.haloGatherOffsets, ogs.haloGatherIds,
                       type, op, v, ogs.hostBuf);
    else
      hostGatherKernel(ogs.NhaloScatter, Nentries, Nvectors, stride, ogs.Nhalo,
                       ogs.haloScatterOffsets, ogs.haloScatterIds,
                       type, op, v, ogs.hostBuf);
  }

  // MPI based gather using libgs
  gsGatherScatter(ogs.hostBuf, Nentries, Nvectors, ogs.Nhalo,
                  type, op, trans, ogs.gsh);

  if (ogs.NhaloGather)
    for (int i=0;i<Nvectors;i++)
      memcpy((char*)gv+ogs.NlocalGather*Nbytes*Nentries + gstride*Nbytes*i,
             (char*)ogs.hostBuf+ogs.Nhalo*Nbytes*Nentries*i,
             ogs.NhaloGather*Nentries*Nbytes);

  // gather interior nodes
  if (ogs.Nlocal) {
    if (trans == ogs_notrans)
      hostGatherKernel(ogs.NlocalGather, Nentries, Nvectors, stride, gstride,
                       ogs.localGatherOffsets, ogs.localGatherIds,
                       type, op, v, gv);
    else
      hostGatherKernel(ogs.NlocalScatter, Nentries, Nvectors, stride, gstride,
                       ogs.localScatterOffsets, ogs.localScatterIds,
                       type, op, v, gv);
  }
}

/*------------------------------------------------------------------------------
  The basic gather kernel
------------------------------------------------------------------------------*/
#define DEFINE_GATHER(T,OP)                                                     \
static void hostGatherKernel_##T##_##OP(const dlong N,                          \
                                        const int   Nentries,                   \
                                        const int   Nvectors,                   \
                                        const dlong stride,                     \
                                        const dlong gstride,                    \
                                        const dlong *gatherStarts,              \
                                        const dlong *gatherIds,                 \
                                        const     T *q,                         \
                                                  T *gatherq)                   \
{                                                                               \
  for(dlong n=0;n<N*Nentries*Nvectors;++n){                                     \
    const int m     = n/(N*Nentries);                                           \
    const dlong vid = n%(N*Nentries);                                           \
    const dlong gid = vid/Nentries;                                             \
    const int k     = vid%Nentries;                                             \
    const dlong start = gatherStarts[gid];                                      \
    const dlong end = gatherStarts[gid+1];                                      \
    T gq = init_##T##_##OP;                                                     \
    for(dlong g=start;g<end;++g){                                               \
      const dlong id = gatherIds[g];                                            \
      OGS_DO_##OP(gq,q[k+id*Nentries+m*stride]);                                \
    }                                                                           \
    gatherq[k+gid*Nentries+m*gstride] = gq;                                     \
  }                                                                             \
}

#define DEFINE_PROCS(T) \
  OGS_FOR_EACH_OP(T,DEFINE_GATHER)

OGS_FOR_EACH_TYPE(DEFINE_PROCS)

#define SWITCH_TYPE_CASE(T) case ogs_##T: { WITH_TYPE(T); break; }
#define SWITCH_TYPE(type) switch(type) { \
    OGS_FOR_EACH_TYPE(SWITCH_TYPE_CASE) case ogs_type_n: break; }

#define SWITCH_OP_CASE(T,OP) case ogs_##OP: { WITH_OP(T,OP); break; }
#define SWITCH_OP(T,op) switch(op) { \
    OGS_FOR_EACH_OP(T,SWITCH_OP_CASE) case ogs_op_n: break; }


void hostGatherKernel(const dlong N,
                      const int Nentries,
                      const int Nvectors,
                      const dlong stride,
                      const dlong gstride,
                      const dlong *gatherStarts,
                      const dlong *gatherIds,
                      const ogs_type type,
                      const ogs_op op,
                      const void *v,
                      void *gv) {

#define WITH_OP(T,OP)                         \
  hostGatherKernel_##T##_##OP(N,              \
                              Nentries,       \
                              Nvectors,       \
                              stride,         \
                              gstride,        \
                              gatherStarts,   \
                              gatherIds,      \
                              (T*) v,         \
                              (T*) gv);
#define WITH_TYPE(T) SWITCH_OP(T,op)

  SWITCH_TYPE(type)

#undef  WITH_TYPE
#undef  WITH_OP
}

} //namespace ogs