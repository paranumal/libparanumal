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

void hostGatherScatter(void* v,
                       const int Nentries,
                       const int Nvectors,
                       const dlong stride,
                       const ogs_type type,
                       const ogs_op op,
                       const ogs_transpose trans,
                       ogs_t &ogs){

  const size_t Nbytes = ogs_type_size[type];

  ogs.reallocHostBuffer(Nbytes*Nentries*Nvectors);

  dlong NhaloGather  = (trans == ogs_notrans) ? ogs.NhaloGather : ogs.NhaloScatter;
  dlong NhaloScatter = (trans == ogs_trans)   ? ogs.NhaloGather : ogs.NhaloScatter;

  void* gsh = (trans == ogs_sym) ? ogs.gshSym : ogs.gsh;

  // gather-scatter halo nodes
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

  // MPI based gather scatter using libgs
  gsGatherScatter(ogs.hostBuf, Nentries, Nvectors, ogs.Nhalo,
                  type, op, trans, gsh);

  if (NhaloScatter) {
    if (trans == ogs_trans)
      hostScatterKernel(ogs.NhaloGather, Nentries, Nvectors, stride, ogs.Nhalo,
                        ogs.haloGatherOffsets, ogs.haloGatherIds,
                        type, op, ogs.hostBuf, v);
    else
      hostScatterKernel(ogs.NhaloScatter, Nentries, Nvectors, stride, ogs.Nhalo,
                        ogs.haloScatterOffsets, ogs.haloScatterIds,
                        type, op, ogs.hostBuf, v);
  }

  if (trans == ogs_notrans) {
    if (ogs.NlocalFused)
      hostGatherScatterKernel(ogs.NlocalFused, Nentries, Nvectors, stride,
                              ogs.localFusedGatherOffsets,  ogs.localFusedGatherIds,
                              ogs.localFusedScatterOffsets, ogs.localFusedScatterIds,
                              type, op, v);
  } else if (trans == ogs_trans) {
    if (ogs.NlocalFused)
      hostGatherScatterKernel(ogs.NlocalFused, Nentries, Nvectors, stride,
                              ogs.localFusedScatterOffsets, ogs.localFusedScatterIds,
                              ogs.localFusedGatherOffsets,  ogs.localFusedGatherIds,
                              type, op, v);
  } else { //ogs_sym
    if (ogs.NlocalFusedSym)
      hostGatherScatterKernel(ogs.NlocalFusedSym, Nentries, Nvectors, stride,
                              ogs.localFusedOffsets, ogs.localFusedIds,
                              ogs.localFusedOffsets, ogs.localFusedIds,
                              type, op, v);
  }
}

/*------------------------------------------------------------------------------
  The basic gatherScatter kernel
------------------------------------------------------------------------------*/
#define DEFINE_GATHERSCATTER(T,OP)                                              \
static void hostGatherScatterKernel_##T##_##OP(const dlong N,                   \
                                               const int   Nentries,            \
                                               const int   Nvectors,            \
                                               const dlong stride,              \
                                               const dlong *gatherStarts,       \
                                               const dlong *gatherIds,          \
                                               const dlong *scatterStarts,      \
                                               const dlong *scatterIds,         \
                                                         T *q)                  \
{                                                                               \
  for(dlong n=0;n<N*Nentries*Nvectors;++n){                                     \
    const int m     = n/(N*Nentries);                                           \
    const dlong vid = n%(N*Nentries);                                           \
    const dlong gid = vid/Nentries;                                             \
    const int k     = vid%Nentries;                                             \
    const dlong gstart = gatherStarts[gid];                                     \
    const dlong gend = gatherStarts[gid+1];                                     \
    T gq = init_##T##_##OP;                                                     \
    for(dlong g=gstart;g<gend;++g){                                             \
      const dlong id = gatherIds[g];                                            \
      OGS_DO_##OP(gq,q[k+id*Nentries+m*stride]);                                \
    }                                                                           \
    const dlong sstart = scatterStarts[gid];                                    \
    const dlong send = scatterStarts[gid+1];                                    \
    for(dlong g=sstart;g<send;++g){                                             \
      const dlong id = gatherIds[g];                                            \
      q[k+id*Nentries+m*stride] = gq;                                           \
    }                                                                           \
  }                                                                             \
}

#define DEFINE_PROCS(T) \
  OGS_FOR_EACH_OP(T,DEFINE_GATHERSCATTER)

OGS_FOR_EACH_TYPE(DEFINE_PROCS)

#define SWITCH_TYPE_CASE(T) case ogs_##T: { WITH_TYPE(T); break; }
#define SWITCH_TYPE(type) switch(type) { \
    OGS_FOR_EACH_TYPE(SWITCH_TYPE_CASE) case ogs_type_n: break; }

#define SWITCH_OP_CASE(T,OP) case ogs_##OP: { WITH_OP(T,OP); break; }
#define SWITCH_OP(T,op) switch(op) { \
    OGS_FOR_EACH_OP(T,SWITCH_OP_CASE) case ogs_op_n: break; }


void hostGatherScatterKernel(const dlong N,
                             const int Nentries,
                             const int Nvectors,
                             const dlong stride,
                             dlong* gatherStarts,
                             dlong* gatherIds,
                             dlong* scatterStarts,
                             dlong* scatterIds,
                             const ogs_type type,
                             const ogs_op op,
                             void* v) {

#define WITH_OP(T,OP)                                \
  hostGatherScatterKernel_##T##_##OP(N,              \
                                     Nentries,       \
                                     Nvectors,       \
                                     stride,         \
                                     gatherStarts,   \
                                     gatherIds,      \
                                     scatterStarts,  \
                                     scatterIds,     \
                                     (T*)v);
#define WITH_TYPE(T) SWITCH_OP(T,op)

  SWITCH_TYPE(type)

#undef  WITH_TYPE
#undef  WITH_OP
}

} //namespace ogs