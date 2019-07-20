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

void occaGatherScatterStart(occa::memory& o_v,
                            const int Nentries,
                            const int Nvectors,
                            const dlong stride,
                            const ogs_type type,
                            const ogs_op op,
                            const ogs_transpose trans,
                            ogs_t &ogs){

  const size_t Nbytes = ogs_type_size[type];

  ogs.reallocOccaBuffer(Nbytes*Nentries*Nvectors);

  dlong NhaloGather = (trans == ogs_notrans) ? ogs.NhaloGather : ogs.NhaloScatter;

  // gather halo nodes on device
  if (NhaloGather) {
    if (trans == ogs_notrans)
      occaGatherKernel(ogs.NhaloGather, Nentries, Nvectors, stride, ogs.Nhalo,
                       ogs.o_haloGatherOffsets, ogs.o_haloGatherIds,
                       type, op, o_v, ogs.o_haloBuf);
    else
      occaGatherKernel(ogs.NhaloScatter, Nentries, Nvectors, stride, ogs.Nhalo,
                       ogs.o_haloScatterOffsets, ogs.o_haloScatterIds,
                       type, op, o_v, ogs.o_haloBuf);

    ogs.device.finish();
    occa::stream currentStream = ogs.device.getStream();
    ogs.device.setStream(dataStream);

    for (int i=0;i<Nvectors;i++)
      ogs.o_haloBuf.copyTo((char*)ogs.haloBuf + ogs.Nhalo*Nbytes*Nentries*i,
                           NhaloGather*Nbytes*Nentries,
                           ogs.Nhalo*Nbytes*Nentries*i,
                           "async: true");

    ogs.device.setStream(currentStream);
  }
}


void occaGatherScatterFinish(occa::memory& o_v,
                             const int Nentries,
                             const int Nvectors,
                             const dlong stride,
                             const ogs_type type,
                             const ogs_op op,
                             const ogs_transpose trans,
                             ogs_t &ogs){

  const size_t Nbytes = ogs_type_size[type];

  void* gsh = (trans == ogs_sym) ? ogs.gshSym : ogs.gsh;

  if (trans == ogs_notrans) {
    if(ogs.NlocalFused)
      occaGatherScatterKernel(ogs.NlocalFused, Nentries, Nvectors, stride,
                              ogs.o_localFusedGatherOffsets,  ogs.o_localFusedGatherIds,
                              ogs.o_localFusedScatterOffsets, ogs.o_localFusedScatterIds,
                              type, op, o_v);
  } else if (trans == ogs_trans) {
    if(ogs.NlocalFused)
      occaGatherScatterKernel(ogs.NlocalFused, Nentries, Nvectors, stride,
                              ogs.o_localFusedScatterOffsets, ogs.o_localFusedScatterIds,
                              ogs.o_localFusedGatherOffsets,  ogs.o_localFusedGatherIds,
                              type, op, o_v);
  } else {//ogs_sym
    if(ogs.NlocalFusedSym)
      occaGatherScatterKernel(ogs.NlocalFusedSym, Nentries, Nvectors, stride,
                              ogs.o_localFusedOffsets, ogs.o_localFusedIds,
                              ogs.o_localFusedOffsets, ogs.o_localFusedIds,
                              type, op, o_v);
  }

  occa::stream currentStream = ogs.device.getStream();
  if (ogs.Nhalo) {
    ogs.device.setStream(dataStream);
    ogs.device.finish();
    ogs.device.setStream(currentStream);
  }

  // MPI based gather scatter using libgs
  gsGatherScatter(ogs.haloBuf, Nentries, Nvectors, ogs.Nhalo,
                  type, op, trans, gsh);

  dlong NhaloScatter = (trans == ogs_trans) ? ogs.NhaloGather : ogs.NhaloScatter;

  if (NhaloScatter) {
    ogs.device.setStream(dataStream);

    // copy gatherScattered halo data back from HOST to DEVICE
    for (int i=0;i<Nvectors;i++)
      ogs.o_haloBuf.copyFrom((char*)ogs.haloBuf + ogs.Nhalo*Nbytes*Nentries*i,
                             NhaloScatter*Nbytes*Nentries,
                             ogs.Nhalo*Nbytes*Nentries*i,
                             "async: true");

    ogs.device.finish();
    ogs.device.setStream(currentStream);

    // scatter back to local nodes
    if (trans == ogs_trans)
      occaScatterKernel(ogs.NhaloGather, Nentries, Nvectors, stride, ogs.Nhalo,
                        ogs.o_haloGatherOffsets, ogs.o_haloGatherIds,
                        type, op, ogs.o_haloBuf, o_v);
    else
      occaScatterKernel(ogs.NhaloScatter, Nentries, Nvectors, stride, ogs.Nhalo,
                        ogs.o_haloScatterOffsets, ogs.o_haloScatterIds,
                        type, op, ogs.o_haloBuf, o_v);
  }
}


#define SWITCH_TYPE_CASE(T) case ogs_##T: { WITH_TYPE(T); break; }
#define SWITCH_TYPE(type) switch(type) { \
    OGS_FOR_EACH_TYPE(SWITCH_TYPE_CASE) case ogs_type_n: break; }

#define SWITCH_OP_CASE(T,OP) case ogs_##OP: { WITH_OP(T,OP); break; }
#define SWITCH_OP(T,op) switch(op) { \
    OGS_FOR_EACH_OP(T,SWITCH_OP_CASE) case ogs_op_n: break; }


void occaGatherScatterKernel(const dlong N,
                             const int Nentries,
                             const int Nvectors,
                             const dlong stride,
                             occa::memory& o_gatherStarts,
                             occa::memory& o_gatherIds,
                             occa::memory& o_scatterStarts,
                             occa::memory& o_scatterIds,
                             const ogs_type type,
                             const ogs_op op,
                             occa::memory&  o_v) {

#define WITH_OP(T,OP)                             \
  gatherScatterKernel_##T##_##OP(N,               \
                                 Nentries,        \
                                 Nvectors,        \
                                 stride,          \
                                 o_gatherStarts,  \
                                 o_gatherIds,     \
                                 o_scatterStarts, \
                                 o_scatterIds,    \
                                 o_v);
#define WITH_TYPE(T) SWITCH_OP(T,op)

  SWITCH_TYPE(type)

#undef  WITH_TYPE
#undef  WITH_OP
}

} //namespace ogs