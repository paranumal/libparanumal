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

void occaGatherStart(occa::memory& o_gv,
                     occa::memory& o_v,
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


void occaGatherFinish(occa::memory& o_gv,
                      occa::memory& o_v,
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

  if(ogs.Nlocal) {
    if (trans == ogs_notrans)
      occaGatherKernel(ogs.Nlocal, Nentries, Nvectors, stride, gstride,
                       ogs.o_localGatherOffsets, ogs.o_localGatherIds,
                       type, op, o_v, o_gv);
    else
      occaGatherKernel(ogs.Nlocal, Nentries, Nvectors, stride, gstride,
                       ogs.o_localScatterOffsets, ogs.o_localScatterIds,
                       type, op, o_v, o_gv);
  }

  occa::stream currentStream = ogs.device.getStream();
  if (ogs.Nhalo) {
    ogs.device.setStream(dataStream);
    ogs.device.finish();
    ogs.device.setStream(currentStream);
  }

  // MPI based gather using libgs
  gsGatherScatter(ogs.haloBuf, Nentries, Nvectors, ogs.Nhalo,
                  type, op, trans, ogs.gsh);

  // copy totally gathered halo data back from HOST to DEVICE
  if (ogs.NhaloGather) {
    ogs.device.setStream(dataStream);

    for (int i=0;i<Nvectors;i++)
      o_gv.copyFrom((char*)ogs.haloBuf+ogs.Nhalo*Nbytes*Nentries*i,
                    ogs.NhaloGather*Nbytes*Nentries,
                    ogs.NlocalGather*Nbytes + gstride*Nbytes*i,
                    "async: true");

    ogs.device.finish();
    ogs.device.setStream(currentStream);
  }
}


#define SWITCH_TYPE_CASE(T) case ogs_##T: { WITH_TYPE(T); break; }
#define SWITCH_TYPE(type) switch(type) { \
    OGS_FOR_EACH_TYPE(SWITCH_TYPE_CASE) case ogs_type_n: break; }

#define SWITCH_OP_CASE(T,OP) case ogs_##OP: { WITH_OP(T,OP); break; }
#define SWITCH_OP(T,op) switch(op) { \
    OGS_FOR_EACH_OP(T,SWITCH_OP_CASE) case ogs_op_n: break; }


void occaGatherKernel(const dlong N,
                      const int Nentries,
                      const int Nvectors,
                      const dlong stride,
                      const dlong gstride,
                      occa::memory& o_gatherStarts,
                      occa::memory& o_gatherIds,
                      const ogs_type type,
                      const ogs_op op,
                      occa::memory& o_v,
                      occa::memory& o_gv) {

#define WITH_OP(T,OP)                     \
  gatherKernel_##T##_##OP(N,              \
                          Nentries,       \
                          Nvectors,       \
                          stride,         \
                          gstride,        \
                          o_gatherStarts, \
                          o_gatherIds,    \
                          o_v,            \
                          o_gv);
#define WITH_TYPE(T) SWITCH_OP(T,op)

  SWITCH_TYPE(type)

#undef  WITH_TYPE
#undef  WITH_OP
}

} //namespace ogs