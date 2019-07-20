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

void occaScatterStart(occa::memory& o_v,
                      occa::memory& o_gv,
                      const int Nentries,
                      const int Nvectors,
                      const dlong stride,
                      const dlong gstride,
                      const ogs_type type,
                      const ogs_op op,
                      const ogs_transpose trans,
                      ogs_t &ogs){

  const size_t Nbytes = ogs_type_size[type];

  if (trans == ogs_sym)
    LIBP_ABORT(string("Calling ogs::Scatter in ogs_sym mode not supported."))

  ogs.reallocOccaBuffer(Nbytes*Nentries*Nvectors);

  if (ogs.NhaloGather) {
    occa::stream currentStream = ogs.device.getStream();
    ogs.device.setStream(dataStream);

    for (int i=0;i<Nvectors;i++)
      o_gv.copyTo((char*)ogs.haloBuf+ogs.Nhalo*Nbytes*Nentries*i,
                  ogs.NhaloGather*Nbytes*Nentries,
                  ogs.NlocalGather*Nbytes + gstride*Nbytes*i,
                  "async: true");

    ogs.device.setStream(currentStream);
  }
}


void occaScatterFinish(occa::memory& o_v,
                       occa::memory& o_gv,
                       const int Nentries,
                       const int Nvectors,
                       const dlong stride,
                       const dlong gstride,
                       const ogs_type type,
                       const ogs_op op,
                       const ogs_transpose trans,
                       ogs_t &ogs){

  const size_t Nbytes = ogs_type_size[type];

  if (trans == ogs_sym)
    LIBP_ABORT(string("Calling ogs::Scatter in ogs_sym mode not supported."))

  if(ogs.Nlocal) {
    if (trans == ogs_notrans)
      occaScatterKernel(ogs.Nlocal, Nentries, Nvectors, gstride, stride,
                        ogs.o_localScatterOffsets, ogs.o_localScatterIds,
                        type, op, o_gv, o_v);
    else
      occaScatterKernel(ogs.Nlocal, Nentries, Nvectors, gstride, stride,
                        ogs.o_localGatherOffsets, ogs.o_localGatherIds,
                        type, op, o_gv, o_v);
  }

  occa::stream currentStream = ogs.device.getStream();
  if (ogs.Nhalo) {
    ogs.device.setStream(dataStream);
    ogs.device.finish();
    ogs.device.setStream(currentStream);
  }

  // MPI based scatter using gslib
  // (must use ogs_notrans so the negative ids don't contribute to op)
  gsGatherScatter(ogs.haloBuf, Nentries, Nvectors, ogs.Nhalo,
                  type, op, ogs_notrans, ogs.gsh);

  dlong NhaloScatter = (trans == ogs_notrans) ? ogs.NhaloScatter : ogs.NhaloGather;

  if (NhaloScatter) {
    ogs.device.setStream(dataStream);

    // copy totally scattered halo data back from HOST to DEVICE
    for (int i=0;i<Nvectors;i++)
      ogs.o_haloBuf.copyFrom((char*)ogs.haloBuf + ogs.Nhalo*Nbytes*Nentries*i,
                             NhaloScatter*Nbytes*Nentries*Nvectors,
                             ogs.Nhalo*Nbytes*Nentries*i,
                             "async: true");

    ogs.device.finish();
    ogs.device.setStream(currentStream);

    if (trans == ogs_notrans)
      occaScatterKernel(ogs.NhaloScatter, Nentries, Nvectors, ogs.Nhalo, stride,
                        ogs.o_haloScatterOffsets, ogs.o_haloScatterIds,
                        type, op, ogs.o_haloBuf, o_v);
    else
      occaScatterKernel(ogs.NhaloGather, Nentries, Nvectors, ogs.Nhalo, stride,
                        ogs.o_haloGatherOffsets, ogs.o_haloGatherIds,
                        type, op, ogs.o_haloBuf, o_v);
  }
}


#define SWITCH_TYPE_CASE(T) case ogs_##T: { WITH_TYPE(T); break; }
#define SWITCH_TYPE(type) switch(type) { \
    OGS_FOR_EACH_TYPE(SWITCH_TYPE_CASE) case ogs_type_n: break; }

void occaScatterKernel(const dlong N,
                       const int Nentries,
                       const int Nvectors,
                       const dlong gstride,
                       const dlong stride,
                       occa::memory& o_scatterStarts,
                       occa::memory& o_scatterIds,
                       const ogs_type type,
                       const ogs_op op,
                       occa::memory& o_gv,
                       occa::memory& o_v) {

#define WITH_TYPE(T)                 \
  scatterKernel_##T(N,               \
                    Nentries,        \
                    Nvectors,        \
                    gstride,         \
                    stride,          \
                    o_scatterStarts, \
                    o_scatterIds,    \
                    o_gv,            \
                    o_v);

  SWITCH_TYPE(type)

#undef  WITH_TYPE
}

} //namespace ogs