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

OGS_DEFINE_TYPE_SIZES()

void halo_t::Exchange(occa::memory& o_v, const int Nentries, const ogs_type type) {
  ExchangeStart (o_v, Nentries, type);
  ExchangeFinish(o_v, Nentries, type);
}

void halo_t::ExchangeStart(occa::memory& o_v,
                           const int Nentries,
                           const ogs_type type){

  const size_t Nbytes = ogs_type_size[type];

  reallocOccaBuffer(Nbytes*Nentries);

  if (Nsend) {
    // copy data from outgoing elements into halo buffer
    ogs::occaExtractKernel(Nsend, Nentries, o_haloSendIds, type, o_v, o_haloBuf);

    device.finish();
    occa::stream currentStream = device.getStream();
    device.setStream(ogs::dataStream);
    o_haloBuf.copyTo(haloBuf,
                     Nbytes*Nentries*Nsend, 0, "async: true");
    device.setStream(currentStream);
  }
}

void halo_t::ExchangeFinish(occa::memory& o_v,
                            const int Nentries,
                            const ogs_type type){

  const size_t Nbytes = ogs_type_size[type];

  occa::stream currentStream = device.getStream();
  if (Nsend) {
    device.setStream(ogs::dataStream);
    device.finish();
    device.setStream(currentStream);
  }

  // MPI based scatter using gslib
  ogs::gsScatter(haloBuf, Nentries, 1, 0, type, ogs_add, gshHalo);

  if (Nhalo) {
    device.setStream(ogs::dataStream);
    // copy data from incoming elements into end of source array
    o_v.copyFrom(((char*)haloBuf) + Nbytes*Nentries*Nsend,
                 Nbytes*Nentries*Nhalo,
                 Nbytes*Nentries*Nlocal, "async: true");
    device.finish();
    device.setStream(currentStream);
  }
}

namespace ogs {

#define SWITCH_TYPE_CASE(T) case ogs_##T: { WITH_TYPE(T); break; }
#define SWITCH_TYPE(type) switch(type) { \
    OGS_FOR_EACH_TYPE(SWITCH_TYPE_CASE) case ogs_type_n: break; }

void occaExtractKernel(const dlong Nsend,
                       const int Nentries,
                       occa::memory& o_haloSendIds,
                       const ogs_type type,
                       occa::memory& o_v,
                       occa::memory& o_haloBuf) {

#define WITH_TYPE(T)                    \
  haloExtractKernel_##T(Nsend,          \
                        Nentries,       \
                        o_haloSendIds,  \
                        o_v,            \
                        o_haloBuf);

  SWITCH_TYPE(type)

#undef  WITH_TYPE
}

} //namespace ogs