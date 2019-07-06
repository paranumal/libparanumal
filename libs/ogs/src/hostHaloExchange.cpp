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

void halo_t::Exchange(void *v, const int Nentries, const ogs_type type) {
  ExchangeStart (v, Nentries, type);
  ExchangeFinish(v, Nentries, type);
}

void halo_t::ExchangeStart(void *v, const int Nentries, const ogs_type type){

  const size_t Nbytes = ogs_type_size[type];

  reallocHostBuffer(Nbytes*Nentries);

  // copy data from outgoing elements into halo buffer
  if (Nsend)
    ogs::hostExtractKernel(Nsend, Nentries, haloSendIds, type, v, hostBuf);
}


void halo_t::ExchangeFinish(void *v, const int Nentries, const ogs_type type){

  const size_t Nbytes = ogs_type_size[type];

  // MPI based scatter using gslib
  ogs::gsScatter(hostBuf, Nentries, 1, 0, type, ogs_add, gshHalo);

  if (Nhalo)
    memcpy((char*)v + Nbytes*Nentries*Nlocal,
           (char*)hostBuf + Nbytes*Nentries*Nsend,
           Nbytes*Nentries*Nhalo);
}

namespace ogs {

/*------------------------------------------------------------------------------
  The basic halo extract kernel
------------------------------------------------------------------------------*/
#define DEFINE_EXTRACT(T)                                                       \
static void hostHaloExtractKernel_##T(const dlong Ngather,                      \
                                      const int   Nentries,                     \
                                      const dlong *haloSendIds,                 \
                                      const     T *v,                           \
                                                T *haloBuf)                     \
{                                                                               \
  for(dlong g=0;g<Ngather*Nentries;++g){                                        \
    const dlong gid = g/Nentries;                                               \
    const int k     = g%Nentries;                                               \
    const dlong id  = haloSendIds[gid];                                         \
    haloBuf[g] = v[k+id*Nentries];                                              \
  }                                                                             \
}

#define DEFINE_PROCS(T) \
  DEFINE_EXTRACT(T)

OGS_FOR_EACH_TYPE(DEFINE_PROCS)

#define SWITCH_TYPE_CASE(T) case ogs_##T: { WITH_TYPE(T); break; }
#define SWITCH_TYPE(type) switch(type) { \
    OGS_FOR_EACH_TYPE(SWITCH_TYPE_CASE) case ogs_type_n: break; }

void hostExtractKernel(const dlong Nsend,
                       const int Nentries,
                       dlong *haloSendIds,
                       const ogs_type type,
                       void *v,
                       void *hostBuf) {

#define WITH_TYPE(T)                        \
  hostHaloExtractKernel_##T(Nsend,          \
                            Nentries,       \
                            haloSendIds,    \
                            (T*) v,         \
                            (T*) hostBuf);

  SWITCH_TYPE(type)

#undef  WITH_TYPE
}

} //namespace ogs