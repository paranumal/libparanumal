/*

The MIT License (MIT)

Copyright (c) 2020 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

OGS_DEFINE_TYPE_SIZES()

using namespace ogs;

void ogs_t::GatheredHaloExchangeStart(occa::memory& o_v,
                                  const int k,
                                  const ogs_type type){

  occa::device &device = platform.device;
  const size_t Nbytes = ogs_type_size[type];

  reallocOccaBuffer(Nbytes*k);

  if (haloGather.Nrows) {
    occa::stream currentStream = device.getStream();
    device.finish(); //make sure data is ready to copy
    device.setStream(dataStream);

    o_v.copyTo((char*)haloBuf,
                haloGather.Nrows*Nbytes*k,
                localGather.Nrows*Nbytes*k,
                "async: true");

    device.setStream(currentStream);
  }
}


void ogs_t::GatheredHaloExchangeFinish(occa::memory& o_v,
                                     const int k,
                                     const ogs_type type){

  occa::device &device = platform.device;
  const size_t Nbytes = ogs_type_size[type];

  occa::stream currentStream = device.getStream();
  if (Nhalo) {
    device.setStream(dataStream);
    device.finish();
    device.setStream(currentStream);
  }

  // MPI based scatter using gslib
  // (must use ogs_notrans so the negative ids don't contribute to op)
  gsGatherScatter(haloBuf, k, 1, Nhalo,
                  type, ogs_add, ogs_notrans, gsh);

  if (haloScatter.Nrows) {
    device.setStream(dataStream);

    // copy totally scattered halo data back from HOST to DEVICE
    o_v.copyFrom((char*)haloBuf+haloGather.Nrows*Nbytes*k,
                 (Nhalo-haloGather.Nrows)*Nbytes*k,
                 Ngather*Nbytes*k,
                 "async: true");

    device.finish();
    device.setStream(currentStream);
  }
}

/* Build global to local mapping */
void ogs_t::GatheredHaloExchangeSetup(){
  dlong *ids = (dlong*) malloc((Ngather+NgatherHalo)*sizeof(dlong));

  for (dlong n=0;n<Ngather+NgatherHalo;n++)
    ids[n] = n;

  GlobalToLocal = (dlong*) malloc(N*sizeof(dlong));

  for (dlong n=0;n<N;n++)
    GlobalToLocal[n] = -1;

  for (dlong i=0;i<localScatter.Nrows;i++) {
    const dlong start = localScatter.rowStarts[i];
    const dlong end   = localScatter.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong colId = localScatter.colIds[j];
      GlobalToLocal[colId] = ids[i];
    }
  }
  for (dlong i=0;i<haloScatter.Nrows;i++) {
    const dlong start = haloScatter.rowStarts[i];
    const dlong end   = haloScatter.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong colId = haloScatter.colIds[j];
      GlobalToLocal[colId] = ids[i+localScatter.Nrows];
    }
  }

  free(ids);

  o_GlobalToLocal = platform.malloc(N*sizeof(dlong), GlobalToLocal);
}