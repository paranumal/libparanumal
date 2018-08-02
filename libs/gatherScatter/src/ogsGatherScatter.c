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

#include "ogs.h"
#include "ogsKernels.h"
#include "ogsInterface.h"

void ogsGatherScatter(void *v, 
                      const char *type, 
                      const char *op, 
                      ogs2_t *ogs){
  ogsHostGatherScatter(v, type, op, ogs->hostGsh);
}

void ogsGatherScatter(occa::memory &o_v, 
                      const char *type, 
                      const char *op, 
                      ogs2_t *ogs){
  ogsGatherScatterStart (o_v, type, op, ogs);
  ogsGatherScatterFinish(o_v, type, op, ogs);
}

void ogsGatherScatterStart(occa::memory &o_v, 
                          const char *type, 
                          const char *op, 
                          ogs2_t *ogs){
  const int one = 1;
  const dlong dOne = 1;

  size_t Nbytes;
  if (!strcmp(type, "float")) 
    Nbytes = sizeof(float);
  else if (!strcmp(type, "double")) 
    Nbytes = sizeof(double);
  else if (!strcmp(type, "int")) 
    Nbytes = sizeof(int);
  else if (!strcmp(type, "long long int")) 
    Nbytes = sizeof(long long int);

  if (ogs->NhaloGather) {
    if (ogs::o_haloBuf.size() < ogs->NhaloGather*Nbytes) {
      if (ogs::o_haloBuf.size()) ogs::o_haloBuf.free();
      ogs::o_haloBuf = ogs->device.mappedAlloc(ogs->NhaloGather*Nbytes);
      ogs::haloBuf = ogs::o_haloBuf.getMappedPointer();
    }
  }

  // gather halo nodes on device
  if (ogs->NhaloGather) {
    ogs::gatherKernel(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, one, dOne, o_v, ogs::o_haloBuf);
    
    ogs->device.finish();
    ogs->device.setStream(ogs::dataStream);
    ogs::o_haloBuf.copyTo(ogs::haloBuf,"async: true");
    ogs->device.setStream(ogs::defaultStream);
  }
}


void ogsGatherScatterFinish(occa::memory &o_v, 
                          const char *type, 
                          const char *op, 
                          ogs2_t *ogs){
  const int one = 1;
  const dlong dOne = 1;

  if(ogs->NlocalGather) {
    ogs::gatherScatterKernel(ogs->NlocalGather, ogs->o_localGatherOffsets, ogs->o_localGatherIds, one, dOne, o_v);
  }

  if (ogs->NhaloGather) {
    ogs->device.setStream(ogs::dataStream);
    ogs->device.finish();

    // MPI based gather scatter using libgs
    ogsHostGatherScatter(ogs::haloBuf, type, op, ogs->haloGshSym);

    // copy totally gather halo data back from HOST to DEVICE
    ogs::o_haloBuf.copyFrom(ogs::haloBuf, "async: true");

    // do scatter back to local nodes
    ogs::scatterKernel(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, one, dOne, ogs::o_haloBuf, o_v);
    ogs->device.finish();
    ogs->device.setStream(ogs::defaultStream);
  }
}
