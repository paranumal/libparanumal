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

#include "mesh.hpp"
#include "core.hpp"

// send data from partition boundary elements
// and receive data to ghost elements
void mesh_t::HaloExchange(void  *v, const int Nentries, const char *type){
  HaloExchangeStart(v, Nentries, type);
  HaloExchangeFinish(v, Nentries, type);
}

void mesh_t::HaloExchangeStart(void  *v, const int Nentries, const char *type){

  size_t Nbytes=0;

  if (!strcmp(type, "double"))
    Nbytes = sizeof(double);
  else if (!strcmp(type, "float"))
    Nbytes = sizeof(float);
  else if (!strcmp(type, "int"))
    Nbytes = sizeof(int);
  else if (!strcmp(type, "long long int"))
    Nbytes = sizeof(long long int);

  if (haloBufferSize < Nbytes*Nentries*(NhaloElements+totalHaloPairs)) {
    if (haloBuffer) free(haloBuffer);
    haloBuffer = malloc(Nbytes*Nentries*(NhaloElements+totalHaloPairs));
    haloBufferSize = Nbytes*Nentries*(NhaloElements+totalHaloPairs);
  }

  // copy data from outgoing elements into halo buffer
  for(dlong i=0;i<NhaloElements;++i){
    // outgoing element
    dlong e = haloElementIds[i];
    memcpy(((char*)haloBuffer)+i*Nbytes*Nentries,
           ((char*)v)         +e*Nbytes*Nentries, Nbytes*Nentries);
  }
}

void mesh_t::HaloExchangeFinish(void  *v, const int Nentries, const char *type){

  size_t Nbytes=0;

  if (!strcmp(type, "double"))
    Nbytes = sizeof(double);
  else if (!strcmp(type, "float"))
    Nbytes = sizeof(float);
  else if (!strcmp(type, "int"))
    Nbytes = sizeof(int);
  else if (!strcmp(type, "long long int"))
    Nbytes = sizeof(long long int);

  ogsGatherScatterVec(haloBuffer, Nentries, type, ogsAdd, ogsHalo);

  // copy data from incoming elements into end of source array
  memcpy(((char*)v)          + Nelements*Nbytes*Nentries,
         ((char*)haloBuffer) + NhaloElements*Nbytes*Nentries,
         totalHaloPairs*Nbytes*Nentries);
}

void mesh_t::HaloExchange(occa::memory &o_v, const int Nentries, const char *type){
  HaloExchangeStart(o_v, Nentries, type);
  HaloExchangeFinish(o_v, Nentries, type);
}

void mesh_t::HaloExchangeStart(occa::memory &o_v, const int Nentries, const char *type){
  size_t Nbytes=0;

  if (!strcmp(type, "double"))
    Nbytes = sizeof(double);
  else if (!strcmp(type, "float"))
    Nbytes = sizeof(float);
  else if (!strcmp(type, "int"))
    Nbytes = sizeof(int);
  else if (!strcmp(type, "long long int"))
    Nbytes = sizeof(long long int);

  if (o_haloBuffer.size() < Nbytes*Nentries*(NhaloElements+totalHaloPairs)) {
    if (o_haloBuffer.size()) {
      o_haloBuffer.free();
      h_haloBuffer.free();
    }

    pinnedHaloBuffer = occaHostMallocPinned(device,
                            Nbytes*Nentries*(NhaloElements+totalHaloPairs),
                            NULL, o_haloBuffer, h_haloBuffer);
  }

  if (NhaloElements) {
    // copy data from outgoing elements into halo buffer
    //WARNING: this kernel is currently only supporting dfloat types
    haloExtractKernel(NhaloElements, Nentries, o_haloElementIds, o_v, o_haloBuffer);

    device.finish();
    device.setStream(dataStream);
    o_haloBuffer.copyTo(pinnedHaloBuffer,
                        Nbytes*Nentries*NhaloElements, 0, "async: true");
    device.setStream(defaultStream);
  }
}

void mesh_t::HaloExchangeFinish(occa::memory &o_v, const int Nentries, const char *type){
  size_t Nbytes=0;

  if (!strcmp(type, "double"))
    Nbytes = sizeof(double);
  else if (!strcmp(type, "float"))
    Nbytes = sizeof(float);
  else if (!strcmp(type, "int"))
    Nbytes = sizeof(int);
  else if (!strcmp(type, "long long int"))
    Nbytes = sizeof(long long int);

  if (NhaloElements) {
    device.setStream(dataStream);
    device.finish();
  }

  ogsGatherScatterVec(pinnedHaloBuffer, Nentries, type, ogsAdd, ogsHalo);

  if (NhaloElements) {
    // copy data from incoming elements into end of source array
    o_v.copyFrom(((char*)pinnedHaloBuffer) + Nbytes*Nentries*NhaloElements,
                 Nbytes*Nentries*totalHaloPairs,
                 Nbytes*Nentries*Nelements, "async: true");
    device.finish();
    device.setStream(defaultStream);
  }
}
