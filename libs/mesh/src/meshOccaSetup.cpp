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

void mesh_t::OccaSetup(occa::properties &kernelInfo){

  if(NinternalElements)
    o_internalElementIds    =
      device.malloc(NinternalElements*sizeof(dlong), internalElementIds);

  if(NhaloElements)
    o_haloElementIds = device.malloc(NhaloElements*sizeof(dlong), haloElementIds);

  if(NglobalGatherElements)
    o_globalGatherElementList =
      device.malloc(NglobalGatherElements*sizeof(dlong), globalGatherElementList);

  if(NlocalGatherElements)
    o_localGatherElementList =
      device.malloc(NlocalGatherElements*sizeof(dlong), localGatherElementList);

  o_vmapM = device.malloc(Nelements*Nfp*Nfaces*sizeof(dlong), vmapM);
  o_vmapP = device.malloc(Nelements*Nfp*Nfaces*sizeof(dlong), vmapP);

  o_EToB = device.malloc(Nelements*Nfaces*sizeof(int), EToB);

  defaultStream = device.getStream();
  dataStream    = device.createStream();

  haloExtractKernel = device.buildKernel(LIBP_DIR "/libs/mesh/okl/meshHaloExtract.okl",
                                         "meshHaloExtract",
                                         kernelInfo);

  kernelInfo["defines/" "p_dim"]= dim;
  kernelInfo["defines/" "p_N"]= N;
  kernelInfo["defines/" "p_Nq"]= N+1;
  kernelInfo["defines/" "p_Np"]= Np;
  kernelInfo["defines/" "p_Nfp"]= Nfp;
  kernelInfo["defines/" "p_Nfaces"]= Nfaces;
  kernelInfo["defines/" "p_NfacesNfp"]= Nfp*Nfaces;
  kernelInfo["defines/" "p_Nvgeo"]= Nvgeo;
  kernelInfo["defines/" "p_Nsgeo"]= Nsgeo;
  kernelInfo["defines/" "p_Nggeo"]= Nggeo;

  kernelInfo["defines/" "p_max_EL_nnz"]= max_EL_nnz; // for Bernstein Bezier lift

  kernelInfo["defines/" "p_cubNq"]= cubNq;
  kernelInfo["defines/" "p_cubNp"]= cubNp;
  kernelInfo["defines/" "p_intNfp"]= intNfp;
  kernelInfo["defines/" "p_intNfpNfaces"]= intNfp*Nfaces;
  kernelInfo["defines/" "p_cubNfp"]= cubNfp;
}
