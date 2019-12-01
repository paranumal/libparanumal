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

void mesh_t::OccaSetup(){

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
  o_mapP  = device.malloc(Nelements*Nfp*Nfaces*sizeof(dlong), mapP);

  o_EToB = device.malloc(Nelements*Nfaces*sizeof(int), EToB);

  defaultStream = device.getStream();

  props["defines/" "p_dim"]= dim;
  props["defines/" "p_N"]= N;
  props["defines/" "p_Nq"]= N+1;
  props["defines/" "p_Np"]= Np;
  props["defines/" "p_Nfp"]= Nfp;
  props["defines/" "p_Nfaces"]= Nfaces;
  props["defines/" "p_NfacesNfp"]= Nfp*Nfaces;
  props["defines/" "p_Nvgeo"]= Nvgeo;
  props["defines/" "p_Nsgeo"]= Nsgeo;
  props["defines/" "p_Nggeo"]= Nggeo;

  props["defines/" "p_cubNq"]= cubNq;
  props["defines/" "p_cubNp"]= cubNp;
  props["defines/" "p_intNfp"]= intNfp;
  props["defines/" "p_intNfpNfaces"]= intNfp*Nfaces;
  props["defines/" "p_cubNfp"]= cubNfp;
}
