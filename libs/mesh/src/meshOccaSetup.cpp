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

  // find elements that have all neighbors on this process
  internalElementIds = (dlong*) calloc(Nelements, sizeof(dlong));
  notInternalElementIds = (dlong*) calloc(Nelements, sizeof(dlong));

  dlong Ninterior = 0, NnotInterior = 0;
  for(dlong e=0;e<Nelements;++e){
    int flag = 0;
    for(int f=0;f<Nfaces;++f)
      if(EToP[e*Nfaces+f]!=-1)
        flag = 1;
    if(!flag)
      internalElementIds[Ninterior++] = e;
    else
      notInternalElementIds[NnotInterior++] = e;
  }

  //printf("NinteriorElements = %d, NnotInternalElements = %d\n", Ninterior, NnotInterior);

  NinternalElements = Ninterior;
  NnotInternalElements = NnotInterior;
  if(Ninterior)
    o_internalElementIds    = device.malloc(Ninterior*sizeof(dlong), internalElementIds);

  if(NnotInterior)
    o_notInternalElementIds = device.malloc(NnotInterior*sizeof(dlong), notInternalElementIds);

  o_vmapM = device.malloc(Nelements*Nfp*Nfaces*sizeof(dlong), vmapM);
  o_vmapP = device.malloc(Nelements*Nfp*Nfaces*sizeof(dlong), vmapP);

  o_EToB = device.malloc(Nelements*Nfaces*sizeof(int), EToB);

  if(totalHaloPairs){
    // copy halo element list to DEVICE
    o_haloElementList = device.malloc(totalHaloPairs*sizeof(dlong), haloElementList);

    o_haloGetNodeIds = device.malloc(Nfp*totalHaloPairs*sizeof(dlong), haloGetNodeIds);
    o_haloPutNodeIds = device.malloc(Nfp*totalHaloPairs*sizeof(dlong), haloPutNodeIds);
  }

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
