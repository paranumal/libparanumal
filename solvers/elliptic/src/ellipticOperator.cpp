/*

  The MIT License (MIT)

  Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "elliptic.hpp"

void elliptic_t::Operator(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_Aq){

  if(disc_c0){
    // int mapType = (mesh.elementType==Mesh::HEXAHEDRA &&
    //                mesh.settings.compareSetting("ELEMENT MAP", "TRILINEAR")) ? 1:0;

    // int integrationType = (mesh.elementType==Mesh::HEXAHEDRA &&
    //                        settings.compareSetting("ELLIPTIC INTEGRATION", "CUBATURE")) ? 1:0;

    gHalo.ExchangeStart(o_q, 1);

    if(mesh.NlocalGatherElements/2){
      // if(integrationType==0) { // GLL or non-hex
        // if(mapType==0)
          partialAxKernel(mesh.NlocalGatherElements/2,
                          mesh.o_localGatherElementList,
                          o_GlobalToLocal,
                          mesh.o_wJ, mesh.o_ggeo,
                          mesh.o_D, mesh.o_S,
                          mesh.o_MM, lambda, o_q, o_AqL);
        /* NC: disabling until we re-add treatment of affine elements
        else
          partialAxKernel(mesh.NlocalGatherElements, mesh.o_localGatherElementList,
                          mesh.o_EXYZ, mesh.o_gllzw, mesh.o_D, mesh.o_S, mesh.o_MM, lambda, o_q, o_Aq);
        */
      // } else {
      //   partialCubatureAxKernel(mesh.NlocalGatherElements,
      //                           mesh.o_localGatherElementList,
      //                           mesh.o_cubggeo,
      //                           mesh.o_cubD,
      //                           mesh.o_cubInterpT,
      //                           lambda,
      //                           o_q,
      //                           o_Aq);
      // }
    }

    // finalize halo exchange
    gHalo.ExchangeFinish(o_q, 1);

    if(mesh.NglobalGatherElements) {

      // if(integrationType==0) { // GLL or non-hex
        // if(mapType==0)
          partialAxKernel(mesh.NglobalGatherElements,
                          mesh.o_globalGatherElementList,
                          o_GlobalToLocal,
                          mesh.o_wJ, mesh.o_ggeo,
                          mesh.o_D, mesh.o_S,
                          mesh.o_MM, lambda, o_q, o_AqL);
        /* NC: disabling until we re-add treatment of affine elements
        else
          partialAxKernel(mesh.NglobalGatherElements, mesh.o_globalGatherElementList,
                          mesh.o_EXYZ, mesh.o_gllzw, mesh.o_D, mesh.o_S, mesh.o_MM, lambda, o_q, o_Aq);
        */
      // } else {
      //   partialCubatureAxKernel(mesh.NglobalGatherElements,
      //                           mesh.o_globalGatherElementList,
      //                           mesh.o_cubggeo,
      //                           mesh.o_cubD,
      //                           mesh.o_cubInterpT,
      //                           lambda, o_q, o_Aq);
      // }
    }

    //gather result to Aq
    ogsMasked.GatherStart(o_Aq, o_AqL, 1, ogs::Add, ogs::Trans);

    if((mesh.NlocalGatherElements+1)/2){
      partialAxKernel((mesh.NlocalGatherElements+1)/2,
                      mesh.o_localGatherElementList+(mesh.NlocalGatherElements/2),
                      o_GlobalToLocal,
                      mesh.o_wJ, mesh.o_ggeo,
                      mesh.o_D, mesh.o_S,
                      mesh.o_MM, lambda, o_q, o_AqL);
    }

    ogsMasked.GatherFinish(o_Aq, o_AqL, 1, ogs::Add, ogs::Trans);

  } else if(disc_ipdg) {

    if(mesh.Nelements) {
      dlong offset = 0;
      partialGradientKernel(mesh.Nelements,
                            offset,
                            mesh.o_vgeo,
                            mesh.o_D,
                            o_q,
                            o_grad);
    }

    // dfloat4 storage -> 4 entries
    traceHalo.ExchangeStart(o_grad, 4);

    if(mesh.NinternalElements)
      partialIpdgKernel(mesh.NinternalElements,
                        mesh.o_internalElementIds,
                        mesh.o_vmapM,
                        mesh.o_vmapP,
                        lambda,
                        tau,
                        mesh.o_vgeo,
                        mesh.o_sgeo,
                        o_EToB,
                        mesh.o_D,
                        mesh.o_LIFT,
                        mesh.o_MM,
                        o_grad,
                        o_Aq);

    traceHalo.ExchangeFinish(o_grad, 4);

    if(mesh.NhaloElements) {
      partialIpdgKernel(mesh.NhaloElements,
                        mesh.o_haloElementIds,
                        mesh.o_vmapM,
                        mesh.o_vmapP,
                        lambda,
                        tau,
                        mesh.o_vgeo,
                        mesh.o_sgeo,
                        o_EToB,
                        mesh.o_D,
                        mesh.o_LIFT,
                        mesh.o_MM,
                        o_grad,
                        o_Aq);
    }
  }
}

