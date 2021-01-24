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

#include "elliptic.hpp"

void elliptic_t::Operator(occa::memory &o_q, occa::memory &o_Aq){

  if(disc_c0){

    int useCubature = (mesh.elementType==HEXAHEDRA || mesh.elementType==QUADRILATERALS) &&
      settings.compareSetting("ELLIPTIC INTEGRATION", "CUBATURE");
  
    ogsMasked->GatheredHaloExchangeStart(o_q, 1, ogs_dfloat);

    int useGlobalToLocal = 1;
    
    if(mesh.NlocalGatherElements){
      if(useCubature==0)
	partialAxKernel(mesh.NlocalGatherElements,
			mesh.o_localGatherElementList,
			ogsMasked->o_GlobalToLocal,
			useGlobalToLocal,
			mesh.o_ggeo, mesh.o_D, mesh.o_S,
			mesh.o_MM, lambda, o_q, o_AqL);
      else
	partialCubatureAxKernel(mesh.NlocalGatherElements,
				mesh.o_localGatherElementList,
				ogsMasked->o_GlobalToLocal,
				useGlobalToLocal,
				mesh.o_cubggeo,
				mesh.o_cubD,
				mesh.o_cubInterp,
				lambda,
				o_q,
				o_Aq);
      
    }

    // finalize halo exchange
    ogsMasked->GatheredHaloExchangeFinish(o_q, 1, ogs_dfloat);

    if(mesh.NglobalGatherElements) {
      if(useCubature==0)
	partialAxKernel(mesh.NglobalGatherElements,
			mesh.o_globalGatherElementList,
			ogsMasked->o_GlobalToLocal,
			useGlobalToLocal,
			mesh.o_ggeo, mesh.o_D, mesh.o_S,
			mesh.o_MM, lambda, o_q, o_AqL);
      else
	partialCubatureAxKernel(mesh.NglobalGatherElements,
				mesh.o_globalGatherElementList,
				ogsMasked->o_GlobalToLocal,
				useGlobalToLocal,
				mesh.o_cubggeo,
				mesh.o_cubD,
				mesh.o_cubInterp,
				lambda, o_q, o_Aq);
      
    }

    //gather result to Aq
    ogsMasked->Gather(o_Aq, o_AqL, ogs_dfloat, ogs_add, ogs_trans);

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
    traceHalo->ExchangeStart(o_grad, 4, ogs_dfloat);

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

    traceHalo->ExchangeFinish(o_grad, 4, ogs_dfloat);

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

