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

#include "esdg.hpp"

//evaluate ODE rhs = f(q,t)
dfloat esdg_t::integrateEntropyChange(occa::memory& o_Q, occa::memory& o_RHS){

  dfloat totalEntropy = 0;
  
  
  // compute entropy
  esIntegrateEntropyChangeKernel(mesh.Nelements,
				 gamma,
				 o_esIqfT,
				 o_esFqT,
				 mesh.o_vgeo,
				 o_esWq,
				 o_Q,
				 o_RHS,
				 o_esTotalEntropy);
  
  o_esTotalEntropy.copyTo(esTotalEntropy);
  
  for(dlong e=0;e<mesh.Nelements;++e){
    totalEntropy += esTotalEntropy[e];
  }

  return totalEntropy;
}


dfloat esdg_t::integrateEntropy(occa::memory& o_Q){
  
  dfloat totalEntropy = 0;
  
  // compute entropy
  esIntegrateEntropyKernel(mesh.Nelements,
			   gamma,
			   o_esIqfT,
			   mesh.o_vgeo,
			   o_esWq,
			   o_Q,
			   o_esTotalEntropy);
  
  o_esTotalEntropy.copyTo(esTotalEntropy);
  
  for(dlong e=0;e<mesh.Nelements;++e){
    totalEntropy += esTotalEntropy[e];
  }
  
  return totalEntropy;
}
