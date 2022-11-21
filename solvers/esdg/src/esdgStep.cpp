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
void esdg_t::rhsf(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

  // use this to determine cutoff parameter for entropy generation when using DODGES
  static int entropyStep=0;
#if 1
  if(entropyStep%100==0){
    
    dfloat totalEntropy = integrateEntropy(o_Q);
    
    printf("%d %d %d %15.14lg % 15.14lg  %%%% N, Ncub, Nelements, Time, Total Entropy, cutoff\n",
	   mesh.N, mesh.cubN, mesh.Nelements, T, totalEntropy);
  }
  ++entropyStep;

#endif

  // entropy stable
  esInterpolateKernel(mesh.Nelements,
		      gamma,
		      o_EToE,
		      mesh.o_EToB,
		      o_esIqfT,
		      o_esFqT,
		      mesh.o_vgeo,
		      o_esWq,
		      o_Q,
		      o_esQc,
		      o_esQe);

  // fix elements with entropy violations using Hadamard product
  esVolumeKernel(mesh.Nelements,
		 gamma,
		 mesh.o_vgeo,
		 o_esQNrT,
		 o_esQNsT,
		 o_esPqLfT,
		 o_esQe,
		 o_RHS);
  
  // surface terms
  esSurfaceKernel(mesh.Nelements,
		  T,
		  gamma,
		  mesh.o_sgeo,
		  o_esVmapM,
		  o_esVmapP,
		  o_EToB,
		  o_esLfT,
		  o_esX,
		  o_esY,
		  o_esZ,
		  o_esQc,
		  o_esQe,
		  o_esQp,
		  o_esQcrecon,
		  o_RHS);
  
}
