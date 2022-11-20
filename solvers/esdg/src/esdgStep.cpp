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
void esdg_t::rhsf(occa::memory& o_Q, occa::memory& o_RHS, const dfloat T){

  // modifications of based schemes
  dfloat lambda = 0, mu = 0;
  dfloat tau = 0;

  settings.getSetting("LAME LAMBDA", lambda);
  settings.getSetting("LAME MU", mu);
  settings.getSetting("DIFFUSION PENALTY TAU", tau);

  // use this to determine cutoff parameter for entropy generation when using DODGES
  static int entropyStep=0;
#if 0
  if(timeStepper.rkStage==0){
    if(entropyStep%100==0){

      dfloat totalEntropy = integrateEntropy(o_Q);
      
      printf("%d %d %d %15.14lg % 15.14lg  %%%% N, Ncub, Nelements, Time, Total Entropy, cutoff\n",
	     mesh.N, mesh.cubN, mesh.Nelements, T, totalEntropy);
    }
    ++entropyStep;
  }
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
  
  // try cubature volume term
  esVolumeCubatureKernel(mesh.Nelements,
			 gamma,
			 mesh.o_vgeo,
#if 0
			 // FIX THIS LATER
			 mesh.o_cubDrWT, 
			 mesh.o_cubDsWT,
#endif
			 o_esPqLfT,
			 o_esItMT,
			 o_esQc,
			 o_esQe,
			 o_RHS,
			 o_entropyChange);
  
  esSurfaceCubatureKernel(mesh.Nelements,
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
			  o_RHS,
			  o_esWf,
			  o_entropyChange);


#if 0
  // DO DIFFUSION LATER
    // add diffusion terms
  // ADD CNS diffusion using entropy  variables
  esVolumeGradientKernel(mesh.Nelements,
			 mesh.o_vgeo,
			 o_esIqfDrPqT,
			 o_esIqfDsPqT, 
			 o_esQe,
			 o_esdQedx,
			 o_esdQedy,
			 o_esdQedz);
  
  esSurfaceGradientKernel(mesh.Nelements,
			  T,
			  mesh.o_sgeo, // affine ftm
			  o_esX,
			  o_esY,
			  o_esZ,
			  o_esVmapM,
			  o_esVmapP,
			  o_EToB,
			  o_esIqfLfT,
			  o_esQe,
			  o_esdQedx,
			  o_esdQedy,
			  o_esdQedz);
      
  esDiffusionFluxesKernel(mesh.Nelements,
			  o_esR,
			  o_esS,
			  o_esRq,
			  o_esSq,
			  o_esMu, // mu
			  mesh.o_vgeo,
			  o_esFqT,
			  o_esQe,
			  o_esdQedx,
			  o_esdQedy,
			  o_esdQedz);
  
  esVolumeDivergenceKernel(mesh.Nelements,
			   mesh.o_vgeo, // affine ftm
			   mesh.o_cubDrWT, 
			   mesh.o_cubDsWT, 
			   o_esQe,
			   o_esdQedx,
			   o_esdQedy,
			   o_esdQedz,
			   o_RHS);
  
  esSurfaceDivergenceKernel(mesh.Nelements,
			    T,
			    tau,
			    mesh.o_sgeo, // affine ftm
			    o_esX,
			    o_esY,
			    o_esZ,
			    o_esVmapM,
			    o_esVmapP,
			    o_EToB,
			    o_esLfT,
			    o_esQc,
			    o_esQe,
			    o_esdQedx,
			    o_esdQedy,
			    o_esdQedz,
			    o_RHS);

#endif
  
#if 0
  // do this in report
  if(!(step%1000)){
    dfloat totalEntropy = integrateEntropyChange(o_Q, o_RHS);
    if(!(step%1000))
      printf("T = %g, Total Entropy Produced in step = %17.15lg\n",
	     T, totalEntropy);
  }

  ++step;
  
  static int stepCounter = 0;
  
  if(!(stepCounter%1000)){
    
    printf("Step: %d, Time: %lg \n",
	   stepCounter, T, );
  }
  
  ++stepCounter;
#endif
}
