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

#include "esdg.hpp"

dfloat esdg_t::MaxWaveSpeed(deviceMemory<dfloat>& o_Q, const dfloat T){

  deviceMemory<dfloat> o_maxSpeed = platform.reserve<dfloat>(mesh.Nelements);

  std::cout << "WHYARE WE IN MAXWAVESPEED" << std::endl;
  
  maxWaveSpeedKernel(mesh.Nelements,
                     mesh.o_vgeo,
                     mesh.o_sgeo,
                     mesh.o_vmapM,
                     mesh.o_EToB,
                     gamma,
                     mu,
                     T,
                     mesh.o_x,
                     mesh.o_y,
                     mesh.o_z,
                     o_Q,
                     o_maxSpeed);

  const dfloat vmax = platform.linAlg().max(mesh.Nelements, o_maxSpeed, mesh.comm);

  return vmax;
}

//evaluate ODE rhs = f(q,t)
void esdg_t::rhsf(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
  
  int useArtificialViscosity =1;
  
  dfloat lambda = 0, mu = 0;
  dfloat tau = 0;
  settings.getSetting("LAME LAMBDA", lambda);
  settings.getSetting("LAME MU", mu);
  settings.getSetting("DIFFUSION PENALTY TAU", tau);

#if 0
  dfloat dt = timeStepper->dt;

  dfloat dtPrevious = timeStepper->dtPrevious;
  
  dfloat cutoff = (timeStepper->outputStep) ? pow(dtPrevious,3): pow(dt,3); // experimental
#endif
  
#if 0
  // ESDG
  static int entropyStep=0;
  if(timeStepper->rkStage==0){
    if(entropyStep%100==0){
      
      dfloat totalEntropy = integrateEntropy(o_Q);
      
      printf("%d %d %d %15.14lg % 15.14lg % 15.14lg %%%% N, Ncub, Nelements, Time, Total Entropy, cutoff\n",
	     mesh.N, mesh.cubN, mesh.Nelements, T, totalEntropy, cutoff);
    }
    ++entropyStep;
  }
#endif

  //  if(timeStepper.stage==0) // TW: warning need to initialize
  //    o_zeroFlag.copyTo(o_projFlag);
  
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
		      o_esQe,
		      o_projectionError);

  // fix elements with entropy violations using Hadamard product
  esVolumeKernel(mesh.Nelements,
		 gamma,
		 maxLVToE,
		 o_NLVToE,
		 o_LVToE,
		 o_esRq,
		 o_esSq,
		 mesh.o_vgeo,
		 o_esQNrT,
		 o_esQNsT,
		 o_esPqLfT,
		 o_esQe,
		 o_RHS,
		 o_projFlag,
		 o_projectionError);
  
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
		  o_RHS,
		  o_projFlag);

  // add diffusion terms
  
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
			  useArtificialViscosity,
			  maxLVToE,
			  o_NLVToE,
			  o_LVToE,
			  o_esR,
			  o_esS,
			  o_esRq,
			  o_esSq,
			  o_artificialViscosity,
			  o_esMu, // mu
			  mesh.o_vgeo,
			  o_esFqT,
			  o_projectionError,
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
			    o_projectionError,
			    o_RHS);

  static int stepCounter = 0;
  
  if(!(stepCounter%1000)){
    
    printf("Step: %d, Time: %lg \n",
	   stepCounter, T);
  }

  ++stepCounter;



#if 0
  
  dlong NlocalGrads = mesh.Nelements*mesh.Np*Ngrads;
  dlong NhaloGrads  = mesh.totalHaloPairs*mesh.Np*Ngrads;
  deviceMemory<dfloat> o_gradq = platform.reserve<dfloat>(NlocalGrads+NhaloGrads);

  // extract q trace halo and start exchange
  fieldTraceHalo.ExchangeStart(o_Q, 1);

  // compute volume contributions to gradients
  gradVolumeKernel(mesh.Nelements,
                   mesh.o_vgeo,
                   mesh.o_D,
                   o_Q,
                   o_gradq);

  // complete trace halo exchange
  fieldTraceHalo.ExchangeFinish(o_Q, 1);

  // compute surface contributions to gradients
  gradSurfaceKernel(mesh.Nelements,
                    mesh.o_sgeo,
                    mesh.o_LIFT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    T,
                    mu,
                    gamma,
                    o_Q,
                    o_gradq);

  // extract viscousStresses trace halo and start exchange
  gradTraceHalo.ExchangeStart(o_gradq, 1);

  // compute volume contribution to esdg RHS
  if (cubature) {
    cubatureVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubD,
                         mesh.o_cubPDT,
                         mesh.o_cubInterp,
                         mesh.o_cubProject,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         T,
                         mu,
                         gamma,
                         o_Q,
                         o_gradq,
                         o_RHS);
  } else {
    volumeKernel(mesh.Nelements,
                 mesh.o_vgeo,
                 mesh.o_D,
                 mesh.o_x,
                 mesh.o_y,
                 mesh.o_z,
                 T,
                 mu,
                 gamma,
                 o_Q,
                 o_gradq,
                 o_RHS);
  }

  // complete trace halo exchange
  gradTraceHalo.ExchangeFinish(o_gradq, 1);

  if (cubature) {
      cubatureSurfaceKernel(mesh.Nelements,
                            mesh.o_vgeo,
                            mesh.o_cubsgeo,
                            mesh.o_vmapM,
                            mesh.o_vmapP,
                            mesh.o_EToB,
                            mesh.o_intInterp,
                            mesh.o_intLIFT,
                            mesh.o_intx,
                            mesh.o_inty,
                            mesh.o_intz,
                            T,
                            mu,
                            gamma,
                            o_Q,
                            o_gradq,
                            o_RHS);
    } else {
      surfaceKernel(mesh.Nelements,
                    mesh.o_sgeo,
                    mesh.o_LIFT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    T,
                    mu,
                    gamma,
                    o_Q,
                    o_gradq,
                    o_RHS);
    }
#endif
}
