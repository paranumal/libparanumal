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

#include "maxwell.hpp"

dfloat maxwell_t::MaxWaveSpeed(){
  //wavespeed is constant 1 everywhere
  const dfloat vmax = 1.0;
  return vmax;
}

//evaluate ODE rhs = f(q,t)
void maxwell_t::rhsf(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

  // extract q halo on DEVICE
  traceHalo.ExchangeStart(o_Q, 1);

  rhsVolume(mesh.NnonPmlElements, mesh.o_nonPmlElements, o_Q, o_RHS);

  rhsSurface(mesh.NinternalElements, mesh.o_internalElementIds, o_Q, o_RHS, T);

  traceHalo.ExchangeFinish(o_Q, 1);

  rhsSurface(mesh.NhaloElements, mesh.o_haloElementIds, o_Q, o_RHS, T); 
  
  if(materialType==HETEROGENEOUS){    
    heterogeneousProjectKernel(mesh.Nelements, mesh.o_cubInterp, mesh.o_cubProject, o_materialInverseWeights, o_RHS);
  }
}


//evaluate ODE rhs = f(q,t)
void maxwell_t::rhsf_pml(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_pmlQ,
			 deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_pmlRHS, const dfloat T){

  // extract q trace halo and start exchange
  traceHalo.ExchangeStart(o_Q, 1);

  // compute volume contribution to maxwell RHS
  rhsVolume(mesh.NnonPmlElements, mesh.o_nonPmlElements, o_Q, o_RHS);
  
  rhsPmlVolume(mesh.NpmlElements, mesh.o_pmlElements, mesh.o_pmlIds,
               o_Q, o_pmlQ, o_RHS, o_pmlRHS);
  
  // complete trace halo exchange
  traceHalo.ExchangeFinish(o_Q, 1);
  
  // compute surface contribution to maxwell RHS
  rhsSurface(mesh.NnonPmlElements, mesh.o_nonPmlElements, o_Q, o_RHS, T);
  rhsPmlSurface(mesh.NpmlElements, mesh.o_pmlElements, mesh.o_pmlIds,
                o_Q, o_pmlQ, o_RHS, o_pmlRHS, T);
  

  if(materialType==HETEROGENEOUS){
    //printf("Calling heterogeneousProjectKernel\n");
    heterogeneousProjectKernel(mesh.Nelements, mesh.o_cubInterp, mesh.o_cubProject, o_materialInverseWeights, o_RHS);
  }

  if(mesh.NpmlElements){
    if(pmlcubature){
      //printf("Calling pmlCubatureTermsKernel\n");
      pmlCubatureTermsKernel(mesh.NpmlElements, mesh.o_pmlElements, mesh.o_pmlIds,
			     o_pmlSigma, o_pmlBeta, mesh.o_cubInterp, mesh.o_cubProject, o_Q, o_pmlQ, o_RHS, o_pmlRHS);
    }else{
      pmlTermsKernel(mesh.NpmlElements, mesh.o_pmlElements, mesh.o_pmlIds,
		     o_pmlSigma, o_pmlBeta, o_Q, o_pmlQ, o_RHS, o_pmlRHS);
    }

#if 0
    if(pmlcubature){
      //printf("Calling anisotropicCubatureProjectKernel\n");
      anisotropicCubatureProjectKernel(mesh.NpmlElements, mesh.o_pmlElements, mesh.o_pmlIds, mesh.o_cubInterp, mesh.o_cubProject, o_pmlInvWeights, o_RHS);
    }else{
      //      printf("Calling anisotropicProjectKernel\n");
      anisotropicProjectKernel(mesh.NpmlElements, mesh.o_pmlElements, mesh.o_pmlIds, o_pmlInvWeights, o_RHS);
    }
#else
    dfloat TOL = 1e-8;
    int maxNit = 100;
    massMatrixSolveKernel(TOL, maxNit, mesh.NpmlElements, mesh.o_pmlElements, mesh.o_pmlIds, mesh.o_MM, mesh.o_cubInterp,  o_pmlWeights, o_RHS);
#endif
  }
}

void maxwell_t::rhsVolume(dlong N, deviceMemory<dlong>& o_ids,
			  deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS){
  
  // compute volume contribution to maxwell RHS
  if (N)
    volumeKernel(N,
                 o_ids,
                 mesh.o_vgeo,
                 mesh.o_D,
                 o_Q,
                 o_RHS);
}

void maxwell_t::rhsSurface(dlong N, deviceMemory<dlong>& o_ids,
			   deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

  // compute volume contribution to maxwell RHS
  if (N){
    if(materialType==ISOTROPIC  || materialType==ANISOTROPIC){
      surfaceKernel(N,
		    o_ids,
		    mesh.o_sgeo,
		    mesh.o_LIFT,
		    mesh.o_vmapM,
		    mesh.o_vmapP,
		    mesh.o_EToB,
		    T,
		    mesh.o_x,
		    mesh.o_y,
		    mesh.o_z,
		    o_Q,
		    o_RHS);
    }
    else{
      heterogeneousSurfaceKernel(N, 
				 o_ids,
				 mesh.o_sgeo,
				 mesh.o_LIFT,
				 mesh.o_vmapM,
				 mesh.o_vmapP,
				 mesh.o_EToB,
				 T,
				 mesh.o_x,
				 mesh.o_y,
				 mesh.o_z,
				 mesh.o_intInterp,
				 mesh.o_intLIFT,
				 o_materialUpwindWeights,
				 o_Q,
				 o_RHS);
      
    }
  }
}

void maxwell_t::rhsPmlVolume(dlong N, deviceMemory<dlong>& o_ids, deviceMemory<dlong>& o_pmlids,
			     deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_pmlQ,
			     deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_pmlRHS){
  
  // compute volume contribution to maxwell RHS
  if (N) {
    pmlVolumeKernel(N,
		    o_ids,
		    o_pmlids,
		    mesh.o_vgeo,
		    mesh.o_D,
		    o_Q,
		    o_pmlQ,
		    o_RHS,
		    o_pmlRHS);
  }
}


void maxwell_t::rhsPmlSurface(dlong N, deviceMemory<dlong>& o_ids, deviceMemory<dlong>& o_pmlids,
			      deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_pmlQ,
			      deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_pmlRHS, const dfloat T){

  // compute volume contribution to maxwell RHS
  if (N)
    pmlSurfaceKernel(N,
		     o_ids,
		     o_pmlids,
		     mesh.o_sgeo,
		     mesh.o_LIFT,
		     mesh.o_vmapM,
		     mesh.o_vmapP,
		     mesh.o_EToB,
		     T,
		     mesh.o_x,
		     mesh.o_y,
		     mesh.o_z,
		     o_Q,
		     o_RHS,
		     o_pmlRHS);
}
