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

#include "ins.hpp"

//evaluate ODE rhs = f(q,t)
void subcycler_t::rhsf(deviceMemory<dfloat>& o_U, deviceMemory<dfloat>& o_RHS, const dfloat T){

  static int Nsubsteps = 0;

  ++Nsubsteps;
  if(!(Nsubsteps%1000)) std::cout << ", " << "Nsubsteps=" << Nsubsteps << std::endl;
  
  dlong Ntotal = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*NVfields;
  deviceMemory<dfloat> o_Ue = platform.reserve<dfloat>(Ntotal);

  //interpolate velocity history for advective field (halo elements first)
  if(mesh.NhaloElements)
    subCycleAdvectionKernel(mesh.NhaloElements,
                           mesh.o_haloElementIds,
                           shiftIndex,
                           order,
                           maxOrder,
                           mesh.Nelements*mesh.Np*NVfields,
                           T,
                           T0,
                           dt,
                           o_Uh,
                           o_Ue);

  // extract Ue halo
  vTraceHalo.ExchangeStart(o_Ue, 1);

  if(mesh.NinternalElements)
    subCycleAdvectionKernel(mesh.NinternalElements,
                           mesh.o_internalElementIds,
                           shiftIndex,
                           order,
                           maxOrder,
                           mesh.Nelements*mesh.Np*NVfields,
                           T,
                           T0,
                           dt,
                           o_Uh,
                           o_Ue);

  // finish exchange of Ue
  vTraceHalo.ExchangeFinish(o_Ue, 1);

  // (lumped) project Ue to C0
  ins->Project(o_Ue, mesh.dim);

#if 0
  {
    // this checks for C0
    dlong Nlocal = mesh.Nelements*mesh.Np;
    dlong Nhalo  = mesh.totalHaloPairs*mesh.Np;
    memory<dfloat> Ue(mesh.dim*(Nlocal+Nhalo));

    o_Ue.copyTo(Ue);

    for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.Nfp*mesh.Nfaces;++n){
	dlong vidM = mesh.vmapM[e*mesh.Nfp*mesh.Nfaces+n];
	dlong vidP = mesh.vmapP[e*mesh.Nfp*mesh.Nfaces+n];
	if(vidM>=0 && vidP>=0){
	  dlong eM = vidM/mesh.Np;
	  dlong eP = vidP/mesh.Np;
	  dlong nM = vidM%mesh.Np;
	  dlong nP = vidP%mesh.Np;
	  
	  dlong lidM = mesh.Np*2*eM + nM;
	  dlong lidP = mesh.Np*2*eP + nP;
	  
	  dfloat du = Ue[lidM]-Ue[lidP];
	  dfloat dv = Ue[lidM+mesh.Np]-Ue[lidP+mesh.Np];
	  if(du || dv)
	    printf("e: %d, n: %d, (x,y)=(%e,%e), du=%e, dv=%e\n",
		   e, n, mesh.x[vidM], mesh.y[vidM], du, dv);
	}
      }
    }
  }
#endif
  
  
  // extract u halo on DEVICE
  vTraceHalo.ExchangeStart(o_U, 1);

  if (cubature)
    advectionVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubD,
                         mesh.o_cubPDT,
                         mesh.o_cubInterp,
                         mesh.o_cubProject,
                         o_Ue,
                         o_U,
                         o_RHS);
  else
    advectionVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_D,
                         o_Ue,
                         o_U,
                         o_RHS);

  vTraceHalo.ExchangeFinish(o_U, 1);

  if (cubature)
    advectionSurfaceKernel(mesh.Nelements,
                          mesh.o_vgeo,
                          mesh.o_cubsgeo,
                          mesh.o_intInterp,
                          mesh.o_intLIFT,
                          mesh.o_vmapM,
                          mesh.o_vmapP,
                          mesh.o_EToB,
                          T,
                          mesh.o_intx,
                          mesh.o_inty,
                          mesh.o_intz,
                          nu,
                          o_Ue,
                          o_U,
                          o_RHS);
  else
    advectionSurfaceKernel(mesh.Nelements,
                          mesh.o_sgeo,
                          mesh.o_LIFT,
                          mesh.o_vmapM,
                          mesh.o_vmapP,
                          mesh.o_EToB,
                          T,
                          mesh.o_x,
                          mesh.o_y,
                          mesh.o_z,
                          nu,
                          o_Ue,
                          o_U,
                          o_RHS);

#if 0
  relaxationFilterKernel(mesh.Nelements, o_FILT, o_U, o_RHS);
#endif
}
