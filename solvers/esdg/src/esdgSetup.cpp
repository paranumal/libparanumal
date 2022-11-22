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

void esdg_t::Setup(platform_t& _platform, mesh_t& _mesh, esdgSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

   
  // port EToE to device
  o_EToE = platform.malloc<dlong>(mesh.Nelements*mesh.Nfaces, mesh.EToE);
  o_EToB = platform.malloc<dlong>(mesh.Nelements*mesh.Nfaces, mesh.EToB);
  
  //get physical parameters
  settings.getSetting("GAMMA", gamma);

  // do element specific set up stff
  switch(mesh.elementType){
  case Mesh::QUADRILATERALS:
    SetupQuad2D();
    break;
  case Mesh::TRIANGLES:
    SetupTri2D();
    break;
  default:
    break;
  }
    
  o_esTotalEntropy = platform.malloc<dfloat>(mesh.Nelements);
  esTotalEntropy.malloc(mesh.Nelements);

  dlong NlocalFields = mesh.Nelements*mesh.Np*Nfields;
  dlong NhaloFields  = mesh.totalHaloPairs*mesh.Np*Nfields;
  //setup timeStepper
  if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    timeStepper.Setup<TimeStepper::ab3>(mesh.Nelements,
					mesh.totalHaloPairs,
					mesh.Np, Nfields, platform, comm);
    
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    timeStepper.Setup<TimeStepper::lserk4>(mesh.Nelements,
					   mesh.totalHaloPairs,
					   mesh.Np, Nfields, platform, comm);
    
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    timeStepper.Setup<TimeStepper::dopri5>(mesh.Nelements,
					   mesh.totalHaloPairs,
					   mesh.Np, Nfields, platform, comm);
  }

  // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat cfl = 0.02; // depends on the stability region size

  dfloat dtAdv  = hmin/((mesh.N+1.)*(mesh.N+1.)*gamma); //just an estimate

  dfloat dt = cfl*dtAdv;
  timeStepper.SetTimeStep(dt);

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd"});

  /*setup trace halo exchange */
  fieldTraceHalo = mesh.HaloTraceSetup(Nfields);
  gradTraceHalo  = mesh.HaloTraceSetup(Ngrads);

  // compute samples of q at interpolation nodes
  q.malloc(NlocalFields+NhaloFields);
  o_q = platform.malloc<dfloat>((NlocalFields+NhaloFields),
				  q);

  qexact.malloc(NlocalFields+NhaloFields);
  o_qexact = platform.malloc<dfloat>((NlocalFields+NhaloFields),
				       q);

  Vort.malloc(mesh.dim*mesh.Nelements*mesh.Np);
  o_Vort = platform.malloc<dfloat>((mesh.dim*mesh.Nelements*mesh.Np),
				     Vort);

  //storage for M*q during reporting
  o_Mq = platform.malloc<dfloat>((NlocalFields+NhaloFields), q);

  // linf from error kernel
  linfError.malloc(Nfields*mesh.Nelements);
  l2Error.malloc(Nfields*mesh.Nelements);
  
  o_linfError = platform.malloc<dfloat>(Nfields*mesh.Nelements);
  o_l2Error   = platform.malloc<dfloat>(Nfields*mesh.Nelements);

  o_cubw = platform.malloc<dfloat>(mesh.cubNp, mesh.cubw);
  
}

esdg_t::~esdg_t() {

}
