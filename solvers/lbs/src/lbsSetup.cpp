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

#include "lbs.hpp"
#define D2Q9 1
#define D3Q15 2


lbs_t& lbs_t::Setup(platform_t& platform, mesh_t& mesh,
                    lbsSettings_t& settings){

  lbs_t* lbs = new lbs_t(platform, mesh, settings);

  // Set reference lattice-Boltzmann data  
  lbs->latticeSetup();
  
  lbs->Npmlfields = mesh.dim*lbs->Nfields;

  // AK: not in use yet ... Setup PML
  // lbs->PmlSetup();

  //setup timeStepper
  dlong Nlocal = mesh.Nelements*mesh.Np*lbs->Nfields;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np*lbs->Nfields;

  //make array of time step estimates for each element
  dfloat *EtoDT = (dfloat *) calloc(mesh.Nelements,sizeof(dfloat));
  dfloat vmax = lbs->MaxWaveSpeed();
  for(dlong e=0;e<mesh.Nelements;++e){
    dfloat h = mesh.ElementCharacteristicLength(e);
    dfloat dtAdv  = h/(vmax*(mesh.N+1.)*(mesh.N+1.));
    EtoDT[e] = dtAdv;
  }



  if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    lbs->timeStepper = new TimeStepper::lserk4(mesh.Nelements, mesh.totalHaloPairs,
					       mesh.Np, lbs->Nfields, *lbs);
  }else {
    LIBP_ABORT(string("Requested TIME INTEGRATOR not found."));
  }

  
  //setup linear algebra module
  platform.linAlg.InitKernels({"innerProd"});

  /*setup trace halo exchange */
  lbs->traceHalo = mesh.HaloTraceSetup(lbs->Nfields);

  // compute samples of q at interpolation nodes
  lbs->q = (dfloat*) calloc(Nlocal+Nhalo, sizeof(dfloat));
  lbs->o_q = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat), lbs->q);

  lbs->F = (dfloat*) calloc(Nlocal+Nhalo, sizeof(dfloat));
  lbs->o_F = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat), lbs->F);


  lbs->Vort = (dfloat*) calloc(mesh.dim*mesh.Nelements*mesh.Np, sizeof(dfloat));
  lbs->o_Vort = platform.malloc((mesh.dim*mesh.Nelements*mesh.Np)*sizeof(dfloat),
				lbs->Vort);

  // Hold macro quantites i.e. density + velocities
  lbs->U = (dfloat*) calloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*lbs->Nmacro, sizeof(dfloat));
  lbs->o_U = platform.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*lbs->Nmacro*sizeof(dfloat), lbs->U);
  
  // Lattice-Boltzmann Model
  lbs->o_LBM = platform.malloc(lbs->Nfields*lbs->Nmacro*sizeof(dfloat), lbs->LBM);
  lbs->o_LMAP = platform.malloc(lbs->Nfields*sizeof(int), lbs->LMAP);


  
  //storage for M*q during reporting
  lbs->o_Mq = platform.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*lbs->Nmacro*sizeof(dfloat), lbs->U);
  mesh.MassMatrixKernelSetup(lbs->Nmacro); // mass matrix operator

  // // OCCA build stuff
  occa::properties kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= lbs->Nfields;
  // kernelInfo["defines/" "p_Npmlfields"]= lbs->Npmlfields;
  kernelInfo["defines/" "p_Nmacro"] = lbs->Nmacro;  

  kernelInfo["defines/" "p_c"] = lbs->c;  
  kernelInfo["defines/" "p_ic2"] = 1.0/ pow(lbs->c,2);  
  kernelInfo["defines/" "p_ic4"] = 1.0/ pow(lbs->c,4);   

  int maxNodes = mymax(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode()=="CUDA") blockMax = 512;

  int NblockV = mymax(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = mymax(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  // set kernel name suffix
  char *suffix;
  if(mesh.elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(mesh.elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(mesh.elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(mesh.elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  if (mesh.dim==2) {
    sprintf(fileName, DLBS "/okl/lbsInitialCondition2D.okl");
    sprintf(kernelName, "lbsInitialCondition2D");
  } else {
    sprintf(fileName, DLBS "/okl/lbsInitialCondition3D.okl");
    sprintf(kernelName, "lbsInitialCondition3D");
  }
  lbs->initialConditionKernel = platform.buildKernel(fileName, kernelName,
						     kernelInfo);
  
  // kernels from volume file
  sprintf(fileName, DLBS "/okl/lbsCollision%s.okl", suffix);  

  sprintf(kernelName, "lbsCollision%s", suffix);
  lbs->collisionKernel =  platform.buildKernel(fileName, kernelName,
					       kernelInfo);

  sprintf(kernelName, "lbsForcing%s", suffix);
  lbs->forcingKernel =  platform.buildKernel(fileName, kernelName,
					     kernelInfo);

  sprintf(kernelName, "lbsMoments%s", suffix);
  lbs->momentsKernel =  platform.buildKernel(fileName, kernelName,
					     kernelInfo);

  sprintf(kernelName, "lbsPhaseField%s", suffix);
  lbs->phaseFieldKernel =  platform.buildKernel(fileName, kernelName,
						kernelInfo);

  // kernels from volume file
  sprintf(fileName, DLBS "/okl/lbsVolume%s.okl", suffix);
  sprintf(kernelName, "lbsVolume%s", suffix);
  lbs->volumeKernel =  platform.buildKernel(fileName, kernelName,
					    kernelInfo);

  // kernels from surface file
  sprintf(fileName, DLBS "/okl/lbsSurface%s.okl", suffix);
  
  sprintf(kernelName, "lbsSurface%s", suffix);
  lbs->surfaceKernel = platform.buildKernel(fileName, kernelName,
					    kernelInfo);

  // vorticity calculation
  sprintf(fileName, DLBS "/okl/lbsVorticity%s.okl", suffix);
  sprintf(kernelName, "lbsVorticity%s", suffix);

  lbs->vorticityKernel = platform.buildKernel(fileName, kernelName,
					      kernelInfo);



  return *lbs;
}

lbs_t::~lbs_t() {
  volumeKernel.free();
  surfaceKernel.free();
  relaxationKernel.free();
  pmlVolumeKernel.free();
  pmlSurfaceKernel.free();
  pmlRelaxationKernel.free();
  vorticityKernel.free();
  initialConditionKernel.free();

  if (timeStepper) delete timeStepper;
  if (traceHalo) traceHalo->Free();

  for (int lev=0;lev<mesh.mrNlevels;lev++)
    if (multirateTraceHalo[lev]) multirateTraceHalo[lev]->Free();
}
