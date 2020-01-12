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

#include "bns.hpp"

bns_t& bns_t::Setup(mesh_t& mesh, linAlg_t& linAlg,
                    bnsSettings_t& settings){

  bns_t* bns = new bns_t(mesh, linAlg, settings);

  //get physical paramters
  settings.getSetting("SPEED OF SOUND", bns->c);
  settings.getSetting("VISCOSITY", bns->nu);
  bns->RT     = bns->c*bns->c;
  bns->tauInv = bns->RT/bns->nu;

  bns->Nfields    = (mesh.dim==3) ? 10:6;
  bns->Npmlfields = mesh.dim*bns->Nfields;

  //setup cubature
  mesh.CubatureSetup();

  //Setup PML
  bns->PmlSetup();

  //setup timeStepper
  dlong Nlocal = mesh.Nelements*mesh.Np*bns->Nfields;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np*bns->Nfields;

  bns->semiAnalytic = 0;
  if (settings.compareSetting("TIME INTEGRATOR","SARK4")
    ||settings.compareSetting("TIME INTEGRATOR","SARK5")
    ||settings.compareSetting("TIME INTEGRATOR","SAAB3")
    ||settings.compareSetting("TIME INTEGRATOR","MRSAAB3"))
    bns->semiAnalytic = 1;

  //semi-analytic exponential coefficients
  dfloat lambda[bns->Nfields];
  for (int i=0;i<mesh.dim+1;i++) lambda[i] = 0.0;
  for (int i=mesh.dim+1;i<bns->Nfields;i++) lambda[i] = -bns->tauInv;

  //make array of time step estimates for each element
  dfloat *EtoDT = (dfloat *) calloc(mesh.Nelements,sizeof(dfloat));
  for(dlong e=0;e<mesh.Nelements;++e){
    dfloat dtAdv  = mesh.ElementCharacteristicLength(e)/((mesh.N+1.)*(mesh.N+1.)*sqrt(3.0)*bns->c); //just an estimate
    dfloat dtVisc = 1.0/bns->tauInv;

    if (bns->semiAnalytic)
      EtoDT[e] = dtAdv;
    else
      EtoDT[e] = mymin(dtAdv, dtVisc);

    /*
    Artificial warping of time step size for multirate testing
    */
#if 1
    dfloat x=0., y=0, z=0;
    for (int v=0;v<mesh.Nverts;v++){
      x += mesh.EX[e*mesh.Nverts+v];
      y += mesh.EY[e*mesh.Nverts+v];
      if (mesh.dim==3)
        z += mesh.EZ[e*mesh.Nverts+v];
    }
    x /=mesh.Nverts;
    y /=mesh.Nverts;
    z /=mesh.Nverts;

    dfloat c = mymin(fabs(x),fabs(y));
    if (mesh.dim==3)
      c = mymin(c,fabs(z));

    c = mymax(0.5, c);
    EtoDT[e] *= c;
#endif
  }

  mesh.mrNlevels=0;
  if (settings.compareSetting("TIME INTEGRATOR","MRAB3") ||
      settings.compareSetting("TIME INTEGRATOR","MRSAAB3")) {
    mesh.MultiRateSetup(EtoDT);
    mesh.MultiRatePmlSetup();
    bns->multirateTraceHalo = mesh.MultiRateHaloTraceSetup(bns->Nfields);
  }

  if (settings.compareSetting("TIME INTEGRATOR","MRAB3")){
    bns->timeStepper = new TimeStepper::mrab3_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, bns->Nfields, bns->Npmlfields, *bns);
  } else if (settings.compareSetting("TIME INTEGRATOR","MRSAAB3")){
    bns->timeStepper = new TimeStepper::mrsaab3_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, bns->Nfields, bns->Npmlfields, lambda, *bns);
  } else if (settings.compareSetting("TIME INTEGRATOR","SAAB3")) {
    bns->timeStepper = new TimeStepper::saab3_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, bns->Nfields, bns->Npmlfields, lambda, *bns);
  } else if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    bns->timeStepper = new TimeStepper::ab3_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, bns->Nfields, bns->Npmlfields, *bns);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    bns->timeStepper = new TimeStepper::lserk4_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, bns->Nfields, bns->Npmlfields, *bns);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    bns->timeStepper = new TimeStepper::dopri5_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, bns->Nfields, bns->Npmlfields, *bns);
  } else if (settings.compareSetting("TIME INTEGRATOR","SARK4")) {
    bns->timeStepper = new TimeStepper::sark4_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, bns->Nfields, bns->Npmlfields, lambda, *bns);
  } else if (settings.compareSetting("TIME INTEGRATOR","SARK5")) {
    bns->timeStepper = new TimeStepper::sark5_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, bns->Nfields, bns->Npmlfields, lambda, *bns);
  } else {
    LIBP_ABORT(string("Requested TIME INTEGRATOR not found."));
  }
  free(EtoDT);

  // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat cfl = 0.4; // depends on the stability region size

  dfloat dtAdv  = hmin/((mesh.N+1.)*(mesh.N+1.)*sqrt(3.0)*bns->c); //just an estimate
  dfloat dtVisc = 1.0/bns->tauInv;

  dfloat dt;
  if (bns->semiAnalytic)
    dt = cfl*dtAdv;
  else
    dt = cfl*mymin(dtAdv, dtVisc);

/*
    Artificial warping of time step size for multirate testing
    */
#if 1
  if (settings.compareSetting("TIME INTEGRATOR","MRAB3") ||
      settings.compareSetting("TIME INTEGRATOR","MRSAAB3"))
    dt /= (1<<(mesh.mrNlevels-1));
#endif

  bns->timeStepper->SetTimeStep(dt);

  //setup linear algebra module
  bns->linAlg.InitKernels({"innerProd"}, mesh.comm);

  /*setup trace halo exchange */
  bns->traceHalo = mesh.HaloTraceSetup(bns->Nfields);

  // compute samples of q at interpolation nodes
  bns->q = (dfloat*) calloc(Nlocal+Nhalo, sizeof(dfloat));
  bns->o_q = mesh.device.malloc((Nlocal+Nhalo)*sizeof(dfloat), bns->q);

  bns->Vort = (dfloat*) calloc(mesh.dim*mesh.Nelements*mesh.Np, sizeof(dfloat));
  bns->o_Vort = mesh.device.malloc((mesh.dim*mesh.Nelements*mesh.Np)*sizeof(dfloat),
                                              bns->Vort);

  //storage for M*q during reporting
  bns->o_Mq = mesh.device.malloc((Nlocal+Nhalo)*sizeof(dfloat), bns->q);

  // OCCA build stuff
  occa::properties kernelInfo = bns->props; //copy base occa properties

  //add boundary data to kernel info
  string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= bns->Nfields;
  kernelInfo["defines/" "p_Npmlfields"]= bns->Npmlfields;

  int maxNodes = mymax(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = 512/mesh.Np; // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  int NblockCub = 512/mesh.cubNp; // works for CUDA
  kernelInfo["defines/" "p_NblockCub"]= NblockCub;

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


  if(mesh.elementType==HEXAHEDRA && bns->pmlcubature)
    LIBP_ABORT("PML CUBATURE Not currently supported with Hex meshes.")

  // kernels from volume file
  sprintf(fileName, DBNS "/okl/bnsVolume%s.okl", suffix);
  sprintf(kernelName, "bnsVolume%s", suffix);
  bns->volumeKernel =  buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);

  if (bns->pmlcubature) {
    sprintf(kernelName, "bnsPmlVolumeCub%s", suffix);
    bns->pmlVolumeKernel =  buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);
  } else {
    sprintf(kernelName, "bnsPmlVolume%s", suffix);
    bns->pmlVolumeKernel =  buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);
  }

  // kernels from relaxation file
  sprintf(fileName, DBNS "/okl/bnsRelaxation%s.okl", suffix);
  sprintf(kernelName, "bnsRelaxation%s", suffix);
  bns->relaxationKernel = buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);
  if (bns->pmlcubature) {
    sprintf(kernelName, "bnsPmlRelaxationCub%s", suffix);
    bns->pmlRelaxationKernel = buildKernel(mesh.device, fileName, kernelName,
                                           kernelInfo, mesh.comm);
  } else {
    bns->pmlRelaxationKernel = buildKernel(mesh.device, fileName, kernelName,
                                           kernelInfo, mesh.comm);
  }


  // kernels from surface file
  sprintf(fileName, DBNS "/okl/bnsSurface%s.okl", suffix);
  if (settings.compareSetting("TIME INTEGRATOR","MRAB3") ||
      settings.compareSetting("TIME INTEGRATOR","MRSAAB3")) {
    sprintf(kernelName, "bnsMRSurface%s", suffix);
    bns->surfaceKernel = buildKernel(mesh.device, fileName, kernelName,
                                           kernelInfo, mesh.comm);

    sprintf(kernelName, "bnsMRPmlSurface%s", suffix);
    bns->pmlSurfaceKernel = buildKernel(mesh.device, fileName, kernelName,
                                           kernelInfo, mesh.comm);
  } else {
    sprintf(kernelName, "bnsSurface%s", suffix);
    bns->surfaceKernel = buildKernel(mesh.device, fileName, kernelName,
                                           kernelInfo, mesh.comm);

    sprintf(kernelName, "bnsPmlSurface%s", suffix);
    bns->pmlSurfaceKernel = buildKernel(mesh.device, fileName, kernelName,
                                           kernelInfo, mesh.comm);
  }


  // mass matrix operator
  sprintf(fileName, LIBP_DIR "/core/okl/MassMatrixOperator%s.okl", suffix);
  sprintf(kernelName, "MassMatrixOperator%s", suffix);

  bns->MassMatrixKernel = buildKernel(mesh.device, fileName, kernelName,
                                      kernelInfo, mesh.comm);

  // vorticity calculation
  sprintf(fileName, DBNS "/okl/bnsVorticity%s.okl", suffix);
  sprintf(kernelName, "bnsVorticity%s", suffix);

  bns->vorticityKernel = buildKernel(mesh.device, fileName, kernelName,
                                     kernelInfo, mesh.comm);

  if (mesh.dim==2) {
    sprintf(fileName, DBNS "/okl/bnsInitialCondition2D.okl");
    sprintf(kernelName, "bnsInitialCondition2D");
  } else {
    sprintf(fileName, DBNS "/okl/bnsInitialCondition3D.okl");
    sprintf(kernelName, "bnsInitialCondition3D");
  }

  bns->initialConditionKernel = buildKernel(mesh.device, fileName, kernelName,
                                            kernelInfo, mesh.comm);

  return *bns;
}

bns_t::~bns_t() {
  volumeKernel.free();
  surfaceKernel.free();
  relaxationKernel.free();
  pmlVolumeKernel.free();
  pmlSurfaceKernel.free();
  pmlRelaxationKernel.free();
  vorticityKernel.free();
  MassMatrixKernel.free();
  initialConditionKernel.free();

  if (timeStepper) delete timeStepper;
  if (traceHalo) traceHalo->Free();

  for (int lev=0;lev<mesh.mrNlevels;lev++)
    if (multirateTraceHalo[lev]) multirateTraceHalo[lev]->Free();
}