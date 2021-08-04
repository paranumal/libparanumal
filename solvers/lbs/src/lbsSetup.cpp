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

  //Setup PML
  lbs->PmlSetup();

  //setup timeStepper
  dlong Nlocal = mesh.Nelements*mesh.Np*lbs->Nfields;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np*lbs->Nfields;

  lbs->semiAnalytic = 0;
  if (settings.compareSetting("TIME INTEGRATOR","SARK4")
    ||settings.compareSetting("TIME INTEGRATOR","SARK5")
    ||settings.compareSetting("TIME INTEGRATOR","SAAB3")
    ||settings.compareSetting("TIME INTEGRATOR","MRSAAB3"))
    lbs->semiAnalytic = 1;

  //semi-analytic exponential coefficients
  dfloat lambda[lbs->Nfields];
  for (int i=0;i<mesh.dim+1;i++) lambda[i] = 0.0;
  for (int i=mesh.dim+1;i<lbs->Nfields;i++) lambda[i] = -lbs->tauInv;

  //make array of time step estimates for each element
  dfloat *EtoDT = (dfloat *) calloc(mesh.Nelements,sizeof(dfloat));
  dfloat vmax = lbs->MaxWaveSpeed();
  for(dlong e=0;e<mesh.Nelements;++e){
    dfloat h = mesh.ElementCharacteristicLength(e);
    dfloat dtAdv  = h/(vmax*(mesh.N+1.)*(mesh.N+1.));
    dfloat dtVisc = 1.0/lbs->tauInv;

    if (lbs->semiAnalytic)
      EtoDT[e] = dtAdv;
    else
      EtoDT[e] = mymin(dtAdv, dtVisc);

    /*
    Artificial warping of time step size for multirate testing
    */
#if 0
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
    lbs->multirateTraceHalo = mesh.MultiRateHaloTraceSetup(lbs->Nfields);
  }

  if (settings.compareSetting("TIME INTEGRATOR","MRAB3")){
    lbs->timeStepper = new TimeStepper::mrab3_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, lbs->Nfields, lbs->Npmlfields, *lbs, mesh);
  } else if (settings.compareSetting("TIME INTEGRATOR","MRSAAB3")){
    lbs->timeStepper = new TimeStepper::mrsaab3_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, lbs->Nfields, lbs->Npmlfields, lambda, *lbs, mesh);
  } else if (settings.compareSetting("TIME INTEGRATOR","SAAB3")) {
    lbs->timeStepper = new TimeStepper::saab3_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, lbs->Nfields, lbs->Npmlfields, lambda, *lbs);
  } else if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    lbs->timeStepper = new TimeStepper::ab3_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, lbs->Nfields, lbs->Npmlfields, *lbs);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    lbs->timeStepper = new TimeStepper::lserk4_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, lbs->Nfields, lbs->Npmlfields, *lbs);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    lbs->timeStepper = new TimeStepper::dopri5_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, lbs->Nfields, lbs->Npmlfields, *lbs, mesh.comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","SARK4")) {
    lbs->timeStepper = new TimeStepper::sark4_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, lbs->Nfields, lbs->Npmlfields, lambda, *lbs, mesh.comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","SARK5")) {
    lbs->timeStepper = new TimeStepper::sark5_pml(mesh.Nelements, mesh.NpmlElements, mesh.totalHaloPairs,
                                              mesh.Np, lbs->Nfields, lbs->Npmlfields, lambda, *lbs, mesh.comm);
  } else {
    LIBP_ABORT(string("Requested TIME INTEGRATOR not found."));
  }
  free(EtoDT);

  //setup linear algebra module
  platform.linAlg.InitKernels({"innerProd"});

  /*setup trace halo exchange */
  lbs->traceHalo = mesh.HaloTraceSetup(lbs->Nfields);

  // compute samples of q at interpolation nodes
  lbs->q = (dfloat*) calloc(Nlocal+Nhalo, sizeof(dfloat));
  lbs->o_q = platform.malloc((Nlocal+Nhalo)*sizeof(dfloat), lbs->q);

  lbs->Vort = (dfloat*) calloc(mesh.dim*mesh.Nelements*mesh.Np, sizeof(dfloat));
  lbs->o_Vort = platform.malloc((mesh.dim*mesh.Nelements*mesh.Np)*sizeof(dfloat),
                                              lbs->Vort);

  // Hold macro quantites i.e. density + velocities
  lbs->U = (dfloat*) calloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*lbs->Nmacro, sizeof(dfloat));
  lbs->o_U = platform.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*lbs->Nmacro*sizeof(dfloat), lbs->U);
  
  // Lattice-Boltzmann Model
  lbs->o_LBM = platform.malloc(lbs->Nfields*lbs->Nmacro*sizeof(dfloat), lbs->LBM);
  lbs->o_LMAP = platform.malloc(lbs->Nfields*sizeof(int), lbs->LMAP);


  #if 0

  lbs->o_LBM.copyTo(lbs->LBM);
  
  printf("%.4e\n", lbs->c);
  for(int i=0; i<lbs->Nfields; i++)
  printf("%.4e %.4e %4e \n", lbs->LBM[i+ 0*lbs->Nfields], lbs->LBM[i+ 1*lbs->Nfields],lbs->LBM[i+ 2*lbs->Nfields]);  
  #endif 

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
  kernelInfo["defines/" "p_Npmlfields"]= lbs->Npmlfields;
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

  // int NblockCub = mymax(1, blockMax/mesh.cubNp);
  // kernelInfo["defines/" "p_NblockCub"]= NblockCub;

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

  // if (lbs->pmlcubature) {
  //   sprintf(kernelName, "lbsPmlVolumeCub%s", suffix);
  //   lbs->pmlVolumeKernel =  platform.buildKernel(fileName, kernelName,
  //                                        kernelInfo);
  // } else {
  //   sprintf(kernelName, "lbsPmlVolume%s", suffix);
  //   lbs->pmlVolumeKernel =  platform.buildKernel(fileName, kernelName,
  //                                        kernelInfo);
  // }

  // // kernels from relaxation file
  // sprintf(fileName, DLBS "/okl/lbsRelaxation%s.okl", suffix);
  // sprintf(kernelName, "lbsRelaxation%s", suffix);
  // lbs->relaxationKernel = platform.buildKernel(fileName, kernelName,
  //                                        kernelInfo);
  // if (lbs->pmlcubature) {
  //   sprintf(kernelName, "lbsPmlRelaxationCub%s", suffix);
  //   lbs->pmlRelaxationKernel = platform.buildKernel(fileName, kernelName,
  //                                          kernelInfo);
  // } else {
  //   lbs->pmlRelaxationKernel = platform.buildKernel(fileName, kernelName,
  //                                          kernelInfo);
  // }


  // kernels from surface file
  sprintf(fileName, DLBS "/okl/lbsSurface%s.okl", suffix);
  // if (settings.compareSetting("TIME INTEGRATOR","MRAB3") ||
  //     settings.compareSetting("TIME INTEGRATOR","MRSAAB3")) {
  //   sprintf(kernelName, "lbsMRSurface%s", suffix);
  //   lbs->surfaceKernel = platform.buildKernel(fileName, kernelName,
  //                                          kernelInfo);

  //   sprintf(kernelName, "lbsMRPmlSurface%s", suffix);
  //   lbs->pmlSurfaceKernel = platform.buildKernel(fileName, kernelName,
  //                                          kernelInfo);
  // } else {
    sprintf(kernelName, "lbsSurface%s", suffix);
    lbs->surfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

    // sprintf(kernelName, "lbsPmlSurface%s", suffix);
    // lbs->pmlSurfaceKernel = platform.buildKernel(fileName, kernelName,
    //                                        kernelInfo);
  // }

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