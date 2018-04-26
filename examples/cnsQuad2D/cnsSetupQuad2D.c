#include "cnsQuad2D.h"

cns_t *cnsSetupQuad2D(mesh2D *mesh, setupAide &newOptions, char* boundaryHeaderFileName){

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  long int hostId = gethostid();

  long int* hostIds = (long int*) calloc(size,sizeof(long int));
  MPI_Allgather(&hostId,1,MPI_LONG,hostIds,1,MPI_LONG,MPI_COMM_WORLD);

  int deviceID = 0;
  for (int r=0;r<rank;r++) {
    if (hostIds[r]==hostId) deviceID++;
  }

  // read thread model/device/platform from newOptions
  if(newOptions.compareArgs("THREAD MODEL", "CUDA")){
    sprintf(deviceConfig, "mode = CUDA, deviceID = %d", deviceID);
  }
  else if(newOptions.compareArgs("THREAD MODEL", "OpenCL")){
    int plat;
    newOptions.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "mode = OpenCL, deviceID = %d, platformID = %d",
	    deviceID, plat);
  }
  else if(newOptions.compareArgs("THREAD MODEL", "OpenMP")){
    sprintf(deviceConfig, "mode = OpenMP");
  }
  else{
    sprintf(deviceConfig, "mode = Serial");
  }

  int check;
  
  cns_t *cns = (cns_t*) calloc(1, sizeof(cns_t));

  mesh->Nfields = 3;
  cns->Nfields = mesh->Nfields;
  
  cns->Nstresses = 3;
  cns->mesh = mesh;

  dlong Ntotal = mesh->Nelements*mesh->Np*mesh->Nfields;
  cns->Nblock = (Ntotal+blockSize-1)/blockSize;
  
  hlong localElements = (hlong) mesh->Nelements;
  MPI_Allreduce(&localElements, &(cns->totalElements), 1, MPI_HLONG, MPI_SUM, MPI_COMM_WORLD);

  // speed of sound (assuming isothermal unit bulk flow) = sqrt(RT)
  cns->RT = 1.3;

  // viscosity
  cns->mu = 2e-5;

  // mean flow
  cns->rbar = 1;
  cns->ubar = 0.2;
  cns->vbar = 0;

  check = newOptions.getArgs("RBAR", cns->rbar);
  if(!check) printf("WARNING setup file does not include RBAR\n");
  
  check = newOptions.getArgs("UBAR", cns->ubar);
  if(!check) printf("WARNING setup file does not include UBAR\n");

  check = newOptions.getArgs("VBAR", cns->vbar);
  if(!check) printf("WARNING setup file does not include VBAR\n");

  check = newOptions.getArgs("VISCOSITY", cns->mu);
  if(!check) printf("WARNING setup file does not include VISCOSITY\n");

  dfloat mach = 0.17;
  check = newOptions.getArgs("MACH NUMBER", mach);
  if(!check) printf("WARNING setup file does not include MACH\n");

  // speed of sound (assuming isothermal unit bulk flow) = sqrt(RT)
  cns->RT = cns->ubar*cns->ubar/(mach*mach);
  
  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields, sizeof(dfloat));
  cns->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  if (newOptions.compareArgs("TIME INTEGRATOR","LSERK4")){
    cns->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
		  		sizeof(dfloat));
  }

  if (newOptions.compareArgs("TIME INTEGRATOR","DOPRI5")){
    int NrkStages = 7;
    cns->rkq  = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
          sizeof(dfloat));
    cns->rkrhsq = (dfloat*) calloc(NrkStages*mesh->Nelements*mesh->Np*mesh->Nfields,
          sizeof(dfloat));
    cns->rkerr  = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
          sizeof(dfloat));

    cns->errtmp = (dfloat*) calloc(cns->Nblock, sizeof(dfloat));

    // Dormand Prince -order (4) 5 with PID timestep control
    int Nrk = 7;
    dfloat rkC[7] = {0.0, 0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0};
    dfloat rkA[7*7] ={             0.0,             0.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                                   0.2,             0.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                              3.0/40.0,        9.0/40.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                             44.0/45.0,      -56.0/15.0,       32.0/9.0,          0.0,             0.0,       0.0, 0.0,
                        19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0,             0.0,       0.0, 0.0,
                         9017.0/3168.0,     -355.0/33.0, 46732.0/5247.0,   49.0/176.0, -5103.0/18656.0,       0.0, 0.0, 
                            35.0/384.0,             0.0,   500.0/1113.0,  125.0/192.0,  -2187.0/6784.0, 11.0/84.0, 0.0 };
    dfloat rkE[7] = {71.0/57600.0,  0.0, -71.0/16695.0, 71.0/1920.0, -17253.0/339200.0, 22.0/525.0, -1.0/40.0 }; 

    cns->Nrk = Nrk;
    cns->rkC = (dfloat*) calloc(cns->Nrk, sizeof(dfloat));
    cns->rkE = (dfloat*) calloc(cns->Nrk, sizeof(dfloat));
    cns->rkA = (dfloat*) calloc(cns->Nrk*cns->Nrk, sizeof(dfloat));

    memcpy(cns->rkC, rkC, cns->Nrk*sizeof(dfloat));
    memcpy(cns->rkE, rkE, cns->Nrk*sizeof(dfloat));
    memcpy(cns->rkA, rkA, cns->Nrk*cns->Nrk*sizeof(dfloat));
    
    cns->dtMIN = 1E-7; //minumum allowed timestep
    cns->ATOL = 1E-7;  //absolute error tolerance
    cns->RTOL = 1E-5;  //relative error tolerance
    cns->safe = 0.9;   //safety factor

    //error control parameters
    cns->beta = 0.05;
    cns->factor1 = 0.2;
    cns->factor2 = 10.0;

    cns->exp1 = 0.2 - 0.75*cns->beta;
    cns->invfactor1 = 1.0/cns->factor1;
    cns->invfactor2 = 1.0/cns->factor2;
    cns->facold = 1E-4;
    
  }

  
  cns->viscousStresses = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*cns->Nstresses,
					   sizeof(dfloat));
  
  cns->Vort = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  
  // fix this later (initial conditions)
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat t = 0;
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];

      dlong qbase = e*mesh->Np*mesh->Nfields + n;

#if 0
      cnsGaussianPulse2D(x, y, t,
			 mesh->q+qbase,
			 mesh->q+qbase+mesh->Np,
			 mesh->q+qbase+2*mesh->Np);
#else
      mesh->q[qbase+0*mesh->Np] = cns->rbar;
      //      mesh->q[qbase+1*mesh->Np] = cns->rbar*cns->ubar*y*(6-y)/9.;
      const dfloat alpha = 20;
      mesh->q[qbase+1*mesh->Np] = cns->rbar*cns->ubar;//*tanh(alpha*y)*tanh(alpha*(6-y));
      mesh->q[qbase+2*mesh->Np] = cns->rbar*cns->vbar;
#endif
    }
  }

  // set penalty parameter
  mesh->Lambda2 = 0.5;

  // set time step
  dfloat hmin = 1e9;
  for(dlong e=0;e<mesh->Nelements;++e){  

    for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      dlong sid = mesh->Nsgeo*(mesh->Nfaces*mesh->Nfp*e + n);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      if(invJ<0) printf("invJ = %g\n", invJ);
      
      // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
      // h = 0.5/(sJ/J)
      
      dfloat hest = .5/(sJ*invJ);

      hmin = mymin(hmin, hest);
    }
  }

  // need to change cfl and defn of dt
  dfloat cfl = 1; // depends on the stability region size

  dfloat dtAdv  = hmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(cns->RT));
  dfloat dtVisc = pow(hmin, 2)/(pow(mesh->N+1,4)*cns->mu);

  dfloat dt = cfl*mymin(dtAdv, dtVisc);
  dt = cfl*dtAdv;
  
  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  //

  newOptions.getArgs("FINAL TIME", mesh->finalTime);

  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  if (newOptions.compareArgs("TIME INTEGRATOR","LSERK4")){
    mesh->dt = mesh->finalTime/mesh->NtimeSteps;
  }

  if (rank ==0)
    printf("dtAdv = %lg (before cfl), dtVisc = %lg (before cfl), dt = %lg\n",
	   dtAdv, dtVisc, dt);

  cns->frame = 0;
  
  // errorStep
  mesh->errorStep = 1000;

  if (rank ==0) printf("dt = %g\n", mesh->dt);

  
  // OCCA build stuff
  
  occa::kernelInfo kernelInfo;
  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

  //add boundary data to kernel info
  kernelInfo.addInclude(boundaryHeaderFileName);
 
  cns->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);

  cns->o_saveq =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);
  
  cns->o_viscousStresses =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*cns->Nstresses*sizeof(dfloat),
			cns->viscousStresses);
  
  cns->o_rhsq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), cns->rhsq);

  cout << "TIME INTEGRATOR (" <<
    newOptions.getArgs("TIME INTEGRATOR") << ")" << endl;

  if (newOptions.compareArgs("TIME INTEGRATOR","LSERK4")){
    cns->o_resq =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), cns->resq);
  }

  if (newOptions.compareArgs("TIME INTEGRATOR","DOPRI5")){
    printf("setting up DOPRI5\n");
    int NrkStages = 7;
    cns->o_rkq =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), cns->rkq);
    cns->o_rkrhsq =
      mesh->device.malloc(NrkStages*mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), cns->rkrhsq);
    cns->o_rkerr =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), cns->rkerr);
  
    cns->o_errtmp = mesh->device.malloc(cns->Nblock*sizeof(dfloat), cns->errtmp);

    cns->o_rkA = mesh->device.malloc(cns->Nrk*cns->Nrk*sizeof(dfloat), cns->rkA);
    cns->o_rkE = mesh->device.malloc(  cns->Nrk*sizeof(dfloat), cns->rkE);
  }
  
  cns->o_Vort = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), cns->Vort);
  
  mesh->o_LIFT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
			mesh->LIFT);

  cns->LIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Nfp*mesh->Nfaces;++m){
      cns->LIFTT[n + m*mesh->Np] = mesh->LIFT[n*mesh->Nfaces*mesh->Nfp+m];
    }
  }
  
  cns->o_LIFTT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
			cns->LIFTT);
  

  if(mesh->totalHaloPairs>0){
    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat));

    cns->o_haloStressesBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*cns->Nstresses*sizeof(dfloat));

    // MPI send buffer
    cns->haloBytes = mesh->totalHaloPairs*mesh->Np*cns->Nfields*sizeof(dfloat);
    occa::memory o_sendBuffer = mesh->device.mappedAlloc(cns->haloBytes, NULL);
    occa::memory o_recvBuffer = mesh->device.mappedAlloc(cns->haloBytes, NULL);
    cns->o_haloBuffer = mesh->device.malloc(cns->haloBytes);
    cns->sendBuffer = (dfloat*) o_sendBuffer.getMappedPointer();
    cns->recvBuffer = (dfloat*) o_recvBuffer.getMappedPointer();

    cns->haloStressesBytes = mesh->totalHaloPairs*mesh->Np*cns->Nstresses*sizeof(dfloat);
    occa::memory o_sendStressesBuffer = mesh->device.mappedAlloc(cns->haloStressesBytes, NULL);
    occa::memory o_recvStressesBuffer = mesh->device.mappedAlloc(cns->haloStressesBytes, NULL);
    cns->o_haloStressesBuffer = mesh->device.malloc(cns->haloStressesBytes);
    cns->sendStressesBuffer = (dfloat*) o_sendBuffer.getMappedPointer();
    cns->recvStressesBuffer = (dfloat*) o_recvBuffer.getMappedPointer();
  }

  //  p_RT, p_rbar, p_ubar, p_vbar
  // p_half, p_two, p_third, p_Nstresses
  
  kernelInfo.addDefine("p_Nfields", mesh->Nfields);
  kernelInfo.addDefine("p_Nstresses", cns->Nstresses);

  kernelInfo.addDefine("p_RT", cns->RT);

  dfloat sqrtRT = sqrt(cns->RT);
  kernelInfo.addDefine("p_sqrtRT", sqrtRT);
  
  kernelInfo.addDefine("p_rbar", cns->rbar);
  kernelInfo.addDefine("p_ubar", cns->ubar);
  kernelInfo.addDefine("p_vbar", cns->vbar);
  

  const dfloat p_one = 1.0, p_two = 2.0, p_half = 1./2., p_third = 1./3., p_zero = 0;

  kernelInfo.addDefine("p_two", p_two);
  kernelInfo.addDefine("p_one", p_one);
  kernelInfo.addDefine("p_half", p_half);
  kernelInfo.addDefine("p_third", p_third);
  kernelInfo.addDefine("p_zero", p_zero);
  
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 512/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  kernelInfo.addDefine("p_Lambda2", 0.5f);

  kernelInfo.addDefine("p_blockSize", blockSize);
  
  cns->volumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsVolumeQuad2D.okl",
				       "cnsVolumeQuad2D",
				       kernelInfo);

  cns->surfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsSurfaceQuad2D.okl",
				       "cnsSurfaceQuad2D",
				       kernelInfo);

  cns->cubatureVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsCubatureVolumeQuad2D.okl",
               "cnsCubatureVolumeQuad2D",
               kernelInfo);

  cns->cubatureSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsCubatureSurfaceQuad2D.okl",
               "cnsCubatureSurfaceQuad2D",
               kernelInfo);

  cns->stressesVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsVolumeQuad2D.okl",
				       "cnsStressesVolumeQuad2D",
				       kernelInfo);

  cns->stressesSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsSurfaceQuad2D.okl",
				       "cnsStressesSurfaceQuad2D",
				       kernelInfo);
  
  cns->vorticityKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsVorticityQuad2D.okl",
				       "cnsVorticityQuad2D",
				       kernelInfo);  
  
  cns->updateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsUpdate2D.okl",
				       "cnsUpdate2D",
				       kernelInfo);

  cns->rkUpdateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsUpdate2D.okl",
				       "cnsRkUpdate2D",
				       kernelInfo);
  
  cns->rkStageKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsUpdate2D.okl",
				       "cnsRkStage2D",
				       kernelInfo);
  
  cns->rkErrorEstimateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsUpdate2D.okl",
				       "cnsErrorEstimate2D",
				       kernelInfo);
   
  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
				       "meshHaloExtract2D",
				       kernelInfo);
  
  return cns;
}
