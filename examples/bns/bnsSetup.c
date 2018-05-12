#include "bns.h"

bns_t *bnsSetup(mesh_t *mesh, setupAide &options){
  
  // OCCA build 
  char deviceConfig[BUFSIZ];
  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  long int hostId = gethostid();

  long int* hostIds = (long int*) calloc(size,sizeof(long int));
  MPI_Allgather(&hostId,1,MPI_LONG,hostIds,1,MPI_LONG,MPI_COMM_WORLD);

  int deviceID = 0;
  for (int r=0;r<rank;r++) {
    if (hostIds[r]==hostId) deviceID++;
  }

  // read thread model/device/platform from options
  if(options.compareArgs("THREAD MODEL", "CUDA")){
    int dev; 
    options.getArgs("DEVICE NUMBER" ,dev);
    sprintf(deviceConfig, "mode = CUDA, deviceID = %d",dev);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenCL")){
    int dev, plat;
    options.getArgs("DEVICE NUMBER", dev);
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "mode = OpenCL, deviceID = %d, platformID = %d", dev, plat);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenMP")){
    sprintf(deviceConfig, "mode = OpenMP");
  }
  else{
    sprintf(deviceConfig, "mode = Serial");
  }
  
  // BNS build
  bns_t *bns = (bns_t*) calloc(1, sizeof(bns_t));

  options.getArgs("MESH DIMENSION", bns->dim);
  options.getArgs("ELEMENT TYPE", bns->elementType);
  

  mesh->Nfields = (bns->dim==3) ? 10:6;
  bns->Nfields = mesh->Nfields; 

  bns->mesh = mesh; 

    
  // Defaulting BNS Input Parameters
  bns->Ma         = 0.1;
  bns->Re         = 100; 
  bns->startTime  = 0.0;
  bns->finalTime  = 1.0;
  bns->cfl        = 0.5; 
  //
  bns->pmlFlag    = 0; 
  bns->probeFlag  = 0; 
  bns->errorFlag  = 0;
  bns->reportFlag = 0;
  bns->errorStep  = 1; 
  bns->reportStep = 1;
  
  // Load input parameters
  int check; 

  check = options.getArgs("REYNOLDS NUMBER", bns->Re);
  if(!check) printf("WARNING setup file does not include REYNOLDS NUMBER\n");

  check = options.getArgs("MACH NUMBER", bns->Ma);
  if(!check) printf("WARNING setup file does not include MACH NUMBER\n");

  check = options.getArgs("START TIME", bns->startTime);
  if(!check) printf("WARNING setup file does not include START TIME\n");

  check = options.getArgs("FINAL TIME", bns->finalTime);
  if(!check) printf("WARNING setup file does not include FINAL TIME\n");

  check = options.getArgs("CFL", bns->cfl);
  if(!check) printf("WARNING setup file does not include CFL\n");

  check = options.getArgs("PROBE FLAG", bns->probeFlag);
  if(!check) printf("WARNING setup file does not include PROBE FLAG\n");

  check = options.getArgs("ERROR FLAG", bns->errorFlag);
  if(!check) printf("WARNING setup file does not include ERROR FLAG\n");

  check = options.getArgs("REPORT FLAG", bns->reportFlag);
  if(!check) printf("WARNING setup file does not include REPORT FLAG\n");

  check = options.getArgs("FIXED TIME STEP", bns->fixed_dt);
  if(!check) printf("WARNING setup file does not include FIXED TIME STEP\n");

  if(options.compareArgs("ABSORBING LAYER", "PML"))
    bns->pmlFlag = 1; 

  if(bns->pmlFlag){
   check = options.getArgs("PML PROFILE ORDER", bns->pmlOrder);
   if(!check) printf("WARNING setup file does not include PML ORDER\n");

   check = options.getArgs("PML SIGMAX MAX", bns->sigmaXmax);
   if(!check) printf("WARNING setup file does not include PML SIGMAX MAX\n");

   check = options.getArgs("PML SIGMAY MAX", bns->sigmaYmax);
   if(!check) printf("WARNING setup file does not include PML SIGMAY MAX\n");

    if(bns->dim==3){
      check = options.getArgs("PML SIGMAZ MAX", bns->sigmaZmax);
      if(!check) printf("WARNING setup file does not include PML SIGMAZ MAX\n");
    }
  }
 
  
  // Set time discretization scheme:fully explicit or not
  bns->fexplicit = 0; 
  if(options.compareArgs("TIME INTEGRATOR", "LSERK") ) // fully explicit schemes
    bns->fexplicit = 1; 
  // Report / error time steps for fixed dt integration
  options.getArgs("TSTEPS FOR SOLUTION OUTPUT", bns->reportStep);
  options.getArgs("TSTEPS FOR ERROR COMPUTE",   bns->errorStep);
  // Output interval for variable dt integration
  options.getArgs("OUTPUT INTERVAL",   bns->outputInterval);

  printf("=============WRITING INPUT PARAMETERS===================\n");

  printf("REYNOLDS NUMBER\t:\t%.2e\n", bns->Re);
  printf("MACH NUMBER\t:\t%.2e\n", bns->Ma);
  printf("CFL NUMBER\t:\t%.2e\n", bns->cfl);
  printf("START TIME\t:\t%.2e\n", bns->startTime);
  printf("FINAL TIME\t:\t%.2e\n", bns->finalTime);
  printf("FIXED DT\t:\t%d\n", bns->fixed_dt);
  printf("PML FORMULATION\t:\t%d\n", bns->pmlFlag);
  if(bns->pmlFlag){
    printf("PML PROFILE N\t:\t%d\n", bns->pmlOrder);
    printf("PML SIGMA X\t:\t%.2e\n", bns->sigmaXmax);
    printf("PML SIGMA Y\t:\t%.2e\n", bns->sigmaYmax);
    printf("PML SIGMA Y\t:\t%.2e\n", bns->sigmaYmax);
  }

    

  // printf("Starting initial conditions\n");
  // Set characteristic length, should be one in a proper problem setting
  dfloat Uref = 1.0;   
  dfloat Lref = 1.0;   
  //
  bns->RT       = Uref*Uref/(bns->Ma*bns->Ma);
  bns->sqrtRT   = sqrt(bns->RT);  
  //
  dfloat nu     = Uref*Lref/bns->Re; 
  bns->tauInv   = bns->RT/nu;

  // Set penalty parameter for flux setting
  bns->Lambda2 = 0.5/(bns->sqrtRT);
  
  // Setting initial conditions
  dfloat rho     = 1.,   u       = 1.,  v       = 0.,   w  = 0.; 
  dfloat sigma11 = 0.f , sigma12 = 0.f, sigma13 = 0.f;
  dfloat sigma22 = 0.f , sigma23 = 0.f;
  dfloat sigma33 = 0.f;

  dfloat q1bar =0., q2bar =0., q3bar =0., q4bar =0., q5bar = 0.; 
  dfloat q6bar =0., q7bar =0., q8bar =0., q9bar =0., q10bar =0.; 
  // Set time step size
  if(bns->dim==2){
     printf("MESH DIMENSION\t:\t%d\n", bns->dim);
    q1bar = rho;
    q2bar = rho*u/bns->sqrtRT;
    q3bar = rho*v/bns->sqrtRT;
    q4bar = (rho*u*v - sigma12)/bns->RT;
    q5bar = (rho*u*u - sigma11)/(sqrt(2.)*bns->RT);
    q6bar = (rho*v*v - sigma22)/(sqrt(2.)*bns->RT);    
  } else{
    printf("MESH DIMENSION\t:\t%d\n", bns->dim);
    
    q1bar  = rho;
    q2bar  = rho*u/bns->sqrtRT;
    q3bar  = rho*v/bns->sqrtRT;
    q4bar  = rho*w/bns->sqrtRT;
    //
    q5bar  = (rho*u*v - sigma12)/bns->RT;
    q6bar  = (rho*u*w - sigma13)/bns->RT;
    q7bar  = (rho*v*w - sigma23)/bns->RT;
    //
    q8bar  = (rho*u*u - sigma11)/(sqrt(2.)*bns->RT);
    q9bar  = (rho*v*v - sigma22)/(sqrt(2.)*bns->RT);
    q10bar = (rho*w*w - sigma33)/(sqrt(2.)*bns->RT);
  }

  dfloat magVelocity = 1.0; 
  if(bns->dim==2)
    magVelocity  = sqrt(q2bar*q2bar+q3bar*q3bar)/(q1bar/bns->sqrtRT);
  else
    magVelocity  = sqrt(q2bar*q2bar+q3bar*q3bar+q4bar*q4bar)/(q1bar/bns->sqrtRT);
  
  // Correction for initial zero velocity
  magVelocity         = mymax(magVelocity,1.0); 

  dfloat ghmin        = 1e9; 
  dfloat dt           = 1e9; 
  dfloat *EtoDT       = (dfloat *) calloc(mesh->Nelements,sizeof(dfloat));

  //Set time step size
  for(dlong e=0;e<mesh->Nelements;++e)
  { 
    dfloat hmin = 1e9, dtmax = 1e9;
    
    EtoDT[e] = dtmax;

    for(int f=0;f<mesh->Nfaces;++f){
      dlong sid   = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];
     
      dfloat htest = 2.0/(sJ*invJ);
      hmin = mymin(hmin, htest); 
    }

    ghmin = mymin(ghmin, hmin);

    dfloat dtex   = bns->cfl*hmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*bns->sqrtRT);
    dfloat dtim   = bns->cfl*1.f/(bns->tauInv);
    dfloat dtest  = 1e9;
    
    if(bns->fexplicit) // fully explicit schemes
      dtest  = mymin(dtex,dtim); 
    else              // implicit-explicit or semi-analytic schemes
      dtest  = dtex;
      
      dt            = mymin(dt, dtest);      // For SR 
      EtoDT[e]      = mymin(EtoDT[e],dtest); // For MR
  }


  


   // printf("dtex = %.5e dtim = %.5e \n", bns->cfl*ghmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*bns->sqrtRT), bns->cfl*1.f/(bns->tauInv));

  // Set multiRate element groups/group  
  if(options.compareArgs("TIME INTEGRATOR", "MRSAAB") ){
    int maxLevels =0; options.getArgs("MAX MRAB LEVELS", maxLevels);

    printf("MR MAX LEVELS\t:\t%d\n", maxLevels);

    if(bns->dim==2)
      bns->dt = meshMRABSetup2D(mesh,EtoDT,maxLevels, (bns->finalTime-bns->startTime));
    else
      bns->dt = meshMRABSetup3D(mesh,EtoDT,maxLevels, (bns->finalTime-bns->startTime));

    bns->NtimeSteps =  (bns->finalTime-bns->startTime)/(pow(2,mesh->MRABNlevels-1)*bns->dt);

    printf("MR  LEVELS\t:\t%d\n", mesh->MRABNlevels);
  }
  else{
    // printf("MESH DIMENSION\t:\t%d\n", bns->dt);  
    //!!!!!!!!!!!!!! Fix time step to compute the error in postprecessing step  
    // MPI_Allreduce to get global minimum dt
    MPI_Allreduce(&dt, &(bns->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
    bns->NtimeSteps = (bns->finalTime-bns->startTime)/bns->dt;
    bns->dt         = (bns->finalTime-bns->startTime)/bns->NtimeSteps;
    //offset index
    bns->shiftIndex = 0;
    // Set element ids for nonPml region, will be modified if PML exists
    mesh->nonPmlNelements = mesh->Nelements; 
    mesh->pmlNelements    = 0; 

    mesh->nonPmlElementIds = (dlong*) calloc(mesh->Nelements, sizeof(dlong));
    for(dlong e=0;e<mesh->Nelements;++e)
     mesh->nonPmlElementIds[e] = e; 

   printf("TIME STEPSIZE\t:\t%.2e\n", bns->dt); 
  }



  if(options.compareArgs("TIME INTEGRATOR", "MRSAAB")){
    bns->Nrhs = 3; 
    // compute samples of q at interpolation nodes
    bns->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rhsq = (dfloat*) calloc(bns->Nrhs*mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->fQM  = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*bns->Nfields, sizeof(dfloat));
  }



  // Initialize  
  if (options.compareArgs("TIME INTEGRATOR","LSERK")){ // LSERK,
    bns->Nrhs = 1; 
    // compute samples of q at interpolation nodes
    bns->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rhsq = (dfloat*) calloc(bns->Nrhs*mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->resq = (dfloat*) calloc(bns->Nrhs*mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
  }

   // Initialize  
  if (options.compareArgs("TIME INTEGRATOR","SARK")){ // SARK for fixed or adaptive time stepping,
    bns->Nrhs  = 1;
    bns->ATOL    = 1.0; options.getArgs("ABSOLUTE TOLERANCE",   bns->ATOL); 
    bns->RTOL    = 1.0; options.getArgs("RELATIVE TOLERANCE",   bns->RTOL);
    bns->dtMIN   = 1.0; options.getArgs("MINUMUM TIME STEP SIZE",   bns->dtMIN); 
    bns->emethod = 0; // 0 PID / 1 PI / 2 P / 3 I    
    bns->rkp     = 4; // order of embedded scheme + 1 

    
    dlong Ntotal  = mesh->Nelements*mesh->Np*bns->Nfields;
    bns->Nblock   = (Ntotal+blockSize-1)/blockSize;

    dlong localElements =  mesh->Nelements;
    MPI_Allreduce(&localElements, &(bns->totalElements), 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    // compute samples of q at interpolation nodes
    bns->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    //
    bns->NrkStages = 7; // 7 stages order(54) or 6 stages order(43) 
    
    bns->rkq      = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rkrhsq   = (dfloat*) calloc(bns->NrkStages*mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rkerr    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->errtmp  =  (dfloat*) calloc(bns->Nblock, sizeof(dfloat)); 
  }

 

  dfloat time = bns->startTime + 0.0; 
  // Define Initial Mean Velocity
  dfloat ramp, drampdt;
  bnsRampFunction(time, &ramp, &drampdt);

 // INITIALIZE PROBLEM 
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat t = 0., x = 0., y = 0., z = 0.;
      x = mesh->x[n + mesh->Np*e];
      y = mesh->y[n + mesh->Np*e];
      if(bns->dim==3)
        z = mesh->z[n + mesh->Np*e];

      const dlong id = e*bns->Nfields*mesh->Np +n; 
      if(bns->dim==2){
        // Uniform Flow
        bns->q[id+0*mesh->Np] = q1bar; 
        bns->q[id+1*mesh->Np] = ramp*q2bar;
        bns->q[id+2*mesh->Np] = ramp*q3bar;
        bns->q[id+3*mesh->Np] = ramp*ramp*q4bar;
        bns->q[id+4*mesh->Np] = ramp*ramp*q5bar;
        bns->q[id+5*mesh->Np] = ramp*ramp*q6bar;   
      }else{

        #if 0
          dfloat Pmax = 0.5; 
          dfloat U0   = 1.0; 
          dfloat r0   = 1.0;
          dfloat gamma= 1.4; 

          dfloat rho = 1.0 + Pmax*exp(-(log(2)* (x*x + y*y + z*z)/r0));

          bns->q[id+0*mesh->Np] = rho; 
          bns->q[id+1*mesh->Np] = rho*U0/bns->sqrtRT;
          bns->q[id+2*mesh->Np] = 0.0;
          bns->q[id+3*mesh->Np] = 0.0;

          bns->q[id+4*mesh->Np] = ramp*ramp*q5bar;
          bns->q[id+5*mesh->Np] = ramp*ramp*q6bar;
          bns->q[id+6*mesh->Np] = ramp*ramp*q7bar;   
          bns->q[id+7*mesh->Np] = ramp*ramp*q8bar;   
          bns->q[id+8*mesh->Np] = ramp*ramp*q9bar;   
          bns->q[id+9*mesh->Np] = ramp*ramp*q10bar; 
        #else
          bns->q[id+0*mesh->Np] = q1bar; 
          bns->q[id+1*mesh->Np] = ramp*q2bar;
          bns->q[id+2*mesh->Np] = ramp*q3bar;
          bns->q[id+3*mesh->Np] = ramp*q4bar;

          bns->q[id+4*mesh->Np] = ramp*ramp*q5bar;
          bns->q[id+5*mesh->Np] = ramp*ramp*q6bar;
          bns->q[id+6*mesh->Np] = ramp*ramp*q7bar;   
          bns->q[id+7*mesh->Np] = ramp*ramp*q8bar;   
          bns->q[id+8*mesh->Np] = ramp*ramp*q9bar;   
          bns->q[id+9*mesh->Np] = ramp*ramp*q10bar; 
        #endif  
      }
       
    }
  }

  
  // Write Problem Info 
  if(rank==0){
    printf("dt   = %g max wave speed = %g",   bns->dt, sqrt(3.)*bns->sqrtRT);
    if(mesh->MRABNlevels)
      printf("Nsteps = %d dt = %.8e MRAB Level: %d  Final Time:%.5e\n", bns->NtimeSteps, bns->dt, mesh->MRABNlevels, bns->startTime+pow(2, mesh->MRABNlevels-1)*(bns->dt*(bns->NtimeSteps+1)));   
    else
     printf("Nsteps = %d dt = %.8e Final Time:%.5e\n", bns->NtimeSteps, bns->dt,  bns->startTime + bns->dt*bns->NtimeSteps);
  }
 
 
  // SET PROBE DATA
  if(bns->probeFlag){
    // mesh->probeNTotal = 3; 
    // dfloat *pX   = (dfloat *) calloc (mesh->probeNTotal, sizeof(dfloat));
    // dfloat *pY   = (dfloat *) calloc (mesh->probeNTotal, sizeof(dfloat));
    // // Fill probe coordinates
    //  pX[0] = 9.00;  pX[1] = 10.00;  pX[2] = 10.50; //pX[3] = 5;
    //  pY[0] = 5.00;  pY[1] =  6.00;  pY[2] =  6.50; //pY[3] = 5*tan(M_PI/6.);
    // meshProbeSetup2D(mesh, pX, pY);
    // free(pX); free(pY);

  }
  
   

  occa::kernelInfo kernelInfo;
  if(bns->dim==3)
    meshOccaSetup3D(mesh, deviceConfig, kernelInfo);
  else
    meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

  kernelInfo.addParserFlag("automate-add-barriers", "disabled");   

 

  // Setup MRAB PML
  if(options.compareArgs("TIME INTEGRATOR","MRSAAB")){
    printf("Preparing Pml for multirate rate\n");
    
    bnsMRABPmlSetup(bns, options);

    mesh->o_MRABelementIds = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    mesh->o_MRABhaloIds    = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    for (int lev=0;lev<mesh->MRABNlevels;lev++) {
      if (mesh->MRABNelements[lev])
        mesh->o_MRABelementIds[lev] = mesh->device.malloc(mesh->MRABNelements[lev]*sizeof(dlong),mesh->MRABelementIds[lev]);
      if (mesh->MRABNhaloElements[lev])
        mesh->o_MRABhaloIds[lev] = mesh->device.malloc(mesh->MRABNhaloElements[lev]*sizeof(dlong), mesh->MRABhaloIds[lev]);
    }
  } 
  else{
   printf("Preparing Pml for single rate integrator\n");
   bnsPmlSetup(bns, options); 

    if (mesh->nonPmlNelements)
        mesh->o_nonPmlElementIds = mesh->device.malloc(mesh->nonPmlNelements*sizeof(dlong), mesh->nonPmlElementIds);

  }
  





// Compute Time Stepper Coefficcients

bnsTimeStepperCoefficients(bns, options);


if(options.compareArgs("TIME INTEGRATOR","MRSAAB")){

  bns->o_q     = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat),bns->q);
  bns->o_rhsq  = mesh->device.malloc(bns->Nrhs*mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->rhsq);
  
  //reallocate halo buffer for trace exchange
  if (mesh->totalHaloPairs) {
    mesh->o_haloBuffer.free();
    mesh->o_haloBuffer = mesh->device.malloc(mesh->totalHaloPairs*mesh->Nfp*bns->Nfields*mesh->Nfaces*sizeof(dfloat));
  }

  bns->o_fQM = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*bns->Nfields*sizeof(dfloat),
                          bns->fQM);
  mesh->o_mapP = mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(int), mesh->mapP);
}


if(options.compareArgs("TIME INTEGRATOR", "LSERK")){
  // 
  bns->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->q);
  bns->o_rhsq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->rhsq);
  bns->o_resq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->resq);
}

if(options.compareArgs("TIME INTEGRATOR","SARK")){

  bns->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->q);
  bns->o_rhsq = 
    mesh->device.malloc(bns->Nrhs*mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->rhsq); 
  
  int Ntotal    = mesh->Nelements*mesh->Np*bns->Nfields;
  
  bns->o_rkq =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->rkq);

  bns->o_saveq =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->rkq);

  bns->o_rkrhsq =
    mesh->device.malloc(bns->NrkStages*mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->rkrhsq);
  bns->o_rkerr =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->rkerr);
  bns->o_errtmp = mesh->device.malloc(bns->Nblock*sizeof(dfloat), bns->errtmp);

  bns->o_rkA = mesh->device.malloc(bns->NrkStages*bns->NrkStages*sizeof(dfloat), bns->rkA);
  bns->o_rkE = mesh->device.malloc(bns->NrkStages*sizeof(dfloat), bns->rkE);

}






  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  int maxCubNodes = mymax(maxNodes,mesh->cubNp);

  kernelInfo.addDefine("p_maxNodes", maxNodes);
  kernelInfo.addDefine("p_maxCubNodes", maxCubNodes);


  int NblockV = 128/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 128/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  int NblockCub = 128/mesh->cubNp; // works for CUDA

  NblockCub = 2; // !!!!!!!!!!!!!!!!!!!!!

  kernelInfo.addDefine("p_NblockCub", NblockCub);

  // physics 
  kernelInfo.addDefine("p_Lambda2", 0.5f);
  kernelInfo.addDefine("p_sqrtRT", bns->sqrtRT);
  kernelInfo.addDefine("p_sqrt2", (dfloat)sqrt(2.));
  kernelInfo.addDefine("p_isq12", (dfloat)sqrt(1./12.));
  kernelInfo.addDefine("p_isq6", (dfloat)sqrt(1./6.));
  kernelInfo.addDefine("p_invsqrt2", (dfloat)sqrt(1./2.));
  kernelInfo.addDefine("p_tauInv", bns->tauInv);

  kernelInfo.addDefine("p_q1bar", q1bar);
  kernelInfo.addDefine("p_q2bar", q2bar);
  kernelInfo.addDefine("p_q3bar", q3bar);
  kernelInfo.addDefine("p_q4bar", q4bar);
  kernelInfo.addDefine("p_q5bar", q5bar);
  kernelInfo.addDefine("p_q6bar", q6bar);
  
  if(bns->dim==3){
    kernelInfo.addDefine("p_q7bar", q7bar);
    kernelInfo.addDefine("p_q8bar", q8bar);
    kernelInfo.addDefine("p_q9bar", q9bar);
    kernelInfo.addDefine("p_q10bar", q10bar);
  }

  kernelInfo.addDefine("p_alpha0", (dfloat).01f);
  kernelInfo.addDefine("p_pmlAlpha", (dfloat)0.1f);
  kernelInfo.addDefine("p_blockSize", blockSize);
  kernelInfo.addDefine("p_NrkStages", bns->NrkStages);


  if(bns->fexplicit) // full explicit or semi-nalaytic
    kernelInfo.addDefine("p_SEMI_ANALYTIC", (int) 0);
  else
    kernelInfo.addDefine("p_SEMI_ANALYTIC", (int) 1);

  if(options.compareArgs("TIME INTEGRATOR", "MRSAAB"))
    kernelInfo.addDefine("p_MRSAAB", (int) 1);
  else
    kernelInfo.addDefine("p_MRSAAB", (int) 0);






  // set kernel name suffix
  char *suffix;
  
  if(bns->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(bns->elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(bns->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(bns->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

char fileName[BUFSIZ], kernelName[BUFSIZ];
for (int r=0;r<size;r++){

    if (r==rank) {
        
     
        // Volume kernels
        sprintf(fileName, "okl/bnsVolume%s.okl", suffix);

        sprintf(kernelName, "bnsVolume%s", suffix);
        bns->volumeKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);
        sprintf(kernelName, "bnsPmlVolume%s", suffix);
        bns->pmlVolumeKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);
       
        // Relaxation kernels
        sprintf(fileName, "okl/bnsRelaxation%s.okl", suffix);
  
        sprintf(kernelName, "bnsRelaxation%s", suffix);
        bns->relaxationKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);
        sprintf(kernelName, "bnsPmlRelaxation%s", suffix);        
        bns->pmlRelaxationKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);
        
        // Surface kernels 
        sprintf(fileName, "okl/bnsSurface%s.okl", suffix);

        if(options.compareArgs("TIME INTEGRATOR","MRSAAB")){
          sprintf(kernelName, "bnsMRSurface%s", suffix);
          bns->surfaceKernel = mesh->device.buildKernelFromSource(fileName,kernelName, kernelInfo);

          sprintf(kernelName, "bnsMRPmlSurface%s", suffix);
          bns->pmlSurfaceKernel = mesh->device.buildKernelFromSource(fileName,kernelName, kernelInfo);
        }else{
          sprintf(kernelName, "bnsSurface%s", suffix);
          bns->surfaceKernel = mesh->device.buildKernelFromSource(fileName,kernelName, kernelInfo);

          sprintf(kernelName, "bnsPmlSurface%s", suffix);
          bns->pmlSurfaceKernel = mesh->device.buildKernelFromSource(fileName,kernelName, kernelInfo);
        }

        
        sprintf(fileName, "okl/bnsUpdate%s.okl", suffix);
        // Update Kernels
        if(options.compareArgs("TIME INTEGRATOR","LSERK")){
          sprintf(kernelName, "bnsLSERKUpdate%s", suffix);
          bns->updateKernel = mesh->device.buildKernelFromSource(fileName, kernelName,kernelInfo);

          sprintf(kernelName, "bnsLSERKPmlUpdate%s", suffix);
          bns->pmlUpdateKernel = mesh->device.buildKernelFromSource(fileName, kernelName,kernelInfo);
        }
        else if(options.compareArgs("TIME INTEGRATOR","SARK")){
          sprintf(kernelName, "bnsSARKUpdateStage%s", suffix);
          bns->updateStageKernel = mesh->device.buildKernelFromSource(fileName,kernelName, kernelInfo);

          sprintf(kernelName, "bnsSARKPmlUpdateStage%s", suffix);
          bns->pmlUpdateStageKernel = mesh->device.buildKernelFromSource(fileName,kernelName, kernelInfo);

          sprintf(kernelName, "bnsSARKUpdate%s", suffix);
          bns->updateKernel = mesh->device.buildKernelFromSource(fileName, kernelName,kernelInfo);

          sprintf(kernelName, "bnsSARKPmlUpdate%s", suffix);
          bns->pmlUpdateKernel = mesh->device.buildKernelFromSource(fileName, kernelName,kernelInfo);
        }
        else if(options.compareArgs("TIME INTEGRATOR","MRSAAB")){
        
        sprintf(kernelName, "bnsMRSAABTraceUpdate%s", suffix);
        bns->traceUpdateKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "bnsMRSAABUpdate%s", suffix);
        bns->updateKernel = mesh->device.buildKernelFromSource(fileName, kernelName,kernelInfo);

        sprintf(kernelName, "bnsMRSAABPmlUpdate%s", suffix);
        bns->pmlUpdateKernel = mesh->device.buildKernelFromSource(fileName, kernelName,kernelInfo);






        }


    }
  MPI_Barrier(MPI_COMM_WORLD);
}



#if 0








 // Volume and Relaxation Kernels
  if(options.compareArgs("RELAXATION TYPE","CUBATURE") && 
     !(options.compareArgs("TIME INTEGRATOR","LSIMEX") || options.compareArgs("TIME INTEGRATOR","IMEXRK"))){ 
          
    printf("Compiling volume kernel for cubature integration\n");
      bns->volumeKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
         "boltzmannVolumeCub2D",
          kernelInfo);

    printf("Compiling PML volume kernel for cubature integration\n");
      bns->pmlVolumeKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
        "boltzmannPmlVolumeCub2D",
          kernelInfo);

    if(options.compareArgs("TIME INTEGRATOR","MRAB") || options.compareArgs("TIME INTEGRATOR","LSERK") || 
       options.compareArgs("TIME INTEGRATOR","SRAB") || options.compareArgs("TIME INTEGRATOR","DOPRI5") ){ 

      printf("Compiling relaxation kernel with cubature integration\n");
       bns->relaxationKernel =
         mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannRelaxation2D.okl",
          "boltzmannRelaxationCub2D",
           kernelInfo); 
    
      //
      printf("Compiling PML relaxation kernel with cubature integration\n");
        bns->pmlRelaxationKernel =
         mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannRelaxation2D.okl",
          "boltzmannPmlRelaxationCub2D",
            kernelInfo);  
    }
    else if(options.compareArgs("TIME INTEGRATOR","MRSAAB")||options.compareArgs("TIME INTEGRATOR","SAAB") || 
            options.compareArgs("TIME INTEGRATOR","SARK")   || options.compareArgs("TIME INTEGRATOR","XDOPRI") ||
            options.compareArgs("TIME INTEGRATOR","SAADRK")){

    printf("Compiling relaxation kernel with cubature integration\n");
     bns->relaxationKernel =
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannRelaxation2D.okl",
        "boltzmannSARelaxationCub2D",
        kernelInfo); 

      //
    printf("Compiling PML relaxation kernel with cubature integration\n");
    bns->pmlRelaxationKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannRelaxation2D.okl",
      "boltzmannSAPmlRelaxationCub2D",
       kernelInfo);  
    }

  }



  if(options.compareArgs("RELAXATION TYPE","COLLOCATION")  && 
    !(options.compareArgs("TIME INTEGRATOR","LSIMEX") || options.compareArgs("TIME INTEGRATOR","IMEXRK"))){ 

     if(options.compareArgs("TIME INTEGRATOR", "MRAB") || options.compareArgs("TIME INTEGRATOR","LSERK") || 
        options.compareArgs("TIME INTEGRATOR","SRAB") || options.compareArgs("TIME INTEGRATOR","DOPRI5")){ 

      printf("Compiling volume kernel with nodal collocation for nonlinear term\n");
      bns->volumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
             "boltzmannVolume2D",
             kernelInfo);

      printf("Compiling PML volume kernel with nodal collocation for nonlinear term\n");
      bns->pmlVolumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
           "boltzmannPmlVolume2D",
           kernelInfo); 
    }

    else if(options.compareArgs("TIME INTEGRATOR","MRSAAB") || options.compareArgs("TIME INTEGRATOR","SAAB") || 
            options.compareArgs("TIME INTEGRATOR","SARK") || options.compareArgs("TIME INTEGRATOR","XDOPRI") || 
            options.compareArgs("TIME INTEGRATOR","SAADRK")){
      printf("Compiling SA volume kernel with nodal collocation for nonlinear term\n");
      bns->volumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
             "boltzmannSAVolume2D",
             kernelInfo);

      printf("Compiling SA PML volume kernel with nodal collocation for nonlinear term\n");
      bns->pmlVolumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
           "boltzmannSAPmlVolume2D",
           kernelInfo); 
      }


  }

 

   // UPDATE Kernels
  if(options.compareArgs("TIME INTEGRATOR","MRAB")){ 
    printf("Compiling MRAB update kernel\n");
    bns->updateKernel =
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannMRABUpdate2D",
            kernelInfo);
    
    printf("Compiling MRAB trace update kernel\n");
    bns->traceUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
               "boltzmannMRABTraceUpdate2D",
               kernelInfo);
    
    printf("Compiling MRAB PML update kernel\n");
    bns->pmlUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannMRABPmlUpdate2D",
            kernelInfo);
   
    printf("Compiling MRAB PML trace update kernel\n");
    bns->pmlTraceUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
               "boltzmannMRABPmlTraceUpdate2D",
                 kernelInfo);
    }

    else if(options.compareArgs("TIME INTEGRATOR","MRSAAB")){

    printf("Compiling MRSAAB update kernel\n");
    bns->updateKernel =
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannMRSAABUpdate2D",
            kernelInfo);

    printf("Compiling MRSAAB trace update kernel\n");
    bns->traceUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
               "boltzmannMRSAABTraceUpdate2D",
               kernelInfo);
    
     printf("Compiling MRSAAB PML update kernel\n");
    bns->pmlUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannMRSAABPmlUpdate2D",
            kernelInfo);
    
    printf("Compiling MRSAAB PML trace update kernel\n");
    bns->pmlTraceUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
               "boltzmannMRSAABPmlTraceUpdate2D",
                 kernelInfo);
     }

     else if(options.compareArgs("TIME INTEGRATOR","SRAB")){
     printf("Compiling SRAB update kernel\n");
     bns->updateKernel =
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannSRABUpdate2D",
            kernelInfo);

      printf("Compiling SRAB PML update kernel\n");
      bns->pmlUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannSRABPmlUpdate2D",
            kernelInfo);

     //   printf("Compiling LSERK update kernel\n");
     // bns->RKupdateKernel =
     // mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
     //     "boltzmannLSERKUpdate2D",
     //        kernelInfo);

     //  printf("Compiling LSERK PML update kernel\n");
     //  bns->RKpmlUpdateKernel =
     //  mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
     //     "boltzmannLSERKPmlUpdate2D",
     //        kernelInfo);


     }


     else if(options.compareArgs("TIME INTEGRATOR","SRSAAB")){
     printf("Compiling SAAB update kernel\n");
     bns->updateKernel =
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannSAABUpdate2D",
            kernelInfo);

      printf("Compiling SAAB PML update kernel\n");
      bns->pmlUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannSAABPmlUpdate2D",
            kernelInfo);
     }

     else if(options.compareArgs("TIME INTEGRATOR","LSERK")){
     printf("Compiling LSERK update kernel\n");
     bns->updateKernel =
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannLSERKUpdate2D",
            kernelInfo);

      printf("Compiling LSERK PML update kernel\n");
      bns->pmlUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannLSERKPmlUpdate2D",
            kernelInfo);
     }



    else if(options.compareArgs("TIME INTEGRATOR","SARK")){
    //SARK STAGE UPDATE
    printf("compiling SARK non-pml  update kernel\n");
    bns->updateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannSARKUpdate2D",
         kernelInfo); 


    printf("compiling SARK non-pml  stage update kernel\n");
    mesh->updateStageKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannSARKStageUpdate2D",
         kernelInfo); 

    //
    printf("compiling SARK Unsplit pml  update kernel\n");
    bns->pmlUpdateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
       "boltzmannSARKPmlUpdate2D",
       kernelInfo); 

    printf("compiling SARK Unsplit pml stage update kernel\n");
    bns->pmlUpdateStageKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
    "boltzmannSARKPmlStageUpdate2D",
    kernelInfo); 
    }

    // DOPRI5 Update Kernels
    else if(options.compareArgs("TIME INTEGRATOR","DOPRI5")){

    bns->updateStageKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannDOPRIRKStage2D",
          kernelInfo); 

    bns->pmlUpdateStageKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannDOPRIRKPmlStage2D",
          kernelInfo); 


    bns->updateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannDOPRIRKUpdate2D",
          kernelInfo); 

    bns->pmlUpdateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannDOPRIRKPmlUpdate2D",
          kernelInfo); 


    bns->errorEstimateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannErrorEstimate2D",
          kernelInfo); 


    }


    else if(options.compareArgs("TIME INTEGRATOR","XDOPRI")){
    // printf("compiling XDOPRI Update kernels\n");
    bns->updateStageKernel =
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
      "boltzmannXDOPRIRKStage2D",
        kernelInfo); 

    bns->pmlUpdateStageKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannXDOPRIRKPmlStage2D",
          kernelInfo);

    bns->updateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannXDOPRIRKUpdate2D",
          kernelInfo); 

    bns->pmlUpdateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannXDOPRIRKPmlUpdate2D",
          kernelInfo); 


    bns->errorEstimateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannErrorEstimate2D",
          kernelInfo);  

     }


     else if(options.compareArgs("TIME INTEGRATOR","SAADRK")){

      printf("compiling SAADRK Update kernels\n");
      bns->updateStageKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannSAADRKStage2D",
          kernelInfo); 

       bns->pmlUpdateStageKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannSAADRKPmlStage2D",
          kernelInfo);

       bns->updateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannSAADRKUpdate2D",
          kernelInfo); 

       bns->pmlUpdateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannSAADRKPmlUpdate2D",
          kernelInfo); 

     
       bns->errorEstimateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannErrorEstimate2D",
          kernelInfo); 





     }

    



    // Surface Kernels
    if(options.compareArgs("TIME INTEGRATOR","MRAB") || options.compareArgs("TIME INTEGRATOR","MRSAAB")){  
    printf("Compiling surface kernel\n");
    bns->surfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
         "boltzmannMRSurface2D",
         kernelInfo);

    printf("Compiling PML surface kernel\n");
    bns->pmlSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
       "boltzmannMRPmlSurface2D",
       kernelInfo);
    }

    else { 
    printf("Compiling surface kernel\n");
    bns->surfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
         "boltzmannSurface2D",
         kernelInfo);

    printf("Compiling PML surface kernel\n");
    bns->pmlSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
       "boltzmannPmlSurface2D",
       kernelInfo);
    }







    // IMEX Kernels
    if(options.compareArgs("TIME INTEGRATOR","LSIMEX")){ 
    // RESIDUAL UPDATE KERNELS
    printf("Compiling LSIMEX non-pml residual update kernel\n");
    bns->residualUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
           "boltzmannLSIMEXResidualUpdate2D",
           kernelInfo);
    
    printf("Compiling LSIMEX non-pml implicit update kernel\n");
    bns->implicitUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
           "boltzmannLSIMEXImplicitUpdate2D",
           kernelInfo);

    printf("Compiling LSIMEX Unsplit pml residual update kernel\n");
   bns->pmlResidualUpdateKernel =
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
          "boltzmannLSIMEXPmlResidualUpdate2D",
          kernelInfo);
     //
   printf("Compiling LSIMEX Unsplit pml implicit update kernel\n");
   bns->pmlImplicitUpdateKernel =
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
          "boltzmannLSIMEXPmlImplicitUpdate2D",
          kernelInfo);
      
     
    if(options.compareArgs("RELAXATION TYPE","CUBATURE")){ 
      printf("Compiling LSIMEX non-pml Implicit Iteration Cubature  kernel\n");

      bns->implicitSolveKernel = 
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannLSIMEXImplicitSolve2D.okl",
             "boltzmannLSIMEXImplicitSolveCub2D",
             kernelInfo); 

          printf("Compiling LSIMEX pml Implicit Iteration Cubature  kernel\n");
      bns->pmlImplicitSolveKernel = 
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannLSIMEXImplicitSolve2D.okl",
             "boltzmannLSIMEXPmlImplicitSolveCub2D",
             kernelInfo); 
    }
    else if(options.compareArgs("RELAXATION TYPE","COLLOCATION")){ 
      //
      printf("Compiling LSIMEX non-pml Implicit Iteration kernel\n");
      bns->implicitSolveKernel = 
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannLSIMEXImplicitSolve2D.okl",
        "boltzmannLSIMEXImplicitSolve2D",
        kernelInfo); 
        
     printf("Compiling LSIMEX Unsplit pml Implicit Iteration  kernel\n");
     bns->pmlImplicitSolveKernel = 
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannLSIMEXImplicitSolve2D.okl",
     "boltzmannLSIMEXImplicitSolve2D",
     kernelInfo);      
      
    }

  printf("Compiling LSIMEX volume kernel integration\n");
    bns->volumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
         "boltzmannVolumeCub2D",
         kernelInfo);

       
    printf("Compiling LSERK Unsplit pml volume kernel with cubature integration\n");
    bns->pmlVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
        "boltzmannIMEXPmlVolumeCub2D",
        kernelInfo);
       
   printf("Compiling surface kernel\n");
    bns->surfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
         "boltzmannSurface2D",
         kernelInfo);
        
    printf("Compiling Unsplit  pml surface kernel\n");
    bns->pmlSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
       "boltzmannPmlSurface2D",
       kernelInfo);

    printf("Compiling LSIMEX non-pml update kernel\n");
    bns->updateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
      "boltzmannLSIMEXUpdate2D",
      kernelInfo);
  //
    printf("Compiling LSIMEX Unsplit pml update kernel\n");
       bns->pmlUpdateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
       "boltzmannLSIMEXPmlUpdate2D",
       kernelInfo);
  }


  if(options.compareArgs("TIME INTEGRATOR","IMEXRK")){

    printf("compiling IMEXRK Update kernels\n");
      bns->updateStageKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannIMEXRK2D.okl",
        "boltzmannIMEXRKStageUpdate2D",
          kernelInfo); 

    printf("compiling IMEXRK Pml Update kernels\n");
     bns->pmlUpdateStageKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannIMEXRK2D.okl",
        "boltzmannIMEXRKPmlStageUpdate2D",
          kernelInfo); 


     printf("compiling IMEXRK Implicit Solve\n");
     bns->implicitSolveKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannIMEXRK2D.okl",
        "boltzmannIMEXRKImplicitSolveCub2D",
          kernelInfo); 


    printf("Compiling volume kernel for cubature integration\n");
      bns->volumeKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
         "boltzmannVolumeCub2D",
          kernelInfo);

    printf("Compiling PML volume kernel for cubature integration\n");
      bns->pmlVolumeKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
        "boltzmannPmlVolumeCub2D",
          kernelInfo);

    printf("Compiling IMEXRK update kernel\n");
     bns->updateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannIMEXRK2D.okl",
        "boltzmannIMEXRKUpdate2D",
          kernelInfo); 
    
    printf("Compiling IMEXRK Pml update kernel\n");
       bns->pmlUpdateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannIMEXRK2D.okl",
        "boltzmannIMEXRKPmlUpdate2D",
          kernelInfo); 


    printf("Compiling IMEXRK Pml damping kernel\n");
       bns->pmlUpdateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannIMEXRK2D.okl",
        "boltzmannIMEXRKPmlUpdate2D",
          kernelInfo); 


     printf("Compiling IMEXRK error kernel\n");
       bns->pmlDampingKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannIMEXRK2D.okl",
        "boltzmannIMEXRKPmlDampingCub2D",
          kernelInfo); 


       bns->errorEstimateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
        "boltzmannErrorEstimate2D",
          kernelInfo); 

  }

    mesh->haloExtractKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
          "meshHaloExtract2D",
             kernelInfo);




#endif


return bns; 
  
}




