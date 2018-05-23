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

  if (size==1) options.getArgs("DEVICE NUMBER" ,deviceID);

  // read thread model/device/platform from options
  if(options.compareArgs("THREAD MODEL", "CUDA")){
    sprintf(deviceConfig, "mode = CUDA, deviceID = %d",deviceID);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenCL")){
    int plat;
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "mode = OpenCL, deviceID = %d, platformID = %d", deviceID, plat);
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

  check = options.getArgs("VISCOSITY", bns->nu);
  if(!check) printf("WARNING setup file does not include VISCOSITY\n");

  check = options.getArgs("SPEED OF SOUND", bns->sqrtRT);
  if(!check) printf("WARNING setup file does not include SPEED OF SOUND\n");

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


  if(rank==0){
    printf("=============WRITING INPUT PARAMETERS===================\n");

    printf("VISCOSITY\t:\t%.2e\n", bns->nu);
    printf("SPEED OF SOUND\t:\t%.2e\n", bns->sqrtRT);
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
    printf("ERROR STEP\t:\t%d\n", bns->errorStep);
  }
    

  
  // bns->RT       = Uref*Uref/(bns->Ma*bns->Ma);
  bns->RT     = bns->sqrtRT*bns->sqrtRT;  
  bns->tauInv = bns->RT/bns->nu;
  bns->Re = 1.0; 
  bns->Ma = 0.2; 

  // Tentative, depends on the reference velocity;
  // bns->Re = Uref*Lref/bns->nu; 
  // bns->Ma = Uref/bns->RT; 

  // Set penalty parameter for flux setting
  bns->Lambda2 = 0.5/(bns->sqrtRT);
  
  // Setting initial conditions
  dfloat rho     = 1.,   u       = 1.,  v       = 0.,   w  = 0.;

  check = options.getArgs("RBAR", rho);
  if(!check) printf("WARNING setup file does not include RBAR\n");
  
  check = options.getArgs("UBAR", u);
  if(!check) printf("WARNING setup file does not include UBAR\n");

  check = options.getArgs("VBAR", v);
  if(!check) printf("WARNING setup file does not include VBAR\n");

  if(bns->dim==3){
    check = options.getArgs("WBAR", w);
    if(!check) printf("WARNING setup file does not include WBAR\n");
  }

  dfloat sigma11 = 0.f , sigma12 = 0.f, sigma13 = 0.f;
  dfloat sigma22 = 0.f , sigma23 = 0.f;
  dfloat sigma33 = 0.f;

  dfloat q1bar =0., q2bar =0., q3bar =0., q4bar =0., q5bar = 0.; 
  dfloat q6bar =0., q7bar =0., q8bar =0., q9bar =0., q10bar =0.; 
  // Set time step size
  if(bns->dim==2){
    if (rank==0) printf("MESH DIMENSION\t:\t%d\n", bns->dim);
    q1bar = rho;
    q2bar = rho*u/bns->sqrtRT;
    q3bar = rho*v/bns->sqrtRT;
    q4bar = (rho*u*v - sigma12)/bns->RT;
    q5bar = (rho*u*u - sigma11)/(sqrt(2.)*bns->RT);
    q6bar = (rho*v*v - sigma22)/(sqrt(2.)*bns->RT);    
  } else{
    if (rank==0) printf("MESH DIMENSION\t:\t%d\n", bns->dim);
    
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
  for(dlong e=0;e<mesh->Nelements;++e){ 
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

    if (rank==0) printf("MR MAX LEVELS\t:\t%d\n", maxLevels);

    if(bns->dim==2)
      bns->dt = meshMRABSetup2D(mesh,EtoDT,maxLevels, (bns->finalTime-bns->startTime));
    else
      bns->dt = meshMRABSetup3D(mesh,EtoDT,maxLevels, (bns->finalTime-bns->startTime));

    bns->NtimeSteps =  (bns->finalTime-bns->startTime)/(pow(2,mesh->MRABNlevels-1)*bns->dt);

    if (rank==0) printf("MR  LEVELS\t:\t%d\n", mesh->MRABNlevels);
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
    bns->rkp     = 5; // order of embedded scheme + 1 

    
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

 
  dfloat Gy = 0.0; options.getArgs("BODYFORCE-Y",Gy);
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
        // dfloat phi = sqrt(x*x + y*y) -0.5; 
        // dfloat signum = 0.5*(1.0 - tanh(M_PI*phi/0.1));
        // bns->q[id+0*mesh->Np] = q1bar*(1.0 + signum); 
        // bns->q[id+1*mesh->Np] = ramp*q2bar;
        // bns->q[id+2*mesh->Np] = ramp*q3bar;
        // bns->q[id+3*mesh->Np] = ramp*ramp*q4bar;
        // bns->q[id+4*mesh->Np] = ramp*ramp*q5bar;
        // bns->q[id+5*mesh->Np] = ramp*ramp*q6bar;   

        // Vortex Problem
        // dfloat r     = sqrt(pow((x-u*time),2) + pow( (y-v*time),2) );
        // dfloat Umax  = 0.5; 
        // dfloat b     = 0.1;

        // dfloat Ur    = Umax/b*r*exp(0.5*(1.0-pow(r/b,2)));

        // dfloat rhor  = rho*exp(-Umax*Umax/(2. * bns->RT) *exp(1.0-r*r/(b*b)));
        // rhor = 1.0;

        // dfloat theta = atan2(y,x);

        // bns->q[id+0*mesh->Np] = rhor*(1./bns->sqrtRT * Gy * y + 1.0); 
        // bns->q[id+1*mesh->Np] = rhor*(-Ur*sin(theta)+u)/bns->sqrtRT;
        // bns->q[id+2*mesh->Np] = rhor*( Ur*cos(theta)+v)/bns->sqrtRT;
        // bns->q[id+3*mesh->Np] = q4bar;
        // bns->q[id+4*mesh->Np] = q5bar;
        // bns->q[id+5*mesh->Np] = q6bar;  


        // Uniform Flow
        bns->q[id+0*mesh->Np] = q1bar; 
        bns->q[id+1*mesh->Np] = ramp*q2bar;
        bns->q[id+2*mesh->Np] = ramp*q3bar;
        bns->q[id+3*mesh->Np] = ramp*ramp*q4bar;
        bns->q[id+4*mesh->Np] = ramp*ramp*q5bar;
        bns->q[id+5*mesh->Np] = ramp*ramp*q6bar;   
      }else{

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
  } else{
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

  
  bns->Vort   = (dfloat*) calloc(3*mesh->Nelements*mesh->Np, sizeof(dfloat));
  bns->o_Vort = mesh->device.malloc(3*mesh->Nelements*mesh->Np*sizeof(dfloat), bns->Vort);
  


  


  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  int maxCubNodes = mymax(maxNodes,mesh->cubNp);

  kernelInfo.addDefine("p_maxNodes", maxNodes);
  kernelInfo.addDefine("p_maxCubNodes", maxCubNodes);


  int NblockV = 128/mesh->Np; // works for CUDA
  NblockV = 1; //!!!!!
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 128/maxNodes; // works for CUDA

  NblockS = 1; // !!!!!
  kernelInfo.addDefine("p_NblockS", NblockS);

  int NblockCub = 128/mesh->cubNp; // works for CUDA

  NblockCub = 1; // !!!!!!!!!!!!!!!!!!!!!

  kernelInfo.addDefine("p_NblockCub", NblockCub);

  // physics 
  kernelInfo.addDefine("p_Lambda2", 0.5f);
  kernelInfo.addDefine("p_sqrtRT", bns->sqrtRT);
  kernelInfo.addDefine("p_sqrt2", (dfloat)sqrt(2.));
  kernelInfo.addDefine("p_isq12", (dfloat)sqrt(1./12.));
  kernelInfo.addDefine("p_isq6", (dfloat)sqrt(1./6.));
  kernelInfo.addDefine("p_invsqrt2", (dfloat)sqrt(1./2.));
  kernelInfo.addDefine("p_tauInv", bns->tauInv);

  dfloat AX = 0, AY = 0, AZ = 0;

  if(options.getArgs("BODYFORCE-X", AX))
    if(AX)
      kernelInfo.addDefine("p_AX", AX/bns->sqrtRT);
  if(options.getArgs("BODYFORCE-Y", AY))
    if(AY)
      kernelInfo.addDefine("p_AY", AY/bns->sqrtRT);

  if(bns->dim==3){
    if(options.getArgs("BODYFORCE-Z", AZ))
      if(AZ)
        kernelInfo.addDefine("p_AZ", AZ/bns->sqrtRT);
   }

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
  char *suffix, *suffixUpdate;
  
  if(bns->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(bns->elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(bns->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(bns->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  if(bns->elementType==TRIANGLES || bns->elementType==QUADRILATERALS)
    suffixUpdate = strdup("2D");
  if(bns->elementType==TETRAHEDRA || bns->elementType==HEXAHEDRA)
    suffixUpdate = strdup("3D");
  



  char fileName[BUFSIZ], kernelName[BUFSIZ];
  for (int r=0;r<size;r++){

    if (r==rank) {

      // Volume kernels
      sprintf(fileName, DBNS "/okl/bnsVolume%s.okl", suffix);

      sprintf(kernelName, "bnsVolume%s", suffix);
      bns->volumeKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);
      sprintf(kernelName, "bnsPmlVolume%s", suffix);
      bns->pmlVolumeKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);
     
      // Relaxation kernels
      sprintf(fileName, DBNS "/okl/bnsRelaxation%s.okl", suffix);

      sprintf(kernelName, "bnsRelaxation%s", suffix);
      bns->relaxationKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);
      sprintf(kernelName, "bnsPmlRelaxation%s", suffix);        
      bns->pmlRelaxationKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);
      
      // Surface kernels 
      sprintf(fileName, DBNS "/okl/bnsSurface%s.okl", suffix);

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

      
      sprintf(fileName, DBNS "/okl/bnsUpdate%s.okl", suffixUpdate);
      // Update Kernels
      if(options.compareArgs("TIME INTEGRATOR","LSERK")){
        sprintf(kernelName, "bnsLSERKUpdate%s", suffixUpdate);
        bns->updateKernel = mesh->device.buildKernelFromSource(fileName, kernelName,kernelInfo);

        sprintf(kernelName, "bnsLSERKPmlUpdate%s", suffixUpdate);
        bns->pmlUpdateKernel = mesh->device.buildKernelFromSource(fileName, kernelName,kernelInfo);
      } else if(options.compareArgs("TIME INTEGRATOR","SARK")){
        sprintf(kernelName, "bnsSARKUpdateStage%s", suffixUpdate);
        bns->updateStageKernel = mesh->device.buildKernelFromSource(fileName,kernelName, kernelInfo);

        sprintf(kernelName, "bnsSARKPmlUpdateStage%s", suffixUpdate);
        bns->pmlUpdateStageKernel = mesh->device.buildKernelFromSource(fileName,kernelName, kernelInfo);

        sprintf(kernelName, "bnsSARKUpdate%s", suffixUpdate);
        bns->updateKernel = mesh->device.buildKernelFromSource(fileName, kernelName,kernelInfo);

        sprintf(kernelName, "bnsSARKPmlUpdate%s", suffixUpdate);
        bns->pmlUpdateKernel = mesh->device.buildKernelFromSource(fileName, kernelName,kernelInfo);
        
        if(bns->fixed_dt==0){
          sprintf(fileName, DBNS "/okl/bnsErrorEstimate.okl");
          sprintf(kernelName, "bnsErrorEstimate");
          bns->errorEstimateKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);
        }
      } else if(options.compareArgs("TIME INTEGRATOR","MRSAAB")){
      
        sprintf(kernelName, "bnsMRSAABTraceUpdate%s", suffixUpdate);
        bns->traceUpdateKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "bnsMRSAABUpdate%s", suffixUpdate);
        bns->updateKernel = mesh->device.buildKernelFromSource(fileName, kernelName,kernelInfo);

        sprintf(kernelName, "bnsMRSAABPmlUpdate%s", suffixUpdate);
        bns->pmlUpdateKernel = mesh->device.buildKernelFromSource(fileName, kernelName,kernelInfo);
      }

      sprintf(fileName, DBNS "/okl/bnsVorticity%s.okl",suffix);
      sprintf(kernelName, "bnsVorticity%s", suffix);
      bns->vorticityKernel = mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);


      // This needs to be unified
      mesh->haloExtractKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract3D.okl",
                                           "meshHaloExtract3D",
                                           kernelInfo);




    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  return bns; 
}




