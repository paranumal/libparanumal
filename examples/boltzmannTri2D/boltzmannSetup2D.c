#include "boltzmann2D.h"

// NBN: toggle use of 2nd stream
#define USE_2_STREAMS

bns_t *boltzmannSetup2D(mesh2D *mesh, char * options){

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  
  // SET SOLVER PARAMETERS
  bns_t *bns = (bns_t*) calloc(1, sizeof(bns_t));
  
  mesh->Nfields = 6;
  bns->Nfields = mesh->Nfields; 

  bns->mesh = mesh; 

  bns->errorStep = 1000; 

    // Initialize
  dfloat Ma      = 0.f,   Re      = 0.f;
  dfloat rho     = 1.f,   u       = 0.f;
  dfloat v       = 0.f,   nu      = 0.f;
  dfloat Uref    = 0.2f,  Lref    = 1.f; 
  dfloat sigma11 = 0.f ,  sigma12 = 0.f, sigma22 = 0.f;  

  if(strstr(options, "PML")){
    printf("Starting initial conditions for PML\n");
    bns->Ma = 0.1;   // 0.17; //MA[mesh->Ntscale]; // Set Mach number
    bns->Re = 1000.;  //1000.; //RE[mesh->Ntscale]; 
    //
    Uref = 1.0; //  0.2;
    Lref = 1.0;   // set Lref
    //
    bns->RT       = Uref*Uref/(bns->Ma*bns->Ma);
    bns->sqrtRT   = sqrt(bns->RT);  
    //
    nu            = Uref*Lref/bns->Re; 
    bns->tauInv   = bns->RT/nu;

    //printf("starting initial conditions\n"); //Zero Flow Conditions
    // rho = 1.0; u = Uref*cos(M_PI/6.0); v = Uref*sin(M_PI/6.0); sigma11 = 0; sigma12 = 0; sigma22 = 0;
    rho = 1.0; u = 1.0; v = 0.0; sigma11 = 0; sigma12 = 0; sigma22 = 0;
    // rho = 1.0; u = Uref*cos(M_PI/6); v = Uref*sin(M_PI/6); sigma11 = 0; sigma12 = 0; sigma22 = 0;
    //
    bns->startTime = 0.0; 
    bns->finalTime = 20.; // Was 200  
  }
  else{
    printf("Starting initial conditions for NONPML\n");
    bns->Ma = 0.1;     //Set Mach number
    bns->Re = 1000.;   // Set Reynolds number
    
    Uref = 1.;   // Set Uref
    Lref = 1.;   // set Lref
    //
    bns->RT  = Uref*Uref/(bns->Ma*bns->Ma);
    bns->sqrtRT = sqrt(bns->RT);  
    //
    nu = Uref*Lref/bns->Re; 
    bns->tauInv = bns->RT/nu;
    
    #if 0
    // Create Periodic Boundaries
    printf("Creating periodic connections if exist \n");
    dfloat xper = 1.0;   dfloat yper = 0.0;
    boltzmannPeriodic2D(mesh,xper,yper);
    #endif

    //printf("starting initial conditions\n"); //Zero Flow Conditions
    rho = 1.0, u = Uref, v = 0.; sigma11 = 0, sigma12 = 0, sigma22 = 0;
    //
    bns->startTime = 0.0; 
    bns->finalTime = 100.;
  }

  // set penalty parameter
  bns->Lambda2 = 0.5/(bns->sqrtRT);


  dfloat time = bns->startTime + 0.0; 
  // Define Initial Mean Velocity
  dfloat ramp, drampdt;
  boltzmannRampFunction2D(time, &ramp, &drampdt);
  
  dfloat q1bar = rho;
  dfloat q2bar = rho*u/bns->sqrtRT;
  dfloat q3bar = rho*v/bns->sqrtRT;
  dfloat q4bar = (rho*u*v - sigma12)/bns->RT;
  dfloat q5bar = (rho*u*u - sigma11)/(sqrt(2.)*bns->RT);
  dfloat q6bar = (rho*v*v - sigma22)/(sqrt(2.)*bns->RT);

  // SET STABLE TIME STEP SIZE
  dfloat cfl = 0.25; 
  if(strstr(options,"LSERK") || strstr(options,"SARK") || strstr(options,"LSIMEX") || 
     strstr(options,"DOPRI5") || strstr(options,"XDOPRI") || strstr(options,"SAADRK") )
    cfl          = 0.75; 


  dfloat magVelocity  = sqrt(q2bar*q2bar+q3bar*q3bar)/(q1bar/bns->sqrtRT);
  magVelocity         = mymax(magVelocity,1.0); // Correction for initial zero velocity

  dfloat ghmin        = 1e9; 
  dfloat dt           = 1e9; 
  dfloat *EtoDT       = (dfloat *) calloc(mesh->Nelements,sizeof(dfloat));

  //Set time step size
  for(int e=0;e<mesh->Nelements;++e)
  { 
    dfloat hmin = 1e9, dtmax = 1e9;
    
    EtoDT[e] = dtmax;

    for(int f=0;f<mesh->Nfaces;++f){
      int sid    = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];
     
      dfloat htest = 2.0/(sJ*invJ);
      hmin = mymin(hmin, htest); 
    }

    ghmin = mymin(ghmin, hmin);

    dfloat dtex   = cfl*hmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*bns->sqrtRT);
    dfloat dtim   = cfl*1.f/(bns->tauInv);
    dfloat dtest = 1e9;
    
    if(strstr(options,"MRAB") || strstr(options,"LSERK") || strstr(options,"SRAB") || strstr(options,"DOPRI5") )
      dtest  = mymin(dtex,dtim); // For fully expliciy schemes
    else if(strstr(options,"MRSAAB") || strstr(options,"SARK")   || strstr(options,"SAAB") || 
            strstr(options,"LSIMEX") || strstr(options,"XDOPRI") || strstr(options,"SAADRK"))
      dtest  = dtex; // For semi analytic  and IMEX
    else
      printf("time discretization method should be choosen\n");
      
      dt            = mymin(dt, dtest);      // For SR 
      EtoDT[e]      = mymin(EtoDT[e],dtest); // For MR
  }

   printf("dtex = %.5e dtim = %.5e \n", cfl*ghmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*bns->sqrtRT), cfl*1.f/(bns->tauInv));

 
  // Set multiRate/singleRate element groups/group  
  if(strstr(options, "MRAB") || strstr(options,"MRSAAB")){
    int maxLevels = 100;
    bns->dt = meshMRABSetup2D(mesh,EtoDT,maxLevels, bns->finalTime- bns->startTime);
    bns->NtimeSteps =  (bns->finalTime-bns->startTime)/(pow(2,mesh->MRABNlevels-1)*bns->dt);

  }
  else{    
    //!!!!!!!!!!!!!! Fix time step to compute the error in postprecessing step  
    // dt = 100e-6; // !!!!!!!!!!!!!!!!
    // MPI_Allreduce to get global minimum dt
    MPI_Allreduce(&dt, &(bns->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
    bns->NtimeSteps = (bns->finalTime-bns->startTime)/bns->dt;
    bns->dt         = (bns->finalTime-bns->startTime)/bns->NtimeSteps;

    //offset index
    bns->shiftIndex = 0;

    // Set element ids for nonPml region, will be modified if PML exists
    mesh->nonPmlNelements = mesh->Nelements; 
    mesh->pmlNelements    = 0; 

    mesh->nonPmlElementIds = (int*) calloc(mesh->Nelements, sizeof(int));
    for(int e=0;e<mesh->Nelements;++e)
     mesh->nonPmlElementIds[e] = e; 
  }


 // INITIALIZE FIELD VARIABLES 
 if(strstr(options, "MRAB") || strstr(options,"MRSAAB")){
    bns->Nrhs = 3;
    // compute samples of q at interpolation nodes
    bns->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rhsq = (dfloat*) calloc(bns->Nrhs*mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->fQM  = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*bns->Nfields, sizeof(dfloat));
  }

  else if(strstr(options,"SRSAAB") || strstr(options,"SRAB") || strstr(options,"SARK")){ //SRAB and SAAB Single rate versions of above
    bns->Nrhs = 3;
    // compute samples of q at interpolation nodes
    bns->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rhsq = (dfloat*) calloc(bns->Nrhs*mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
  }

  else if (strstr(options,"LSERK")){ // LSERK, SARK etc.
    bns->Nrhs = 1; 
    // compute samples of q at interpolation nodes
    bns->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rhsq = (dfloat*) calloc(bns->Nrhs*mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->resq = (dfloat*) calloc(bns->Nrhs*mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
  }


  else if (strstr(options,"LSIMEX")){ // LSERK, SARK etc.
    bns->Nrhs = 1; 
    // compute samples of q at interpolation nodes
    bns->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rhsq = (dfloat*) calloc(bns->Nrhs*mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
  }


  else if (strstr(options,"DOPRI5") || strstr(options,"XDOPRI")){ // LSERK, SARK etc.
    bns->Nrhs = 1; 

    int Ntotal    = mesh->Nelements*mesh->Np*bns->Nfields;
    bns->Nblock   = (Ntotal+blockSize-1)/blockSize;

    int localElements =  mesh->Nelements;
    MPI_Allreduce(&localElements, &(bns->totalElements), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // compute samples of q at interpolation nodes
    bns->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    //
    bns->NrkStages = 7;
    
    bns->rkq      = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rkrhsq   = (dfloat*) calloc(bns->NrkStages*mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rkerr    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->errtmp  =  (dfloat*) calloc(bns->Nblock, sizeof(dfloat)); 
  }


  else if (strstr(options,"SAADRK")){ 
    bns->Nrhs = 1; 

    int Ntotal    = mesh->Nelements*mesh->Np*bns->Nfields;
    bns->Nblock   = (Ntotal+blockSize-1)/blockSize;

    int localElements =  mesh->Nelements;
    MPI_Allreduce(&localElements, &(bns->totalElements), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // compute samples of q at interpolation nodes
    bns->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    //
    bns->NrkStages = 7; // 5 or 7 // SAADRK methods 43 or 53 respectively
    
    bns->rkq      = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rkrhsq   = (dfloat*) calloc(bns->NrkStages*mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rkerr    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->errtmp  =  (dfloat*) calloc(bns->Nblock, sizeof(dfloat)); 
  }





 // INITIALIZE PROBLEM 
  int cnt = 0;
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat t = 0;
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];
#if 0
      // Couette Flow 
      dfloat uex = y ; 
      dfloat sex = 0.0; 
      for(int k=1; k<=10; k++){
        dfloat lamda = k*M_PI;
        // dfloat coef = -bns->RT*bns->tauInv/2. + sqrt(pow((bns->RT*bns->tauInv),2) /4.0 - (lamda*lamda*bns->RT*bns->RT));
        dfloat coef = -bns->tauInv/2. + bns->tauInv/2. * sqrt(1. - 4.*pow(1./ bns->tauInv, 2)* bns->RT*lamda*lamda);
        uex += 2.*pow(-1,k)/(lamda)*exp(coef*time)*sin(lamda*y); //
        sex  += 2.*pow(-1,k)*coef/(lamda*lamda)*exp(coef*time)*cos(lamda*y); 
      }

      bns->q[cnt+0] = rho; // uniform density, zero flow
      bns->q[cnt+1] = rho*uex/bns->sqrtRT;
      bns->q[cnt+2] = 0.0;
      bns->q[cnt+3] = sex/bns->RT;
      bns->q[cnt+4] = 0.0;
      bns->q[cnt+5] = 0.0;  
#endif


#if 0
      // Vortex Problem
      dfloat r     = sqrt(pow((x-u*time),2) + pow( (y-v*time),2) );
      dfloat Umax  = 0.5*u; 
      dfloat b     = 0.2;

      dfloat Ur    = Umax/b*r*exp(0.5*(1.0-pow(r/b,2)));

      dfloat rhor  = rho*exp(-Umax*Umax/(2. * bns->RT) *exp(1.0-r*r/(b*b)));

      dfloat theta = atan2(y,x);

      bns->q[cnt+0] = rhor; 
      bns->q[cnt+1] = rhor*(-Ur*sin(theta) +u)/bns->sqrtRT;
      bns->q[cnt+2] = rhor*( Ur*cos(theta) +v)/bns->sqrtRT;
      bns->q[cnt+3] = q4bar;
      bns->q[cnt+4] = q5bar;
      bns->q[cnt+5] = q6bar;  
    
#endif

#if 1
      // Uniform Flow
      bns->q[cnt+0] = q1bar; 
      bns->q[cnt+1] = ramp*q2bar;
      bns->q[cnt+2] = ramp*q3bar;
      bns->q[cnt+3] = ramp*ramp*q4bar;
      bns->q[cnt+4] = ramp*ramp*q5bar;
      bns->q[cnt+5] = ramp*ramp*q6bar;  
#endif

      cnt += bns->Nfields;
    }
  }

  
  // Write Problem Info 
  if(rank==0){
    printf("dt   = %g ",   bns->dt);
    printf("max wave speed = %g\n", sqrt(3.)*bns->sqrtRT);
    if(mesh->MRABNlevels)
      printf("Nsteps = %d dt = %.8e MRAB Level: %d  Final Time:%.5e\n", bns->NtimeSteps, bns->dt, mesh->MRABNlevels, bns->startTime+pow(2, mesh->MRABNlevels-1)*(bns->dt*(bns->NtimeSteps+1)));   
    else
     printf("Nsteps = %d dt = %.8e Final Time:%.5e\n", bns->NtimeSteps, bns->dt,  bns->startTime + bns->dt*bns->NtimeSteps);
  }
 
 
  // SET PROBE DATA
  if(strstr(options, "PROBE")){
    mesh->probeNTotal = 3; 
    dfloat *pX   = (dfloat *) calloc (mesh->probeNTotal, sizeof(dfloat));
    dfloat *pY   = (dfloat *) calloc (mesh->probeNTotal, sizeof(dfloat));
    // Fill probe coordinates
     pX[0] = 9.00;  pX[1] = 10.00;  pX[2] = 10.50; //pX[3] = 5;
     pY[0] = 5.00;  pY[1] =  6.00;  pY[2] =  6.50; //pY[3] = 5*tan(M_PI/6.);

    meshProbeSetup2D(mesh, pX, pY);

    free(pX); free(pY);

  }
  


  //OCCCA SETUP
  char deviceConfig[BUFSIZ];
  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", rank%2);
  //sprintf(deviceConfig, "mode = OpenCL, deviceID = 1, platformID = 0");
  //sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");  
 

  occa::kernelInfo kernelInfo;
  meshOccaSetup2D(mesh, deviceConfig,  kernelInfo);  

  kernelInfo.addParserFlag("automate-add-barriers", "disabled");   

  // Setup MRAB PML
  if(strstr(options, "MRAB") || strstr(options,"MRSAAB")){
     printf("Preparing Pml for multirate rate\n");
    boltzmannMRABPmlSetup2D(bns, options);

    mesh->o_MRABelementIds = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    mesh->o_MRABhaloIds    = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    for (int lev=0;lev<mesh->MRABNlevels;lev++) {
      if (mesh->MRABNelements[lev])
        mesh->o_MRABelementIds[lev] = mesh->device.malloc(mesh->MRABNelements[lev]*sizeof(int),
        mesh->MRABelementIds[lev]);
      if (mesh->MRABNhaloElements[lev])
        mesh->o_MRABhaloIds[lev] = mesh->device.malloc(mesh->MRABNhaloElements[lev]*sizeof(int),
        mesh->MRABhaloIds[lev]);
    }
  } 
  else{
   printf("Preparing Pml for single rate\n");
   boltzmannPmlSetup2D(bns, options); 

    if (mesh->nonPmlNelements)
        mesh->o_nonPmlElementIds = mesh->device.malloc(mesh->nonPmlNelements*sizeof(int), mesh->nonPmlElementIds);

  }
  


#if 0
mesh2D *meshsave = mesh;
int fld = 0; 
  for(int lev=0; lev<mesh->MRABNlevels; lev++){

        for(int m=0; m<mesh->MRABNelements[lev]; m++){
          int e = mesh->MRABelementIds[lev][m];
          for(int n=0; n<mesh->Np; n++){
             bns->q[bns->Nfields*(n + e*mesh->Np) + fld] = lev;
          }
        }


        for(int m=0; m<mesh->MRABpmlNelements[lev]; m++){
          int e = mesh->MRABpmlElementIds[lev][m];
          for(int n=0; n<mesh->Np; n++){
             bns->q[bns->Nfields*(n + e*mesh->Np) + fld] = lev;
          }
        }
  }

 char fname[BUFSIZ];
 sprintf(fname, "ElementGroupsMRSAAB.dat");
 // meshPlotVTU2D(mesh, fname,fld);
 boltzmannPlotFieldTEC2D(mesh, fname,0, fld);


mesh = meshsave; 

#endif



#if 1

// Compute Time Stepper Coefficcients

boltzmannTimeStepperCoefficients(bns, options);


if(strstr(options,"MRSAAB") || strstr(options,"MRAB") || 
   strstr(options,"SRAB")   || strstr(options,"SAAB") ){

  //printf("Setting Rhs for AB\n");
  // bns->o_rhsq.free();
  bns->o_q     = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat),bns->q);
  bns->o_rhsq  = mesh->device.malloc(bns->Nrhs*mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->rhsq);

  if(strstr(options,"MRSAAB") || strstr(options, "MRAB")){

    if (mesh->totalHaloPairs) {
      //reallocate halo buffer for trace exchange
      mesh->o_haloBuffer.free();
      mesh->o_haloBuffer = mesh->device.malloc(mesh->totalHaloPairs*mesh->Nfp*bns->Nfields*mesh->Nfaces*sizeof(dfloat));
    }

  bns->o_fQM = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*bns->Nfields*sizeof(dfloat),
                                  bns->fQM);
  mesh->o_mapP = mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(int), mesh->mapP);
  }
}


if(strstr(options, "LSERK")){
  // 
  bns->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->q);
  bns->o_rhsq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->rhsq);
  bns->o_resq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->resq);
}

if(strstr(options,"DOPRI5") || strstr(options,"XDOPRI") || strstr(options,"SAADRK")){

  bns->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->q);
  bns->o_rhsq = 
    mesh->device.malloc(bns->Nrhs*mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->rhsq); 
  
  int Ntotal    = mesh->Nelements*mesh->Np*bns->Nfields;
  
  bns->o_rkq =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->rkq);
  bns->o_rkrhsq =
    mesh->device.malloc(bns->NrkStages*mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->rkrhsq);
  bns->o_rkerr =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->rkerr);
  bns->o_errtmp = mesh->device.malloc(bns->Nblock*sizeof(dfloat), bns->errtmp);

  bns->o_rkA = mesh->device.malloc(bns->NrkStages*bns->NrkStages*sizeof(dfloat), bns->rkA);
  bns->o_rkE = mesh->device.malloc(bns->NrkStages*sizeof(dfloat), bns->rkE);

}

if(strstr(options, "SARK")){
  bns->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->q);
  bns->o_rhsq = 
    mesh->device.malloc(bns->Nrhs*mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->rhsq); 
  bns->o_resq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->resq);
  bns->o_qS =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat),bns->q);
  }
 
  


else if(strstr(options, "LSIMEX")){ 
  bns->Nimex = 4;
  bns->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->q);
  bns->o_rhsq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->rhsq);
  // bns->o_resq =
  //   mesh->device.malloc(mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->resq);
  bns->o_qY =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->q);

  bns->o_qZ =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->q);
}  



#endif 


  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  int maxCubNodes = mymax(maxNodes,mesh->cubNp);

  kernelInfo.addDefine("p_maxNodes", maxNodes);
  kernelInfo.addDefine("p_maxCubNodes", maxCubNodes);


  int NblockV = 128/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 128/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  int NblockCub = 128/mesh->cubNp; // works for CUDA
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
  kernelInfo.addDefine("p_alpha0", (dfloat).01f);
  kernelInfo.addDefine("p_pmlAlpha", (dfloat)0.1f);
  kernelInfo.addDefine("p_blockSize", blockSize);

  kernelInfo.addDefine("p_NrkStages", bns->NrkStages);



 // Volume and Relaxation Kernels
  if(strstr(options, "CUBATURE") && !strstr(options,"LSIMEX")){ 
          
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

    if(strstr(options, "MRAB") || strstr(options, "LSERK") || 
       strstr(options,"SRAB")  || strstr(options,"DOPRI5") ){ 

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
    else if(strstr(options, "MRSAAB")||strstr(options, "SAAB") || 
            strstr(options,"SARK")   || strstr(options,"XDOPRI") || strstr(options,"SAADRK")){

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



  if(strstr(options, "COLLOCATION")  && !strstr(options,"LSIMEX")){ 

     if(strstr(options, "MRAB") || strstr(options, "LSERK") || strstr(options,"SRAB") || strstr(options,"DOPRI5")){ 

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

    else if(strstr(options, "MRSAAB") || strstr(options, "SAAB") || 
            strstr(options,"SARK") || strstr(options,"XDOPRI") || strstr(options,"SAADRK")){
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
  if(strstr(options, "MRAB")){ 
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

    else if(strstr(options, "MRSAAB")){

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

     else if(strstr(options,"SRAB")){
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


     else if(strstr(options,"SRSAAB")){
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

     else if(strstr(options,"LSERK")){
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



    else if(strstr(options,"SARK")){
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
    else if(strstr(options,"DOPRI5")){

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


    else if(strstr(options,"XDOPRI")){
    printf("compiling XDOPRI Update kernels\n");
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


     else if(strstr(options,"SAADRK")){

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
    if(strstr(options,"MRAB") || strstr(options,"MRSAAB")){  
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
    if(strstr(options, "LSIMEX")){ 
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
      
     
    if(strstr(options, "CUBATURE")){ 
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
    else if(strstr(options, "COLLOCATION")){ 
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

      

    mesh->haloExtractKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
          "meshHaloExtract2D",
             kernelInfo);


return bns; 
  
}




