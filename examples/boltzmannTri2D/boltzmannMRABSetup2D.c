#include "boltzmann2D.h"

#define USE_2_STREAMS

void boltzmannMRABSetup2D(mesh2D *mesh, char * options){

  iint rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  // SET SOLVER PHYSICAL PARAMETERS
  
  mesh->Nfields = 6;

 // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields, sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(3*mesh->Nelements*mesh->Np*mesh->Nfields, sizeof(dfloat));

  mesh->fQM  = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*mesh->Nfields, sizeof(dfloat));
  mesh->fQP  = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*mesh->Nfields, sizeof(dfloat));




 // Initialize
  dfloat Ma      = 0.f,   Re     = 0.f;
  dfloat rho     = 1.f,   u      = 0.f;
  dfloat v       = 0.f,   nu     = 0.f;
  dfloat Uref    = 1.f,   Lref   = 1.f; 
  dfloat sigma11 = 0.f , sigma12 = 0.f, sigma22 = 0.f;  

  if(strstr(options, "PML")){
    printf("Starting initial conditions for PML\n");
    Ma = 0.01;    //Set Mach number
    Re = 1000.;  // Set Reynolds number
    //
    Uref = 1.;   // Set Uref
    Lref = 1.;   // set Lref
    //
    mesh->RT      = Uref*Uref/(Ma*Ma);
    mesh->sqrtRT  = sqrt(mesh->RT);  
    //
    nu            = Uref*Lref/Re; 
    mesh->tauInv  = mesh->RT/nu;

    //printf("starting initial conditions\n"); //Zero Flow Conditions
    rho = 1., u = Uref; v = 0.; sigma11 = 0, sigma12 = 0, sigma22 = 0;
    //
    mesh->finalTime = 30.0;
  }
  else{
    printf("Starting initial conditions for NONPML\n");
    Ma = 0.1;     //Set Mach number
    Re = 1000.;   // Set Reynolds number
    //
    Uref = 1.;   // Set Uref
    Lref = 1.;   // set Lref
    //
    mesh->RT  = Uref*Uref/(Ma*Ma);
    mesh->sqrtRT = sqrt(mesh->RT);  
    //
    nu = Uref*Lref/Re; 
    mesh->tauInv = mesh->RT/nu;
    
    #if 0
    // Create Periodic Boundaries
    printf("Creating periodic connections if exist \n");
    dfloat xper = 1.0;   dfloat yper = 0.0;
    boltzmannPeriodic2D(mesh,xper,yper);
    #endif

    //printf("starting initial conditions\n"); //Zero Flow Conditions
    rho = 1., u = 0., v = 0.; sigma11 = 0, sigma12 = 0, sigma22 = 0;
    //
    mesh->finalTime = 200.;
  }

  // set penalty parameter
  mesh->Lambda2 = 0.5/(mesh->sqrtRT);

  // Define Initial Mean Velocity
  dfloat ramp, drampdt;
  boltzmannRampFunction2D(0, &ramp, &drampdt);


  dfloat q1bar = rho;
  dfloat q2bar = rho*u/mesh->sqrtRT;
  dfloat q3bar = rho*v/mesh->sqrtRT;
  dfloat q4bar = (rho*u*v - sigma12)/mesh->RT;
  dfloat q5bar = (rho*u*u - sigma11)/(sqrt(2.)*mesh->RT);
  dfloat q6bar = (rho*v*v - sigma22)/(sqrt(2.)*mesh->RT);


  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      dfloat t = 0;
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];

      mesh->q[cnt+0] = q1bar; // uniform density, zero flow
      mesh->q[cnt+1] = ramp*q2bar;
      mesh->q[cnt+2] = ramp*q3bar;
      mesh->q[cnt+3] = ramp*ramp*q4bar;
      mesh->q[cnt+4] = ramp*ramp*q5bar;
      mesh->q[cnt+5] = ramp*ramp*q6bar;  
      cnt += mesh->Nfields;
    }
  }



  
  // Set stable time step size for each element
  dfloat cfl          = 0.3; 
  dfloat magVelocity  = sqrt(q2bar*q2bar+q3bar*q3bar)/(q1bar/mesh->sqrtRT);
  magVelocity         = mymax(magVelocity,1.0); // Correction for initial zero velocity

  dfloat *EtoDT       = (dfloat *) calloc(mesh->Nelements,sizeof(dfloat));

  //Set time step size
  for(iint e=0;e<mesh->Nelements;++e)
  { 
    dfloat hmin = 1e9, dtmax = 1e9;
    
    EtoDT[e] = dtmax;

    for(iint f=0;f<mesh->Nfaces;++f){
      iint sid    = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];
     
      dfloat htest = 0.5/(sJ*invJ);
      hmin = mymin(hmin, htest); 
    }

    dfloat dtex   = cfl*hmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*mesh->sqrtRT);
    dfloat dtim   = 1.f/(mesh->tauInv);

    dfloat dtest = 1e9;
    
    if(strstr(options,"MRAB"))
      dtest  = mymin(dtex,dtim); // For fully expliciy schemes
    else if(strstr(options,"MRSAAB"))
      dtest  = dtex; // For semi analytic and IMEX 
    else
      printf("MRAB or MRSAAB time discretization methods should be choosen\n");
      

    EtoDT[e]      = mymin(EtoDT[e],dtest);
  }


  int maxLevels = 100;
  meshMRABSetup2D(mesh,EtoDT,maxLevels);
  
 


  

     
  // if(strstr(options, "MR_GROUPS")){
  //   iint cnt1 = 0;
  //  for (iint e=0;e<mesh->Nelements;e++){
  //     for(iint n=0;n<mesh->Np;n++){
  //     mesh->q[cnt1+0] = mesh->MRABlevel[e]; 
  //     cnt1 += mesh->Nfields;
  //     }
  //    }
  //   mesh->sqrtRT = 1.0; iint zero = 0;
  //   char fname[BUFSIZ];
  //   sprintf(fname, "MR_GROUPS_%04d_%04d.vtu", rank, 0);
  //   meshPlotVTU2D(mesh, fname, zero);
  // }
 

  

  // dfloat hmin =0.0, hmax= 0.0; 
  // MPI_Allreduce(&hgmax, &hmax, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);    
  // MPI_Allreduce(&hgmin, &hmin, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);    

  mesh->errorStep = 10000;
  if(rank==0){
  printf("dt   = %g ",   mesh->dt);
  printf("max wave speed = %g\n", sqrt(3.)*mesh->sqrtRT);
  printf("Nsteps = %d dt = %.8e MRAB Level: %d  Final Time:%.5e\n", mesh->NtimeSteps, mesh->dt, mesh->MRABNlevels, pow(2, mesh->MRABNlevels-1)*(mesh->dt*mesh->NtimeSteps+1));

  }
 


  if(strstr(options,"MRSAAB")){

  mesh->MRSAAB_A = (dfloat *) calloc(3*3*mesh->MRABNlevels,sizeof(dfloat));
  mesh->MRSAAB_B = (dfloat *) calloc(3*3*mesh->MRABNlevels,sizeof(dfloat));
  mesh->MRSAAB_C = (dfloat *) calloc(    mesh->MRABNlevels,sizeof(dfloat));
  mesh->MRAB_A   = (dfloat *) calloc(3*3*mesh->MRABNlevels,sizeof(dfloat));
  mesh->MRAB_B   = (dfloat *) calloc(3*3*mesh->MRABNlevels,sizeof(dfloat));
  mesh->MRAB_C   = (dfloat *) calloc(    mesh->MRABNlevels,sizeof(dfloat));

  iint MRABorder = 3; 

  for(iint l = 0; l<mesh->MRABNlevels; ++l){
    // MRSAAB coefficients
    dfloat cc = -mesh->tauInv;
    dfloat h  = mesh->dt * pow(2,l); 
    //
    for (iint order=0; order<3; ++order){
      // computation of coefficients based on magnitude
      const iint id = order*mesh->MRABNlevels*3 + l*3;

      if(order==0){
        if(fabs(cc*h)>1e-2){  // Use exponentials
           // Full dt coeeficients
          mesh->MRSAAB_A[id + 0] = (exp(cc*h) - 1.)/cc;
          mesh->MRSAAB_A[id + 1] = 0.f;
          mesh->MRSAAB_A[id + 2] = 0.f;
          // Half coefficients
          mesh->MRSAAB_B[id + 0] = (exp((cc*h)/2.) - 1.)/cc;
          mesh->MRSAAB_B[id + 1] = 0.f;
          mesh->MRSAAB_B[id + 2] = 0.f;
        }
        else{
          // Full dt coeeficients
          mesh->MRSAAB_A[id + 0] = (pow(cc,3)*pow(h,4))/24. + (pow(cc,2)*pow(h,3))/6. + (cc*pow(h,2))/2. + h; 
          mesh->MRSAAB_A[id + 1] = 0.f;
          mesh->MRSAAB_A[id + 2] = 0.f;
          // Half coefficients
          mesh->MRSAAB_B[id + 0] = (pow(cc,3)*pow(h,4))/384. + (pow(cc,2)*pow(h,3))/48. + (cc*pow(h,2))/8. + h/2.;
          mesh->MRSAAB_B[id + 1] = 0.f;
          mesh->MRSAAB_B[id + 2] = 0.f;
          }

        // MRAB coefficients
        mesh->MRAB_A[id + 0]   =  h ;
        mesh->MRAB_A[id + 1]   =  0.f ;
        mesh->MRAB_A[id + 2]   =  0.f ;

        mesh->MRAB_B[id+0]     =  h/2. ;
        mesh->MRAB_B[id+1]     =  0.f ;
        mesh->MRAB_B[id+2]     =  0.f ;
      }

      else if(order==1){
        if(fabs(cc*h)>1e-2){  // Use exponentials
         // Full dt coeeficients
          mesh->MRSAAB_A[id + 0] = (exp(cc*h) - 2.*cc*h + cc*h*exp(cc*h) - 1.)/(pow(cc,2)*h);
          mesh->MRSAAB_A[id + 1] = (cc*h - exp(cc*h) + 1.)/(pow(cc,2)*h);
          mesh->MRSAAB_A[id + 2] = 0.f;
          // Half coefficients
          mesh->MRSAAB_B[id + 0] = (exp((cc*h)/2.) - (3*cc*h)/2. + cc*h*exp((cc*h)/2.) - 1.)/(pow(cc,2)*h);
          mesh->MRSAAB_B[id + 1] = ((cc*h)/2. - exp((cc*h)/2.) + 1.)/(pow(cc,2)*h);
          mesh->MRSAAB_B[id + 2] = 0.f;
        }
        else{
         // Full dt coeeficients
          mesh->MRSAAB_A[id + 0] = (pow(cc,3)*pow(h,4))/20. + (5.*pow(cc,2)*pow(h,3))/24. + (2.*h*pow(h,2))/3. + (3.*h)/2.; 
          mesh->MRSAAB_A[id + 1] = - h/2. - (cc*pow(h,2))/6. - (pow(cc,2)*pow(h,3))/24. - (pow(cc,3)*pow(h,4))/120. ; 
          mesh->MRSAAB_A[id + 2] = 0.f;
          // Half coefficients
          mesh->MRSAAB_B[id + 0] = (11.*pow(cc,3)*pow(h,4))/3840. + (3.*pow(cc,2)*pow(h,3))/128. + (7.*cc*pow(h,2))/48. + (5.*h)/8.;
          mesh->MRSAAB_B[id + 1] = - h/8. - (cc*pow(h,2))/48. - (pow(cc,2)*pow(h,3))/384. - (pow(cc,3)*pow(h,4))/3840.;
          mesh->MRSAAB_B[id + 2] = 0.f; 
        }

        // MRAB coefficients
        mesh->MRAB_A[id + 0]   =  3.*h/2. ;
        mesh->MRAB_A[id + 1]   = -1.*h/2. ;
        mesh->MRAB_A[id + 2]   =  0.f ;

        mesh->MRAB_B[id + 0]   =  5.*h/8. ;
        mesh->MRAB_B[id + 1]   = -1.*h/8. ;
        mesh->MRAB_B[id + 2]   =   0.f ;
      }

      else{
        if(fabs(cc*h)>1e-2){  // Use exponentials
          // Full dt coeeficients
          mesh->MRSAAB_A[id+0] = (exp(cc*h) - (5.*cc*h)/2. - 3.*pow(cc,2)*pow(h,2 )+ pow(cc,2)*pow(h,2)*exp(cc*h) + (3.*cc*h*exp(cc*h))/2. - 1.)/(pow(cc,3)*pow(h,2));
          mesh->MRSAAB_A[id+1] = (4.*cc*h - 2.*exp(cc*h) + 3.*pow(cc,2)*pow(h,2 )- 2.*cc*h*exp(cc*h) + 2.)/(pow(cc,3)*pow(h,2));
          mesh->MRSAAB_A[id+2] = -((3.*cc*h)/2. - exp(cc*h) + pow(cc,2)*pow(h,2 )- (cc*h*exp(cc*h))/2. + 1.)/(pow(cc,3)*pow(h,2));
          // Half coefficients
          mesh->MRSAAB_B[id+0] = (exp((cc*h)/2.) - 2.*cc*h - (15.*pow(cc,2)*pow(h,2))/8.f + pow(cc*h,2)*exp((cc*h)/2.) + (3.*cc*h*exp((cc*h)/2.))/2. - 1.)/(pow(cc,3)*pow(h,2.));
          mesh->MRSAAB_B[id+1] = (3.*cc*h - 2.*exp((cc*h)/2.0) + (5.*pow(cc*h,2))/4. - 2.*cc*h*exp((cc*h)/2.) + 2.)/(pow(cc,3)*pow(h,2));
          mesh->MRSAAB_B[id+2] = -(cc*h - exp((cc*h)/2.) + (3.*pow(cc*h,2))/8. - (cc*h*exp((cc*h)/2.))/2. + 1.)/(pow(cc,3)*pow(h,2));
        }
        else{
          // Full dt coeeficients
          mesh->MRSAAB_A[id+0] = (pow(cc,3)*pow(h,4))/18. + (19.*pow(cc,2)*pow(h,3))/80. + (19.*cc*pow(h,2))/24. + (23.*h)/12.;
          mesh->MRSAAB_A[id+1] = - (4*h)/3. - (5.*cc*pow(h,2))/12 - (pow(cc,2)*pow(h,3))/10. - (7.*pow(cc,3)*pow(h,4))/360.;
          mesh->MRSAAB_A[id+2] = (pow(cc,3)*pow(h,4))/180. + (7.*pow(cc,2)*pow(h,3))/240. + (cc*pow(h,2))/8. + (5.*h)/12.;
          // Half coefficients
          mesh->MRSAAB_B[id+0] = (139.*pow(cc,3)*pow(h,4))/46080. + (pow(cc,2)*pow(h,3))/40. + (61.*cc*pow(h,2))/384. + (17.*h)/24.; 
          mesh->MRSAAB_B[id+1] = - (7.*h)/24. - (3.*cc*pow(h,2))/64. - (11.*pow(cc,2)*pow(h,3))/1920. - (13.*pow(cc,3)*pow(h,4))/23040. ; 
          mesh->MRSAAB_B[id+2] = (7.*pow(cc,3)*pow(h,4))/46080. + (pow(cc,2)*pow(h,3))/640. + (5.*cc*pow(h,2))/384. + h/12.;
        }
      // MRAB coefficients
        mesh->MRAB_A[id+0]   =  23.*h/12. ;
        mesh->MRAB_A[id+1]   = -16.*h/12. ;
        mesh->MRAB_A[id+2]   =  5. *h/12. ;

        mesh->MRAB_B[id+0]   =  17.*h/24. ;
        mesh->MRAB_B[id+1]   = - 7.*h/24. ;
        mesh->MRAB_B[id+2]   =   2.*h/24. ;
      }

    }

    // Exponential part
    mesh->MRSAAB_C[l]    = exp(cc*h);
    mesh->MRAB_C[l]      =   h ;
   }

  }


 //SetupOcca
  char deviceConfig[BUFSIZ];
  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", rank%2);
  //sprintf(deviceConfig, "mode = OpenCL, deviceID = 1, platformID = 0");
  //sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");  

  occa::kernelInfo kernelInfo;
  meshOccaSetup2D(mesh, deviceConfig,  kernelInfo);     

  // Setup MRAB PML
  boltzmannMRABPmlSetup2D(mesh, options); 



  mesh->o_MRABelementIds = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
  mesh->o_MRABhaloIds = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
  for (iint lev=0;lev<mesh->MRABNlevels;lev++) {
    if (mesh->MRABNelements[lev])
      mesh->o_MRABelementIds[lev] = mesh->device.malloc(mesh->MRABNelements[lev]*sizeof(iint),
         mesh->MRABelementIds[lev]);
    if (mesh->MRABNhaloElements[lev])
      mesh->o_MRABhaloIds[lev] = mesh->device.malloc(mesh->MRABNhaloElements[lev]*sizeof(iint),
         mesh->MRABhaloIds[lev]);
  }

  
 

  // Set AB 
  mesh->o_rhsq.free();
  mesh->o_rhsq = mesh->device.malloc(3*mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);

  if (mesh->totalHaloPairs) {
    //reallocate halo buffer for trace exchange
    mesh->o_haloBuffer.free();
    mesh->o_haloBuffer = mesh->device.malloc(mesh->totalHaloPairs*mesh->Nfp*mesh->Nfields*mesh->Nfaces*sizeof(dfloat));
  }
  


  mesh->o_c2 = mesh->device.malloc(mesh->Nelements*mesh->cubNp*sizeof(dfloat), mesh->c2);
  mesh->o_fQM = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*mesh->Nfields*sizeof(dfloat),
         mesh->fQM);
  mesh->o_fQP = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*mesh->Nfields*sizeof(dfloat),
         mesh->fQP);
  mesh->o_mapP = mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint), 
         mesh->mapP);


  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  int maxCubNodes = mymax(maxNodes,mesh->cubNp);

  kernelInfo.addDefine("p_maxNodes", maxNodes);
  kernelInfo.addDefine("p_maxCubNodes", maxCubNodes);


  int NblockV = 512/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  int NblockCub = 512/mesh->cubNp; // works for CUDA
  kernelInfo.addDefine("p_NblockCub", NblockCub);

  // physics 
  kernelInfo.addDefine("p_Lambda2", 0.5f);
  kernelInfo.addDefine("p_sqrtRT", mesh->sqrtRT);
  kernelInfo.addDefine("p_sqrt2", (float)sqrt(2.));
  kernelInfo.addDefine("p_isq12", (float)sqrt(1./12.));
  kernelInfo.addDefine("p_isq6", (float)sqrt(1./6.));
  kernelInfo.addDefine("p_invsqrt2", (float)sqrt(1./2.));
  kernelInfo.addDefine("p_tauInv", mesh->tauInv);




  kernelInfo.addDefine("p_q1bar", q1bar);
  kernelInfo.addDefine("p_q2bar", q2bar);
  kernelInfo.addDefine("p_q3bar", q3bar);
  kernelInfo.addDefine("p_q4bar", q4bar);
  kernelInfo.addDefine("p_q5bar", q5bar);
  kernelInfo.addDefine("p_q6bar", q6bar);
  kernelInfo.addDefine("p_alpha0", (float).01f);
  kernelInfo.addDefine("p_pmlAlpha", (float).2);


   
  // Volume and Relaxation Kernels
  if(strstr(options, "CUBATURE")){ 
          
    printf("Compiling volume kernel for cubature integration\n");
      mesh->volumeKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
         "boltzmannMRABVolumeCub2D",
          kernelInfo);

    printf("Compiling PML volume kernel for cubature integration\n");
      mesh->pmlVolumeKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
        "boltzmannMRABPmlVolumeCub2D",
          kernelInfo);

    if(strstr(options, "MRAB")){ 

      printf("Compiling MRAB relaxation kernel with cubature integration\n");
       mesh->relaxationKernel =
         mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannRelaxation2D.okl",
          "boltzmannMRABRelaxationCub2D",
           kernelInfo); 
    
      //
      printf("Compiling MRAB PML relaxation kernel with cubature integration\n");
        mesh->pmlRelaxationKernel =
         mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannRelaxation2D.okl",
          "boltzmannMRABPmlRelaxationCub2D",
            kernelInfo);  
    }
    else if(strstr(options, "MRSAAB")){

    printf("Compiling MRSAAB relaxation kernel with cubature integration\n");
     mesh->relaxationKernel =
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannRelaxation2D.okl",
        "boltzmannMRSAABRelaxationCub2D",
        kernelInfo); 

      //
    printf("Compiling MRSAAB PML relaxation kernel with cubature integration\n");
    mesh->pmlRelaxationKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannRelaxation2D.okl",
      "boltzmannMRSAABPmlRelaxationCub2D",
       kernelInfo);  
    }
  }



  if(strstr(options, "COLLOCATION")){ 

     if(strstr(options, "MRAB")){ 

      printf("Compiling MRAB volume kernel with nodal collocation for nonlinear term\n");
      mesh->volumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
             "boltzmannABVolume2D",
             kernelInfo);

      printf("Compiling MRAB PML volume kernel with nodal collocation for nonlinear term\n");
      mesh->pmlVolumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
           "boltzmannABPmlVolume2D",
           kernelInfo); 
    }

    else if(strstr(options, "MRSAAB")){
      printf("Compiling MRSAAB volume kernel with nodal collocation for nonlinear term\n");
      mesh->volumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
             "boltzmannMRSAABVolume2D",
             kernelInfo);

      printf("Compiling MRSAAB PML volume kernel with nodal collocation for nonlinear term\n");
      mesh->pmlVolumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
           "boltzmannMRSAABPmlVolume2D",
           kernelInfo); 
      }


  }

 

   // UPDATE Kernels
  if(strstr(options, "MRAB")){ 
    printf("Compiling MRAB update kernel\n");
    mesh->updateKernel =
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannMRABUpdate2D",
            kernelInfo);
    
    printf("Compiling MRAB trace update kernel\n");
    mesh->traceUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
               "boltzmannMRABTraceUpdate2D",
               kernelInfo);
    
    printf("Compiling MRAB PML update kernel\n");
    mesh->pmlUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannMRABPmlUpdate2D",
            kernelInfo);
   
    printf("Compiling MRAB PML trace update kernel\n");
    mesh->pmlTraceUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
               "boltzmannMRABPmlTraceUpdate2D",
                 kernelInfo);
    }

    else if(strstr(options, "MRSAAB")){

    printf("Compiling MRSAAB update kernel\n");
    mesh->updateKernel =
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannMRSAABUpdate2D",
            kernelInfo);

    printf("Compiling MRSAAB trace update kernel\n");
    mesh->traceUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
               "boltzmannMRSAABTraceUpdate2D",
               kernelInfo);
    
     printf("Compiling MRSAAB PML update kernel\n");
    mesh->pmlUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannMRSAABPmlUpdate2D",
            kernelInfo);
    
    printf("Compiling MRSAAB PML trace update kernel\n");
    mesh->pmlTraceUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
               "boltzmannMRSAABPmlTraceUpdate2D",
                 kernelInfo);
     }


    // Surface Kernels
    printf("Compiling surface kernel\n");
    mesh->surfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
         "boltzmannMRABSurface2D",
         kernelInfo);

    printf("Compiling PML surface kernel\n");
    mesh->pmlSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
       "boltzmannMRABPmlSurface2D",
       kernelInfo);
      

    mesh->haloExtractKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
          "meshHaloExtract2D",
             kernelInfo);



  
}




