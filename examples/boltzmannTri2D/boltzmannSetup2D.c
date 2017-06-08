#include "boltzmann2D.h"
//#include <complex.h>

// NBN: toggle use of 2nd stream
#define USE_2_STREAMS

void boltzmannSetup2D(mesh2D *mesh, char * options){

  mesh->Nfields = 8;
  
  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
               sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
               sizeof(dfloat));
  mesh->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
               sizeof(dfloat));

  mesh->pmlqx    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
            sizeof(dfloat));
  mesh->rhspmlqx = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
            sizeof(dfloat));
  mesh->respmlqx = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
            sizeof(dfloat));

  mesh->pmlqy    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
            sizeof(dfloat));
  mesh->rhspmlqy = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
            sizeof(dfloat));
  mesh->respmlqy = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
            sizeof(dfloat));

  
  
  mesh->sigmax = (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
  mesh->sigmay = (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
   
  // SET SOLVER OPTIONS 
  // Initial Conditions, Flow Properties
  dfloat Ma = 0.f , Re= 0.f, rho = 1.f, u = 0.f, v = 0.f, nu = 0.f, Uref=1.f, Lref = 1.f; 
  dfloat sigma11 = 0.f , sigma12 = 0.f, sigma22 = 0.f;  
  if(strstr(options, "PML")){
    printf("Starting initial conditions for PML\n");
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
    
    //printf("starting initial conditions\n"); //Zero Flow Conditions
    rho = 1., u = Uref; v = 0.; sigma11 = 0, sigma12 = 0, sigma22 = 0;
    //
    mesh->finalTime = 20.0;
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


     if(strstr(options, "PML")){
      iint id = mesh->Np*mesh->Nfields*e + n;
      mesh->pmlqx[id+0*mesh->Np] = 0.f*q1bar;
      mesh->pmlqx[id+1*mesh->Np] = 0.f*q2bar;
      mesh->pmlqx[id+2*mesh->Np] = 0.f*q3bar;
      mesh->pmlqx[id+3*mesh->Np] = 0.f*q4bar;
      mesh->pmlqx[id+4*mesh->Np] = 0.f*q5bar;
      mesh->pmlqx[id+5*mesh->Np] = 0.f*q6bar;

      mesh->pmlqy[id+0*mesh->Np] = 0.f*q1bar;
      mesh->pmlqy[id+1*mesh->Np] = 0.f*q2bar;
      mesh->pmlqy[id+2*mesh->Np] = 0.f*q3bar;
      mesh->pmlqy[id+3*mesh->Np] = 0.f*q4bar;
      mesh->pmlqy[id+4*mesh->Np] = 0.f*q5bar;
      mesh->pmlqy[id+5*mesh->Np] = 0.f*q6bar;
    }

    }
  }

  //
  

  printf("starting parameters\n");
  
  // set penalty parameter
  mesh->Lambda2 = 0.5/(mesh->sqrtRT);

  

  dfloat xmin = -4, xmax = 8, ymin = -4, ymax = 4;

  dfloat xsigma  = 80, ysigma  = 80; // For Quadratic Pml Profile
  dfloat cxsigma = 80, cysigma = 80; // For Constant  Pml Profile

    
  iint *pmlElementIds = (iint*) calloc(mesh->Nelements, sizeof(iint));
  iint *nonPmlElementIds = (iint*) calloc(mesh->Nelements, sizeof(iint));
  iint pmlNelements = 0;
  iint nonPmlNelements = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    dfloat cx = 0, cy = 0;
    for(iint n=0;n<mesh->Nverts;++n){
      cx += mesh->EX[e*mesh->Nverts+n];
      cy += mesh->EY[e*mesh->Nverts+n];
    }
    cx /= mesh->Nverts;
    cy /= mesh->Nverts;
    
   iint isPml = 0;
    
for(iint n=0;n<mesh->Np;++n){
  dfloat x = mesh->x[n + e*mesh->Np];
  dfloat y = mesh->y[n + e*mesh->Np];
      //      if(cx<xmax+1 && cx>xmin-1 && cy<ymax+1 && cy>ymin-1){

  if(strstr(options,"PML")){

    if(strstr(options,"QUADRATIC")){

      if(cx>xmax){
        mesh->sigmax[mesh->Np*e + n] = xsigma*pow(x-xmax,2);
        isPml = 1;
      }
      if(cx<xmin){
         mesh->sigmax[mesh->Np*e + n] = xsigma*pow(x-xmin,2);
         isPml = 1;
      }
      if(cy>ymax){
        mesh->sigmay[mesh->Np*e + n] = ysigma*pow(y-ymax,2);
        isPml = 1;
      }
      if(cy<ymin){
        mesh->sigmay[mesh->Np*e + n] = ysigma*pow(y-ymin,2);   
        isPml = 1;
      }
    
    }

    if(strstr(options,"CONSTANT")){

      if(cx>xmax){
        mesh->sigmax[mesh->Np*e + n] = cxsigma;
        isPml = 1;
      }
      if(cx<xmin){
         mesh->sigmax[mesh->Np*e + n]= cxsigma;
        isPml = 1;
      }
      if(cy>ymax){
        mesh->sigmay[mesh->Np*e + n] = cysigma;
        isPml = 1;
      }
      if(cy<ymin){
        mesh->sigmay[mesh->Np*e + n] = cysigma;
        isPml = 1;
      }
    }
    }
  }
    
    if(isPml)
      pmlElementIds[pmlNelements++] = e;
    else
      nonPmlElementIds[nonPmlNelements++] = e;
    
  }



  printf("detected pml: %d pml %d non-pml %d total \n", pmlNelements, nonPmlNelements, mesh->Nelements);

  // set time step
  dfloat hmin = 1e9, hmax = 0;
  for(iint e=0;e<mesh->Nelements;++e){ 

    //printf("global index: %d  pml = %d and Nonpml= %d\n",e, pmlElementIds[e], nonPmlElementIds[e]); 

    for(iint f=0;f<mesh->Nfaces;++f){
      iint sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];
     
      dfloat hest = 2./(sJ*invJ);

      hmin = mymin(hmin, hest);
      hmax = mymax(hmax, hest);
    }
  }


  
  dfloat cfl = 0.5; 
  dfloat magVelocity = sqrt(q2bar*q2bar+q3bar*q3bar)/(q1bar/mesh->sqrtRT);
  magVelocity = mymax(magVelocity,1.0); // Correction for initial zero velocity
  //
  dfloat dtex = hmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*mesh->sqrtRT);
  dfloat dtim = 1./(mesh->tauInv);

  dfloat dt = 0.f;

  // AK: Set time step size
  if(strstr(options, "LSERK")){
    printf("Time discretization method: LSERK with CFL: %.2f \n",cfl);
    dt = cfl*mymin(dtex,dtim);
    printf("dt = %.4e explicit-dt = %.4e , implicit-dt= %.4e  ratio= %.4e\n", dt,dtex,dtim, dtex/dtim);
 
  }
  if(strstr(options, "SARK3")){ 
    printf("Time discretization method: SARK with CFL: %.2f \n",cfl);
    dt = cfl*dtex;
    printf("dt = %.4e explicit-dt = %.4e , implicit-dt= %.4e  ratio= %.4e\n", dt,dtex,dtim, dtex/dtim);

  }
  if(strstr(options, "SAAB3")){ 
    printf("Time discretization method: Low Storage SAAB  with CFL: %.2f \n",cfl);
    dt = 1.0/3.0 *cfl*dtex; 
    printf("dt = %.4e explicit-dt = %.4e , implicit-dt= %.4e  ratio= %.4e\n", dt,dtex,dtim, dtex/dtim);
  }
  if(strstr(options, "LSIMEX")){ 
    printf("Time discretization method: Low Storage IMEX  with CFL: %.2f \n",cfl);
    dt = cfl*dtex;     
    printf("dt = %.4e explicit-dt = %.4e , implicit-dt= %.4e  ratio= %.4e\n", dt,dtex,dtim, dtex/dtim);
  }



  printf("hmin = %g\n", hmin);
  printf("hmax = %g\n", hmax);
  printf("cfl = %g\n", cfl);
  printf("dt = %g ", dt);
  printf("max wave speed = %g\n", sqrt(3.)*mesh->sqrtRT);
  
  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
 
  
  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  mesh->dt = mesh->finalTime/mesh->NtimeSteps;

  // errorStep
  mesh->errorStep = 5000;

  printf("Nsteps = %d with dt = %.8e\n", mesh->NtimeSteps, mesh->dt);

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // use rank to choose DEVICE
 // sprintf(deviceConfig, "mode = CUDA, deviceID = %d", (rank+1)%3);
  //sprintf(deviceConfig, "mode = OpenCL, deviceID = 1, platformID = 0");
  // sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 0");

   sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");  


  occa::kernelInfo kernelInfo;

  meshOccaSetup2D(mesh, deviceConfig,  kernelInfo);




  if(strstr(options, "LSERK")){

   if(strstr(options, "PML")){ 
    mesh->o_pmlqx =    
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);
    mesh->o_rhspmlqx =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqx);
    mesh->o_respmlqx =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->respmlqx);

    mesh->o_pmlqy =    
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
    mesh->o_rhspmlqy =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqy);
    mesh->o_respmlqy =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->respmlqy);
      }  

  }
  else if(strstr(options, "SARK3")){ 
    mesh->o_qold =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);

    mesh->o_rhsq =
        mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);    
    mesh->o_rhsq2 =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);
    mesh->o_rhsq3 =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);

       
    if(strstr(options, "PML")){
    // pml variables
    mesh->o_pmlqx =    
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);
    // pml variables
    mesh->o_pmlqxold =    
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);

    mesh->o_rhspmlqx =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqx);
   
    mesh->o_rhspmlqx2 =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->respmlqx);
    
    mesh->o_rhspmlqx3 =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->respmlqx);

    // pml variables
    // pml variables
    mesh->o_pmlqy =    
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);
    // pml variables
    mesh->o_pmlqyold =    
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);

    mesh->o_rhspmlqy =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqx);
   
    mesh->o_rhspmlqy2 =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->respmlqx);
    
    mesh->o_rhspmlqy3 =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->respmlqx);
    }

    //
    for(int i=0; i<5; i++){
      for(int j=0; j<5; j++){
  mesh->sarka[i][j] = 0.0;
  //mesh->sarkpmla[i][j] = 0.0 ; 
      }
      mesh->sarkb[i] = 0.0;
      mesh->sarke[i] = 0.0;
      // mesh->sarkpmlb[i] = 0.0;
      // mesh->sarkpmle[i] = 0.0;     
    }

    dfloat coef = -mesh->tauInv;
    dfloat  h   = mesh->dt; 



#if 0 // First version of SARK
    // Base Runge-Kutta Method
    dfloat a21 = 1./3.;   dfloat a31 = -3./16. ;    dfloat a32 = 15./16.;
    dfloat b1 = 1./6.;    dfloat b2 = 3./10.;       dfloat b3 = 8./15.; 
    dfloat c1 = 0.;       dfloat c2 = 1./3.;        dfloat c3 = 3./4.; 
      
#else // Second Version
    // Base Runge-Kutta Method
    dfloat a21 = 1.f/2.f;   dfloat a31 = -1.f ;    dfloat a32 = 2.f;
    dfloat b1 = 1.f/6.f;    dfloat b2 = 2./3.;       dfloat b3 = 1./6.; 
    dfloat c1 = 0.f;       dfloat c2 = 1./2.;        dfloat c3 = 1.; 
    
#endif

    // Base Method
    mesh->rk3a[0][0] = 0.;  mesh->rk3a[1][0] = a21;  mesh->rk3a[2][0] = a31;   mesh->rk3a[2][1] = a32; 
    mesh->rk3b[0] = b1;     mesh->rk3b[1] = b2;      mesh->rk3b[2] = b3; 
    mesh->rk3c[0] = c1;     mesh->rk3c[1] = c2;      mesh->rk3c[2] = c3; 


    printf("Coef*dt = %.8e", coef*dt);

    if(fabs(coef*h)>1e-2){
#if 0  // First Version

      //  Exponential Coefficients
      mesh->sarka[1][0] = -(a21*exp(c2*coef*h)*(exp(-c2*coef*h) - 1.))/(c2*coef*h); // a21
      mesh->sarka[2][0] = -(a31*exp(c3*coef*h)*(exp(-c3*coef*h) - 1.))/(c3*coef*h); // a31
      mesh->sarka[2][1] = -(a32*exp(c3*coef*h)*(exp(-c3*coef*h) - 1.))/(c3*coef*h); // a32 

      // If 1/tau*h is too small say <<1, need to write in terms of Taylor coefficients
      mesh->sarkb[0] =   (exp(coef*h)*((exp(-coef*h)*(c2 + c3 - c2*c3 + c2*c3*exp(coef*h) - 1.))
               /(coef*(c1 - c2)*(c1 - c3)) + (exp(-coef*h)*(2.*exp(coef*h) - 2.) - coef*h*exp(-coef*h)*(c2*exp(coef*h) - c3 - c2 + c3*exp(coef*h) + 2.))
               /(pow(coef,3.)*pow(h,2)*(c1 - c2)*(c1 - c3))))/h;
      mesh->sarkb[1] =  -(exp(coef*h)*((exp(-coef*h)*(c1 + c3 - c1*c3 + c1*c3*exp(coef*h) - 1.))
               /(coef*(c1 - c2)*(c2 - c3)) + (exp(-coef*h)*(2.*exp(coef*h) - 2.) - coef*h*exp(-coef*h)*(c1*exp(coef*h) - c3 - c1 + c3*exp(coef*h) + 2.))
               /(pow(coef,3)*pow(h,2)*(c1 - c2)*(c2 - c3))))/h;
      mesh->sarkb[2] =   (exp(coef*h)*((exp(-coef*h)*(c1 + c2 - c1*c2 + c1*c2*exp(coef*h) - 1.))/(coef*(c1 - c3)*(c2 - c3)) + (exp(-coef*h)*(2.*exp(coef*h) - 2.) - coef*h*exp(-coef*h)*(c1*exp(coef*h) - c2 - c1 + c2*exp(coef*h) + 2.))
               /(pow(coef,3)*pow(h,2)*(c1 - c3)*(c2 - c3))))/h;
      //
      //
      mesh->sarke[0] = exp(coef*h*c2); 
      mesh->sarke[1] = exp(coef*h*c3); 
      mesh->sarke[2] = exp(coef*h*1.0);
#else //Second formulation 

      //  Exponential Coefficients
      mesh->sarka[1][0] = (exp(coef*h/2.) - 1.)/(coef*h); // a21
      mesh->sarka[2][0] = -1.0*(exp(coef*h)-1.0)/(coef*h); // a31
      mesh->sarka[2][1] =  2.0*(exp(coef*h)-1.0)/(coef*h);// a32 

      // If 1/tau*h is too small say <<1, need to write in terms of Taylor coefficients
      mesh->sarkb[0] =   (-4. -coef*h + exp(coef*h)*(4.-3.*coef*h+pow(h*coef,2.)))/ (pow(coef*h,3)) ;
      mesh->sarkb[1] =  4.*(2. + coef*h + exp(coef*h)*(-2. + coef*h)) / (pow(coef*h,3)) ;
      mesh->sarkb[2] =   (-4. -3.*coef*h - pow(coef*h,2)+ exp(coef*h)*(4. - coef*h))/ (pow(coef*h,3)) ;
      //
      //
      mesh->sarke[0] = exp(coef*h*c2); 
      mesh->sarke[1] = exp(coef*h*c3); 
      mesh->sarke[2] = exp(coef*h*1.0);
#endif
       

    }
    else{

      printf("Computing SARK coefficients  with 3th order Taylor series expansion\n");

#if 0 // First formulation

      //  fifth Order Taylor Series Expansion
      mesh->sarka[1][0] = (a21*c2*(pow(c2,4)*pow(coef,4)*pow(h,4) + 5.*pow(c2,3)*pow(coef,3)*pow(h,3) 
           + 20.*pow(c2,2)*pow(coef,2)*pow(h,2) + 60.*c2*coef*h + 120.))/120. ;  // a21
      mesh->sarka[2][0] = (a31*c3*(pow(c3,4)*pow(coef,4)*pow(h,4) + 5.*pow(c3,3)*pow(coef,3)*pow(h,3) 
           + 20.*pow(c3,2)*pow(coef,2)*pow(h,2) + 60.*c3*coef*h + 120.))/120. ; // a31
      mesh->sarka[2][1] = (a32*c3*(pow(c3,4)*pow(coef,4)*pow(h,4) + 5.*pow(c3,3)*pow(coef,3)*pow(h,3) 
           + 20.*pow(c3,2)*pow(coef,2)*pow(h,2) + 60.*c3*coef*h + 120.))/120. ; // a32 

      // If 1/tau*h is too small say <<1, need to write in terms of Taylor coefficients
      mesh->sarkb[0] =   -(2520.*c2 + 2520.*c3 - 5040.*c2*c3 - 420.*coef*h - 84.*pow(coef,2)*pow(h,2) - 14.*pow(coef,3)*pow(h,3) 
         - 2.*pow(coef,4)*pow(h,4) + 210.*c2*pow(coef,2)*pow(h,2) + 210.*c3*pow(coef,2)*pow(h,2) + 42.*c2*pow(coef,3)*pow(h,3) + 42.*c3*pow(coef,3)*pow(h,3) 
         + 7.*c2*pow(coef,4)*pow(h,4) + 7.*c3*pow(coef,4)*pow(h,4) + 840.*c2*coef*h + 840.*c3*coef*h - 840.*c2*c3*pow(coef,2)*pow(h,2) 
         - 210.*c2*c3*pow(coef,3)*pow(h,3) - 42.*c2*c3*pow(coef,4)*pow(h,4) - 2520.*c2*c3*coef*h - 1680.)/(5040.*(c1 - c2)*(c1 - c3));

      mesh->sarkb[1] =  (2520.*c1 + 2520.*c3 - 5040.*c1*c3 - 420.*coef*h - 84.*pow(coef,2)*pow(h,2) - 14.*pow(coef,3)*pow(h,3) - 2*pow(coef,4)*pow(h,4) 
       + 210*c1*pow(coef,2)*pow(h,2) + 42*c1*pow(coef,3)*pow(h,3) + 210*c3*pow(coef,2)*pow(h,2) + 7.*c1*pow(coef,4)*pow(h,4) + 42.*c3*pow(coef,3)*pow(h,3) 
       + 7.*c3*pow(coef,4)*pow(h,4) + 840.*c1*coef*h + 840.*c3*coef*h - 840.*c1*c3*pow(coef,2)*pow(h,2) - 210.*c1*c3*pow(coef,3)*pow(h,3) 
       - 42.*c1*c3*pow(coef,4)*pow(h,4) - 2520.*c1*c3*coef*h - 1680.)/(5040.*(c1 - c2)*(c2 - c3));

      mesh->sarkb[2] =   -(2520.*c1 + 2520.*c2 - 5040.*c1*c2 - 420.*coef*h - 84.*pow(coef,2)*pow(h,2) - 14.*pow(coef,3)*pow(h,3) - 2.*pow(coef,4)*pow(h,4) 
         + 210.*c1*pow(coef,2)*pow(h,2) + 210.*c2*pow(coef,2)*pow(h,2) + 42.*c1*pow(coef,3)*pow(h,3) + 42.*c2*pow(coef,3)*pow(h,3) + 7.*c1*pow(coef,4)*pow(h,4) 
         + 7.*c2*pow(coef,4)*pow(h,4) + 840.*c1*coef*h + 840.*c2*coef*h - 840.*c1*c2*pow(coef,2)*pow(h,2) - 210.*c1*c2*pow(coef,3)*pow(h,3) 
         - 42.*c1*c2*pow(coef,4)*pow(h,4) - 2520.*c1*c2*coef*h - 1680.)/(5040.*(c1 - c3)*(c2 - c3));

      mesh->sarke[0] = exp(coef*h*c2); 
      mesh->sarke[1] = exp(coef*h*c3); 
      mesh->sarke[2] = exp(coef*h*1.0);           
#else  // Second Formulation

      //  Exponential Coefficients
      mesh->sarka[1][0] = (pow(coef,4)*pow(h,4))/3840. + (pow(coef,3)*pow(h,3))/384. + (pow(coef,2)*pow(h,2))/48. + (coef*h)/8 + 1./2.;
      mesh->sarka[2][0] =  -1.0 *  ((pow(coef,4)*pow(h,4))/120. + (pow(coef,3)*pow(h,3))/24. + (pow(coef,2)*pow(h,2))/6. + (coef*h)/2. + 1.);
      mesh->sarka[2][1] =  2.0*((pow(coef,4)*pow(h,4))/120. + (pow(coef,3)*pow(h,3))/24. + (pow(coef,2)*pow(h,2))/6. + (coef*h)/2. + 1.);

      mesh->sarkb[0] =  (5.*pow(coef,4)*pow(h,4))/1008. + (pow(coef,3)*pow(h,3))/45. + (3.*pow(coef,2)*pow(h,2))/40. + (coef*h)/6. + 1./6. ; 
      mesh->sarkb[1] =  (pow(coef,4)*pow(h,4))/252. + (pow(coef,3)*pow(h,3))/45. + (pow(coef,2)*pow(h,2))/10. + (coef*h)/3. + 2./3. ;
      mesh->sarkb[2] =  1./6. - (pow(coef,3)*pow(h,3))/360. - (pow(coef,4)*pow(h,4))/1680. - (pow(coef,2)*pow(h,2))/120. ;

      //
      mesh->sarke[0] = exp(coef*h*c2); 
      mesh->sarke[1] = exp(coef*h*c3); 
      mesh->sarke[2] = exp(coef*h*1.0);


#endif                



    }
         
  }
  else if(strstr(options, "SAAB3")){ 
    mesh->o_rhsq2 =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);
    mesh->o_rhsq3 =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);

    mesh->o_pmlqx =    
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);
    mesh->o_rhspmlqx =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqx);
    mesh->o_rhspmlqx2 =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqx);
    mesh->o_rhspmlqx3 =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqx);  

    mesh->o_pmlqy =    
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
    mesh->o_rhspmlqy =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqy);
    mesh->o_rhspmlqy2 =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqy);
    mesh->o_rhspmlqy3 =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqy);  


    // Classical Adams Bashforth Coefficients
    mesh->mrab[0] = 23.*dt/12. ;
    mesh->mrab[1] = -4.*dt/3. ;
    mesh->mrab[2] =  5.*dt/12. ;

    // SAAB NONPML Coefficients, expanded to fix very small tauInv case
    dfloat cc = -mesh->tauInv;
    dfloat h  = mesh->dt; 
    //
    if(fabs(cc*h)>1e-2){
      mesh->saab[0] = (exp(cc*h) - (5.*cc*h)/2. - 3.*pow(cc,2)*pow(h,2 )+ pow(cc,2)*pow(h,2)*exp(cc*h) + (3.*cc*h*exp(cc*h))/2. - 1.)/(pow(cc,3)*pow(h,2));
      mesh->saab[1] = (4.*cc*h - 2.*exp(cc*h) + 3.*pow(cc,2)*pow(h,2 )- 2.*cc*h*exp(cc*h) + 2.)/(pow(cc,3)*pow(h,2));
      mesh->saab[2] = -((3.*cc*h)/2. - exp(cc*h) + pow(cc,2)*pow(h,2 )- (cc*h*exp(cc*h))/2. + 1.)/(pow(cc,3)*pow(h,2));
      mesh->saabexp = exp(cc*h);
    }
    else{

      mesh->saab[0] = (pow(cc,3)*pow(h,4))/18. + (19.*pow(cc,2)*pow(h,3))/80. + (19.*cc*pow(h,2))/24. + (23.*h)/12.;
      mesh->saab[1] = -(4.*h)/3. - (5.*cc*pow(h,2))/12. - (pow(cc,2)*pow(h,3))/10. - (7.*pow(cc,3)*pow(h,4))/360.;
      mesh->saab[2] = (pow(cc,3)*pow(h,4))/180. + (7.*pow(cc,2)*pow(h,3))/240. + (cc*pow(h,2))/8. + (5.*h)/12.;
      mesh->saabexp = exp(cc*h);
    }         
  }


  else if(strstr(options, "LSIMEX")){ 
    mesh->o_qY =    
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
    // pml variables
    mesh->o_qZ =    
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
    // 
    mesh->o_qS =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);


    if(strstr(options, "PML")){ 
      // pml variables
      mesh->o_pmlqx =    
  mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);
      mesh->o_qSx =
  mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqx);
      mesh->o_qYx =
  mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqx);
      mesh->o_qZx =
  mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->respmlqx);

      mesh->o_pmlqy =    
  mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_qSy =
  mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqx);
      mesh->o_qYy =
  mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqy);
      mesh->o_qZy =
  mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->respmlqy);
    }
  } 
  
  



  mesh->o_sigmax =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), mesh->sigmax);

  mesh->o_sigmay =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), mesh->sigmay);

  mesh->nonPmlNelements = nonPmlNelements;
  mesh->pmlNelements    = pmlNelements;

  if(mesh->nonPmlNelements)
    mesh->o_nonPmlElementIds = 
      mesh->device.malloc(nonPmlNelements*sizeof(iint), nonPmlElementIds);

  if(mesh->pmlNelements)
    mesh->o_pmlElementIds = 
      mesh->device.malloc(pmlNelements*sizeof(iint), pmlElementIds);

  // specialization for Boltzmann

  kernelInfo.addDefine("p_maxNodesVolume", mymax(mesh->cubNp,mesh->Np));
  
  kernelInfo.addDefine("p_pmlAlpha", (float).2);
  
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 128/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 128/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  printf("Np: %d \t Ncub: %d \n", mesh->Np, mesh->cubNp);

  printf("maxNodes: %d \t NblockV: %d \t NblockS: %d  \n", maxNodes, NblockV, NblockS);

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

  
  if(strstr(options, "LSERK")){ 

    if(strstr(options, "CUBATURE")){ 
          
      printf("Compiling LSERK volume kernel with cubature integration\n");
      mesh->volumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
           "boltzmannVolumeCub2D",
           kernelInfo);

      printf("Compiling LSERK relaxation kernel with cubature integration\n");
       mesh->relaxationKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannRelaxation2D.okl",
          "boltzmannRelaxationCub2D",
          kernelInfo); 
      printf("Compiling LSERK Unsplit pml volume kernel with cubature integration\n");
      mesh->pmlVolumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
          "boltzmannPmlVolumeCub2D",
          kernelInfo);
      //
      printf("Compiling LSERK Unsplit pml relaxation kernel with cubature integration\n");
      mesh->pmlRelaxationKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannRelaxation2D.okl",
      "boltzmannPmlRelaxationCub2D",
      kernelInfo);  
    }

    if(strstr(options, "COLLOCATION")){ 
      printf("Compiling pml volume kernel with nodal collocation for nonlinear term\n");
      mesh->volumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
             "boltzmannVolume2D",
             kernelInfo);

      printf("Compiling Unsplit pml volume kernel with nodal collocation for nonlinear term\n");
      mesh->pmlVolumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
           "boltzmannPmlVolume2D",
           kernelInfo); 
    }


    printf("Compiling surface kernel\n");
    mesh->surfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
         "boltzmannSurface2D",
         kernelInfo);

    printf("Compiling Unsplit  pml surface kernel\n");
    mesh->pmlSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
       "boltzmannPmlSurface2D",
       kernelInfo);
      

    printf("Compiling update kernel\n");
    mesh->updateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannLSERKUpdate2D",
         kernelInfo);

        
    printf("Compiling Unsplit pml update kernel\n");
    mesh->pmlUpdateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
     "boltzmannLSERKPmlUpdate2D",
     kernelInfo);
      
    mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
    "meshHaloExtract2D",
    kernelInfo);

  }
  
  else if(strstr(options, "SARK3")){ 
    if(strstr(options, "CUBATURE")){ 
      printf("Compiling SA volume kernel with cubature integration\n");
      mesh->volumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
               "boltzmannVolumeCub2D",
               kernelInfo);

      printf("Compiling SA relaxation kernel with cubature integration\n");
      mesh->relaxationKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannRelaxation2D.okl",
           "boltzmannSARelaxationCub2D",
           kernelInfo); 

        printf("Compiling SA Unsplit pml volume kernel with cubature integration\n");
        mesh->pmlVolumeKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
           "boltzmannPmlVolumeCub2D",
           kernelInfo);

        printf("Compiling SA Unsplit pml relaxation kernel with cubature integration\n");
        mesh->pmlRelaxationKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannRelaxation2D.okl",
           "boltzmannSAPmlRelaxationCub2D",
           kernelInfo); 
  }
  else if(strstr(options, "COLLOCATION")){ 
    printf("Compiling SA volume kernel\n");
    mesh->volumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
           "boltzmannSAVolume2D",
           kernelInfo);

      printf("Compiling SA Unsplit pml volume kernel\n");
      mesh->pmlVolumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
             "boltzmannSAPmlVolume2D",
             kernelInfo); 
  }
  //
  printf("Compiling surface kernel\n");
  mesh->surfaceKernel =
  mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
       "boltzmannSurface2D",
       kernelInfo);

      printf("Compiling Unsplit  pml surface kernel\n");
      mesh->pmlSurfaceKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
         "boltzmannPmlSurface2D",
         kernelInfo);
  
  //SARK STAGE UPDATE
  printf("compiling SARK non-pml  update kernel\n");
  mesh->updateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannSARK3Update2D",
         kernelInfo); 


  printf("compiling SARK non-pml  stage update kernel\n");
  mesh->updateStageKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
         "boltzmannSARK3StageUpdate2D",
         kernelInfo); 

  //
  printf("compiling SARK Unsplit pml  update kernel\n");
  mesh->pmlUpdateKernel =
  mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
       "boltzmannSARK3PmlUpdate2D",
       kernelInfo); 

  printf("compiling SARK Unsplit pml stage update kernel\n");
  mesh->pmlUpdateStageKernel =
  mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
  "boltzmannSARK3PmlStageUpdate2D",
  kernelInfo); 
     
  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
         "meshHaloExtract2D",
         kernelInfo);   


  }

  else if(strstr(options, "LSIMEX")){ 
    // RESIDUAL UPDATE KERNELS
    printf("Compiling LSIMEX non-pml residual update kernel\n");
    mesh->residualUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
           "boltzmannLSIMEXResidualUpdate2D",
           kernelInfo);
    
    printf("Compiling LSIMEX non-pml implicit update kernel\n");
    mesh->implicitUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
           "boltzmannLSIMEXImplicitUpdate2D",
           kernelInfo);

    printf("Compiling LSIMEX Unsplit pml residual update kernel\n");
   mesh->pmlResidualUpdateKernel =
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
          "boltzmannLSIMEXPmlResidualUpdate2D",
          kernelInfo);
     //
   printf("Compiling LSIMEX Unsplit pml implicit update kernel\n");
   mesh->pmlImplicitUpdateKernel =
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
          "boltzmannLSIMEXPmlImplicitUpdate2D",
          kernelInfo);
      
     
    if(strstr(options, "CUBATURE")){ 
      printf("Compiling LSIMEX non-pml Implicit Iteration Cubature  kernel\n");

       mesh->implicitSolveKernel = 
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannLSIMEXImplicitSolve2D.okl",
             "boltzmannLSIMEXImplicitSolveCub2D",
             kernelInfo); 

          printf("Compiling LSIMEX pml Implicit Iteration Cubature  kernel\n");
      mesh->pmlImplicitSolveKernel = 
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannLSIMEXImplicitSolve2D.okl",
             "boltzmannLSIMEXImplicitSolveCub2D",
             kernelInfo); 
    }
    else if(strstr(options, "COLLOCATION")){ 
      //
      printf("Compiling LSIMEX non-pml Implicit Iteration kernel\n");
      mesh->implicitSolveKernel = 
      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannLSIMEXImplicitSolve2D.okl",
        "boltzmannLSIMEXImplicitSolve2D",
        kernelInfo); 
        
     printf("Compiling LSIMEX Unsplit pml Implicit Iteration  kernel\n");
     mesh->pmlImplicitSolveKernel = 
     mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannLSIMEXImplicitSolve2D.okl",
     "boltzmannLSIMEXImplicitSolve2D",
     kernelInfo);      
      
    }

  printf("Compiling LSIMEX volume kernel integration\n");
    mesh->volumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
         "boltzmannVolumeCub2D",
         kernelInfo);

       
    printf("Compiling LSERK Unsplit pml volume kernel with cubature integration\n");
    mesh->pmlVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
        "boltzmannPmlVolumeCub2D",
        kernelInfo);
       
   printf("Compiling surface kernel\n");
    mesh->surfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
         "boltzmannSurface2D",
         kernelInfo);
        
    printf("Compiling Unsplit  pml surface kernel\n");
    mesh->pmlSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
       "boltzmannPmlSurface2D",
       kernelInfo);

    printf("Compiling LSIMEX non-pml update kernel\n");
    mesh->updateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
      "boltzmannLSIMEXUpdate2D",
      kernelInfo);
  //
    printf("Compiling LSIMEX Unsplit pml update kernel\n");
       mesh->pmlUpdateKernel =
       mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
       "boltzmannLSIMEXPmlUpdate2D",
       kernelInfo);
  
    mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
      "meshHaloExtract2D",
      kernelInfo);

  }


  else if(strstr(options, "SAAB3")){ 

  if(strstr(options, "CUBATURE")){ 
    printf("Compiling SA volume kernel with cubature integration\n");
    mesh->volumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
    "boltzmannVolumeCub2D",
    kernelInfo);

    printf("Compiling SA relaxation kernel with cubature integration\n");
    mesh->relaxationKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannRelaxation2D.okl",
    "boltzmannSARelaxationCub2D",
    kernelInfo); 

    printf("Compiling SA Unsplit pml volume kernel with cubature integration\n");
    mesh->pmlVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
     "boltzmannPmlVolumeCub2D",
     kernelInfo);

    printf("Compiling SA Unsplit pml relaxation kernel with cubature integration\n");
    mesh->pmlRelaxationKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannRelaxation2D.okl",
     "boltzmannSAPmlRelaxationCub2D",
     kernelInfo); 

  }
  else if(strstr(options, "COLLOCATION")){ 
    printf("Compiling SA volume kernel\n");
    mesh->volumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
     "boltzmannSAVolume2D",
     kernelInfo);

    printf("Compiling SA Unsplit pml volume kernel\n");
    mesh->pmlVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolume2D.okl",
     "boltzmannSAPmlVolume2D",
     kernelInfo); 
    }

  printf("Compiling surface kernel\n");
  mesh->surfaceKernel =
  mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
   "boltzmannSurface2D",
   kernelInfo);

  printf("Compiling Unsplit  pml surface kernel\n");
  mesh->pmlSurfaceKernel =
  mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurface2D.okl",
   "boltzmannPmlSurface2D",
   kernelInfo);
       
//SARK STAGE UPDATE
  printf("compiling SAAB3 non-pml  update kernel\n");
  mesh->updateKernel =
  mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
    "boltzmannSAAB3Update2D",
    kernelInfo); 


    //SARK STAGE UPDATE
  printf("compiling SAAB3 non-pml  update kernel\n");
  mesh->pmlUpdateKernel =
  mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdate2D.okl",
    "boltzmannSAAB3PmlUpdate2D",
    kernelInfo); 

  mesh->haloExtractKernel =
  mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
    "meshHaloExtract2D",
    kernelInfo);   


  }

  
}




