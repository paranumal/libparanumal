#include "boltzmann2D.h"
//#include <complex.h>

// NBN: toggle use of 2nd stream
#define USE_2_STREAMS
// #undef USE_2_STREAMS

// void  ContourIntegralSARK3(mesh2D *mesh){
// dfloat complex z1 = 1.0 + 3.0 * I;
// dfloat complex z2 = 1.0 - 4.0 * I;

// dfloat a21 = mesh->rk3a[1][0]; 
// dfloat a31 = mesh->rk3a[2][0]; 
// dfloat a32 = mesh->rk3a[2][1]; 
// //
//  dfloat b1 = mesh->rk3b[0];     
//  dfloat b2 = mesh->rk3b[1];      
//  dfloat b3 = mesh->rk3b[2]; 
//  //
//  dfloat c1 = mesh->rk3c[0]; 
//  dfloat c2 = mesh->rk3c[1]; 
//  dfloat c3 = mesh->rk3c[2]; 


// iint M = 64; 

// dfloat coef = -mesh->tauInv;
// dfloat h    =  mesh->dt; 

// // Create  number of M points on cylinder on complex plane
// dfloat complex R[M]; 
// for(int s=0; s<M; s++)
// {

//  R[s] = exp(1.0*I*M_PI*( ((s+1.)-0.5)/(M+1)));

// }


     // //  Exponential Coefficients
     // mesh->sarka[1][0] = -(a21*exp(c2*coef*h)*(exp(-c2*coef*h) - 1.))/(c2*coef*h); // a21
     // mesh->sarka[2][0] = -(a31*exp(c3*coef*h)*(exp(-c3*coef*h) - 1.))/(c3*coef*h); // a31
     // mesh->sarka[2][1] = -(a32*exp(c3*coef*h)*(exp(-c3*coef*h) - 1.))/(c3*coef*h); // a32 

     // // If 1/tau*h is too small say <<1, need to write in terms of Taylor coefficients
     // mesh->sarkb[0] =   (exp(coef*h)*((exp(-coef*h)*(c2 + c3 - c2*c3 + c2*c3*exp(coef*h) - 1.))
     //                    /(coef*(c1 - c2)*(c1 - c3)) + (exp(-coef*h)*(2.*exp(coef*h) - 2.) - coef*h*exp(-coef*h)*(c2*exp(coef*h) - c3 - c2 + c3*exp(coef*h) + 2.))
     //                    /(pow(coef,3.)*pow(h,2)*(c1 - c2)*(c1 - c3))))/h;
     // mesh->sarkb[1] =  -(exp(coef*h)*((exp(-coef*h)*(c1 + c3 - c1*c3 + c1*c3*exp(coef*h) - 1.))
     //                    /(coef*(c1 - c2)*(c2 - c3)) + (exp(-coef*h)*(2.*exp(coef*h) - 2.) - coef*h*exp(-coef*h)*(c1*exp(coef*h) - c3 - c1 + c3*exp(coef*h) + 2.))
     //                    /(pow(coef,3)*pow(h,2)*(c1 - c2)*(c2 - c3))))/h;
     // mesh->sarkb[2] =   (exp(coef*h)*((exp(-coef*h)*(c1 + c2 - c1*c2 + c1*c2*exp(coef*h) - 1.))/(coef*(c1 - c3)*(c2 - c3)) + (exp(-coef*h)*(2.*exp(coef*h) - 2.) - coef*h*exp(-coef*h)*(c1*exp(coef*h) - c2 - c1 + c2*exp(coef*h) + 2.))
     //                    /(pow(coef,3)*pow(h,2)*(c1 - c3)*(c2 - c3))))/h;
     // //
     // mesh->sarke[0] = exp(coef*h*c2); 
     // mesh->sarke[1] = exp(coef*h*c3); 
     // mesh->sarke[2] = exp(coef*h*1.0);
     
     // // PML Region Coefficients
     // coef = 0.5*coef; 
     //  //  Exponential Coefficients
     // mesh->sarkpmla[1][0] = -(a21*exp(c2*coef*h)*(exp(-c2*coef*h) - 1.))/(c2*coef*h); // a21
     // mesh->sarkpmla[2][0] = -(a31*exp(c3*coef*h)*(exp(-c3*coef*h) - 1.))/(c3*coef*h); // a31
     // mesh->sarkpmla[2][1] = -(a32*exp(c3*coef*h)*(exp(-c3*coef*h) - 1.))/(c3*coef*h); // a32 

     // // If 1/tau*h is too small say <<1, need to write in terms of Taylor coefficients
     // mesh->sarkpmlb[0] =   (exp(coef*h)*((exp(-coef*h)*(c2 + c3 - c2*c3 + c2*c3*exp(coef*h) - 1.))
     //                    /(coef*(c1 - c2)*(c1 - c3)) + (exp(-coef*h)*(2.*exp(coef*h) - 2.) - coef*h*exp(-coef*h)*(c2*exp(coef*h) - c3 - c2 + c3*exp(coef*h) + 2.))
     //                    /(pow(coef,3.)*pow(h,2)*(c1 - c2)*(c1 - c3))))/h;
     // mesh->sarkpmlb[1] =  -(exp(coef*h)*((exp(-coef*h)*(c1 + c3 - c1*c3 + c1*c3*exp(coef*h) - 1.))
     //                    /(coef*(c1 - c2)*(c2 - c3)) + (exp(-coef*h)*(2.*exp(coef*h) - 2.) - coef*h*exp(-coef*h)*(c1*exp(coef*h) - c3 - c1 + c3*exp(coef*h) + 2.))
     //                    /(pow(coef,3)*pow(h,2)*(c1 - c2)*(c2 - c3))))/h;
     // mesh->sarkpmlb[2] =   (exp(coef*h)*((exp(-coef*h)*(c1 + c2 - c1*c2 + c1*c2*exp(coef*h) - 1.))/(coef*(c1 - c3)*(c2 - c3)) + (exp(-coef*h)*(2.*exp(coef*h) - 2.) - coef*h*exp(-coef*h)*(c1*exp(coef*h) - c2 - c1 + c2*exp(coef*h) + 2.))
     //                    /(pow(coef,3)*pow(h,2)*(c1 - c3)*(c2 - c3))))/h;
     // //
     // mesh->sarkpmle[0] = exp(coef*h*c2); 
     // mesh->sarkpmle[1] = exp(coef*h*c3); 
     // mesh->sarkpmle[2] = exp(coef*h*1.0);


// }

void boltzmannSplitPmlSetup2D(mesh2D *mesh){

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

  mesh->pmlNT    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				    sizeof(dfloat));
  mesh->rhspmlNT = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				    sizeof(dfloat));
  mesh->respmlNT = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				    sizeof(dfloat));
  
  mesh->sigmax= (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
  mesh->sigmay= (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
  
  // set temperature, gas constant, wave speeds

   

#if PML_ENABLED
    printf("Starting initial conditions\n");
    dfloat Ma = 0.1;
    mesh->RT  = 9.0;
    mesh->sqrtRT = sqrt(mesh->RT);  

    dfloat Re = 1000/mesh->sqrtRT; 
    mesh->tauInv = mesh->sqrtRT * Re / Ma;
    dfloat nu = mesh->RT/mesh->tauInv; 
    // Uniform Flow Condition
    dfloat rho = 1, u = 1., v = 0.; 
    dfloat sigma11 = 0, sigma12 = 0, sigma22 = 0;
    //
    //
    mesh->finalTime = 20.;

 #else
    dfloat Ma = 0.1;
    mesh->RT  = 9.0;
    mesh->sqrtRT = sqrt(mesh->RT);  

    dfloat Re = 100/mesh->sqrtRT; 
    mesh->tauInv = mesh->sqrtRT * Re / Ma;
    dfloat nu = mesh->RT/mesh->tauInv; 
    // Create Periodic Boundaries
    printf("Creating periodic connections if exist \n");
    dfloat xper = 1.0;   dfloat yper = 0.0;
    boltzmannPeriodic2D(mesh,xper,yper);

    printf("starting initial conditions\n"); //Zero Flow Conditions
    dfloat rho = 1, u = 0., v = 0.; 
    dfloat sigma11 = 0, sigma12 = 0, sigma22 = 0;
     //
    mesh->finalTime = 10.;

 #endif 

  
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

#if 0
      boltzmannCavitySolution2D(x, y, t,
				mesh->q+cnt, mesh->q+cnt+1, mesh->q+cnt+2);
#endif

#if 0
      boltzmannGaussianPulse2D(x, y, t,
			       mesh->q+cnt,
			       mesh->q+cnt+1,
			       mesh->q+cnt+2,
			       mesh->q+cnt+3,
			       mesh->q+cnt+4,
			       mesh->q+cnt+5);
#endif
      mesh->q[cnt+0] = q1bar; // uniform density, zero flow
      mesh->q[cnt+1] = ramp*q2bar;
      mesh->q[cnt+2] = ramp*q3bar;
      mesh->q[cnt+3] = ramp*ramp*q4bar;
      mesh->q[cnt+4] = ramp*ramp*q5bar;
      mesh->q[cnt+5] = ramp*ramp*q6bar;
    
      cnt += mesh->Nfields;

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

//
  

  printf("starting parameters\n");
  
  // set penalty parameter
  //  mesh->Lambda2 = 0.5/(sqrt(3.)*mesh->sqrtRT);
  mesh->Lambda2 = 0.5/(mesh->sqrtRT);

  // find elements with center inside PML zone
  dfloat xmin = -4, xmax = 8, ymin = -4, ymax = 4;
  dfloat xsigma = 80, ysigma = 80;
  //    dfloat xsigma = 0, ysigma = 0;
  
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
    
    // add element outside [xmin,xmax]x[ymin,ymax] to pml
#if 0
    if(cx<xmin || cx>xmax)
      mesh->sigmax[e] = xsigma;
    if(cy<ymin || cy>ymax)
      mesh->sigmay[e] = ysigma;
#endif

    iint isPml = 0;
    
    for(iint n=0;n<mesh->Np;++n){
      dfloat x = mesh->x[n + e*mesh->Np];
      dfloat y = mesh->y[n + e*mesh->Np];
      //      if(cx<xmax+1 && cx>xmin-1 && cy<ymax+1 && cy>ymin-1){

      if(cx>xmax){
	//  mesh->sigmax[mesh->Np*e + n] = xsigma;
	mesh->sigmax[mesh->Np*e + n] = xsigma*pow(x-xmax,2);
	isPml = 1;
      }
      if(cx<xmin){
	//  mesh->sigmax[mesh->Np*e + n] = xsigma;
	mesh->sigmax[mesh->Np*e + n] = xsigma*pow(x-xmin,2);
	isPml = 1;
      }
      if(cy>ymax){
	//	  mesh->sigmay[mesh->Np*e + n] = ysigma;
	  mesh->sigmay[mesh->Np*e + n] = ysigma*pow(y-ymax,2);
	  isPml = 1;
      }
      if(cy<ymin){
	//  mesh->sigmay[mesh->Np*e + n] = ysigma;
	mesh->sigmay[mesh->Np*e + n] = ysigma*pow(y-ymin,2);
	isPml = 1;
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

      // A = 0.5*h*L
      // => J*2 = 0.5*h*sJ*2
      // => h = 2*J/sJ
      
      dfloat hest = 2./(sJ*invJ);

      hmin = mymin(hmin, hest);
      hmax = mymax(hmax, hest);
    }
  }


  
dfloat cfl = 0.5; 

//
  dfloat magVelocity = sqrt(q2bar*q2bar+q3bar*q3bar)/(q1bar/mesh->sqrtRT);
  magVelocity = mymax(magVelocity,1.0); // Correction for initial zero velocity
  //
  dfloat dtex = hmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*mesh->sqrtRT);
  dfloat dtim = 4./(mesh->tauInv*magVelocity);

  // AK: Set time step size
#if TIME_DISC==LSERK
      printf("Time discretization method: LSERK with CFL: %.2f \n",cfl);
      dfloat dt = mesh->dtfactor*cfl*mymin(dtex,dtim);

      printf("dt = %.4e explicit-dt = %.4e , implicit-dt= %.4e  ratio= %.4e\n", dt,dtex,dtim, dtex/dtim);
 
#elif TIME_DISC==SARK3 
      printf("Time discretization method: SARK with CFL: %.2f \n",cfl);
      dfloat dt = mesh->dtfactor*cfl*mymin(dtex,dtim);
      printf("dt = %.4e explicit-dt = %.4e , implicit-dt= %.4e  ratio= %.4e\n", dt,dtex,dtim, dtex/dtim);

#elif TIME_DISC==LSIMEX
      printf("Time discretization method: Low Storage IMEX  with CFL: %.2f \n",cfl);
      dfloat dt = mesh->dtfactor*cfl*mymin(dtex,dtim);
      
      printf("dt = %.4e explicit-dt = %.4e , implicit-dt= %.4e  ratio= %.4e\n", dt,dtex,dtim, dtex/dtim);
#endif



  printf("hmin = %g\n", hmin);
  printf("hmax = %g\n", hmax);
  printf("cfl = %g\n", cfl);
  printf("dt = %g ", dt);
  printf("max wave speed = %g\n", sqrt(3.)*mesh->sqrtRT);
  
  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
 

   //mesh->dt = 1e-4;
  
  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  mesh->dt = mesh->finalTime/mesh->NtimeSteps;

  // errorStep
  mesh->errorStep = 1000;

  printf("Nsteps = %d with dt = %.8e\n", mesh->NtimeSteps, mesh->dt);

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // use rank to choose DEVICE
  // sprintf(deviceConfig, "mode = CUDA, deviceID = %d", (rank+1)%3);
 sprintf(deviceConfig, "mode = OpenCL, deviceID = 1, platformID = 0");
 // sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 0");

 // sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");	 


  occa::kernelInfo kernelInfo;

  meshOccaSetup2D(mesh, deviceConfig,  kernelInfo);




  #if TIME_DISC==LSERK 
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

  #elif TIME_DISC==SARK3
    mesh->o_qold =
        mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);

    // mesh->o_rhsq =
    //     mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);    
    mesh->o_rhsq2 =
        mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);
    mesh->o_rhsq3 =
        mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);
    
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

        //
      for(int i=0; i<5; i++){
        for(int j=0; j<5; j++){
          mesh->sarka[i][j] = 0.0;
          mesh->sarkpmla[i][j] = 0.0 ; 
        }
        mesh->sarkb[i] = 0.0;
        mesh->sarke[i] = 0.0;
        mesh->sarkpmlb[i] = 0.0;
        mesh->sarkpmle[i] = 0.0;     
      }

     dfloat coef = -mesh->tauInv;
     dfloat  h   = mesh->dt; 
     // Base Runge-Kutta Method
     dfloat a21 = 1./3.;   dfloat a31 = -3./16. ;    dfloat a32 = 15./16.;
     dfloat b1 = 1./6.;    dfloat b2 = 3./10.;       dfloat b3 = 8./15.; 
     dfloat c1 = 0.;       dfloat c2 = 1./3.;        dfloat c3 = 3./4.; 
      // Base Method
     mesh->rk3a[0][0] = 0.;  mesh->rk3a[1][0] = a21;  mesh->rk3a[2][0] = a31;   mesh->rk3a[2][1] = a32; 
     mesh->rk3b[0] = b1;     mesh->rk3b[1] = b2;      mesh->rk3b[2] = b3; 
     mesh->rk3c[0] = c1;     mesh->rk3c[1] = c2;      mesh->rk3c[2] = c3; 

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
     mesh->sarke[0] = exp(coef*h*c2); 
     mesh->sarke[1] = exp(coef*h*c3); 
     mesh->sarke[2] = exp(coef*h*1.0);
     
     // PML Region Coefficients
     coef = 0.5*coef; 
      //  Exponential Coefficients
     mesh->sarkpmla[1][0] = -(a21*exp(c2*coef*h)*(exp(-c2*coef*h) - 1.))/(c2*coef*h); // a21
     mesh->sarkpmla[2][0] = -(a31*exp(c3*coef*h)*(exp(-c3*coef*h) - 1.))/(c3*coef*h); // a31
     mesh->sarkpmla[2][1] = -(a32*exp(c3*coef*h)*(exp(-c3*coef*h) - 1.))/(c3*coef*h); // a32 

     // If 1/tau*h is too small say <<1, need to write in terms of Taylor coefficients
     mesh->sarkpmlb[0] =   (exp(coef*h)*((exp(-coef*h)*(c2 + c3 - c2*c3 + c2*c3*exp(coef*h) - 1.))
                        /(coef*(c1 - c2)*(c1 - c3)) + (exp(-coef*h)*(2.*exp(coef*h) - 2.) - coef*h*exp(-coef*h)*(c2*exp(coef*h) - c3 - c2 + c3*exp(coef*h) + 2.))
                        /(pow(coef,3.)*pow(h,2)*(c1 - c2)*(c1 - c3))))/h;
     mesh->sarkpmlb[1] =  -(exp(coef*h)*((exp(-coef*h)*(c1 + c3 - c1*c3 + c1*c3*exp(coef*h) - 1.))
                        /(coef*(c1 - c2)*(c2 - c3)) + (exp(-coef*h)*(2.*exp(coef*h) - 2.) - coef*h*exp(-coef*h)*(c1*exp(coef*h) - c3 - c1 + c3*exp(coef*h) + 2.))
                        /(pow(coef,3)*pow(h,2)*(c1 - c2)*(c2 - c3))))/h;
     mesh->sarkpmlb[2] =   (exp(coef*h)*((exp(-coef*h)*(c1 + c2 - c1*c2 + c1*c2*exp(coef*h) - 1.))/(coef*(c1 - c3)*(c2 - c3)) + (exp(-coef*h)*(2.*exp(coef*h) - 2.) - coef*h*exp(-coef*h)*(c1*exp(coef*h) - c2 - c1 + c2*exp(coef*h) + 2.))
                        /(pow(coef,3)*pow(h,2)*(c1 - c3)*(c2 - c3))))/h;
     //
     mesh->sarkpmle[0] = exp(coef*h*c2); 
     mesh->sarkpmle[1] = exp(coef*h*c3); 
     mesh->sarkpmle[2] = exp(coef*h*1.0);


    
  #elif TIME_DISC==LSIMEX
    mesh->o_qY =    
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
     // pml variables
   mesh->o_qZ =    
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
     // 
   mesh->o_qS =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);
    

      #if PML_ENABLED
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
      #endif

  
  #elif TIME_DISC==SAAB
   mesh->o_rhsq2 =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->q);
   
   mesh->o_rhsq3 =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->q);
  
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
      dfloat dtc = mesh->dt; 
      //
      mesh->saab[0] = (exp(cc*dtc) - (5*cc*dtc)/2 - 3*pow(cc,2)*pow(dtc,2 )+ pow(cc,2)*pow(dtc,2)*exp(cc*dtc) + (3*cc*dtc*exp(cc*dtc))/2 - 1)/(pow(cc,3)*pow(dtc,2));
      mesh->saab[1] = (4*cc*dtc - 2*exp(cc*dtc) + 3*pow(cc,2)*pow(dtc,2 )- 2*cc*dtc*exp(cc*dtc) + 2)/(pow(cc,3)*pow(dtc,2));
      mesh->saab[2] = -((3*cc*dtc)/2 - exp(cc*dtc) + pow(cc,2)*pow(dtc,2 )- (cc*dtc*exp(cc*dtc))/2 + 1)/(pow(cc,3)*pow(dtc,2));

      // mesh->saab[0] = (pow(cc,3)*pow(dtc,4))/18. + (19.*pow(cc,2)*pow(dtc,3))/80. + (19.*cc*pow(dtc,2))/24. + (23.*dtc)/12.;
      // mesh->saab[1] = -( (7.*pow(cc,3)*pow(dtc,4))/360. + (pow(cc,2)*pow(dtc,3))/10. + (5.*cc*pow(dtc,2))/12.  + (4.*dtc)/3. );
      // mesh->saab[2] = (pow(cc,3)*pow(dtc,4))/180. + (7.*pow(cc,2)*pow(dtc,3))/240. + (cc*pow(dtc,2))/8. + (5.*dtc)/12.;



      //Define exp(tauInv*dt) 
      mesh->saabexp = exp(-mesh->tauInv*dt);

  #endif  
  
  



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

  int NblockV = 256/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 256/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

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

  #if TIME_DISC==LSIMEX // Required Info Local Reduction
   #if CUBATURE_ENABLED==0
  int G = pow(2, ceil( log(mesh->Np * NblockV / 32) /log(2) )) ;
  int S = 32; 
  #else
  int G = pow(2, ceil(log(mesh->cubNp * NblockV / 32) /log(2) )) ;
  int S = 32; 
  #endif

  printf("\n G = %d   S = %d NblockV = %d  Ncub = %d  Np = %d \n",G,S, NblockV, mesh->cubNp,mesh->Np);

  kernelInfo.addDefine("p_G", G);
  kernelInfo.addDefine("p_S", S);

  int imex_iter_max = 25;  
  dfloat imex_tol   = 1.0e-7;  
  dfloat nodetol     = 1e-12; 
  kernelInfo.addDefine("p_LSIMEX_MAXITER", (int) imex_iter_max);
  kernelInfo.addDefine("p_LSIMEX_TOL", (float) imex_tol);
  kernelInfo.addDefine("p_NODETOL", (float) nodetol);
  #endif
 

printf("\n NblockV = %d  Ncub = %d  Np = %d \n", NblockV, mesh->cubNp,mesh->Np);




 #if TIME_DISC==LSERK
    #if CUBATURE_ENABLED
      printf("Compiling LSERK volume kernel with cubature integration\n");
      mesh->volumeKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannVolume2D.okl",
                   "boltzmannVolumeCub2D",
                   kernelInfo);
      
      printf("Compiling LSERK pml volume kernel with cubature integration\n");
      mesh->pmlVolumeKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannVolume2D.okl",
               "boltzmannSplitPmlVolumeCub2D",
               kernelInfo);
      
       printf("Compiling LSERK relaxation kernel with cubature integration\n");
       mesh->relaxationKernel =
       mesh->device.buildKernelFromSource("okl/boltzmannRelaxation2D.okl",
               "boltzmannRelaxationCub2D",
               kernelInfo); 

      printf("Compiling LSERK pml relaxation kernel with cubature integration\n");
       mesh->pmlRelaxationKernel =
       mesh->device.buildKernelFromSource("okl/boltzmannRelaxation2D.okl",
               "boltzmannSplitPmlRelaxationCub2D",
               kernelInfo);   
    #else
      printf("Compiling pml volume kernel with nodal collocation for nonlinear term\n");
      mesh->volumeKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannVolume2D.okl",
                 "boltzmannVolume2D",
                 kernelInfo);

      printf("Compiling pml volume kernel with nodal collocation for nonlinear term\n");
      mesh->pmlVolumeKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannVolume2D.okl",
               "boltzmannSplitPmlVolume2D",
               kernelInfo); 
    #endif


    printf("Compiling surface kernel\n");
    mesh->surfaceKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannSurface2D.okl",
                 "boltzmannSurface2D",
                 kernelInfo);

    printf("Compiling pml surface kernel\n");
    mesh->pmlSurfaceKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannSurface2D.okl",
                 "boltzmannSplitPmlSurface2D",
                 kernelInfo);
 
    printf("Compiling update kernel\n");
    mesh->updateKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
                 "boltzmannLSERKUpdate2D",
                 kernelInfo);
    
    printf("Compiling pml update kernel\n");
    mesh->pmlUpdateKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
                 "boltzmannLSERKSplitPmlUpdate2D",
                 kernelInfo);
  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource("okl/meshHaloExtract2D.okl",
               "meshHaloExtract2D",
               kernelInfo);

#elif TIME_DISC==SARK3 

   #if CUBATURE_ENABLED
      printf("Compiling SA volume kernel with cubature integration\n");
      mesh->volumeKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannVolume2D.okl",
                   "boltzmannVolumeCub2D",
                   kernelInfo);

      printf("Compiling SA pml volume kernel with cubature integration\n");
      mesh->pmlVolumeKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannVolume2D.okl",
               "boltzmannSASplitPmlVolumeCub2D",
               kernelInfo);

       printf("Compiling SA relaxation kernel with cubature integration\n");
       mesh->relaxationKernel =
       mesh->device.buildKernelFromSource("okl/boltzmannRelaxation2D.okl",
               "boltzmannSARelaxationCub2D",
               kernelInfo); 

      printf("Compiling SA pml relaxation kernel with cubature integration\n");
       mesh->pmlRelaxationKernel =
       mesh->device.buildKernelFromSource("okl/boltzmannRelaxation2D.okl",
               "boltzmannSASplitPmlRelaxationCub2D",
               kernelInfo); 


    #else
      printf("Compiling SA volume kernel\n");
      mesh->volumeKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannVolume2D.okl",
                 "boltzmannSAVolume2D",
                 kernelInfo);

      printf("Compiling SAAB pml volume kernel\n");
      mesh->pmlVolumeKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannVolume2D.okl",
               "boltzmannSASplitPmlVolume2D",
               kernelInfo); 
    #endif


  printf("Compiling surface kernel\n");
  mesh->surfaceKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannSurface2D.okl",
               "boltzmannSurface2D",
               kernelInfo);

  printf("Compiling pml surface kernel\n");
  mesh->pmlSurfaceKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannSurface2D.okl",
               "boltzmannSplitPmlSurface2D",
               kernelInfo);


    //SARK STAGE UPDATE
    printf("compiling SARK non-pml  update kernel\n");
    mesh->updateKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
                 "boltzmannSARK3Update2D",
                 kernelInfo); 


   printf("compiling SARK non-pml  stage update kernel\n");
    mesh->updateStageKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
                 "boltzmannSARK3StageUpdate2D",
                 kernelInfo); 

     //SARK STAGE UPDATE
    printf("compiling SARK non-pml  update kernel\n");
    mesh->pmlUpdateKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
                 "boltzmannSARK3SplitPmlUpdate2D",
                 kernelInfo); 


   printf("compiling SARK non-pml  stage update kernel\n");
    mesh->pmlUpdateStageKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
                 "boltzmannSARK3SplitPmlStageUpdate2D",
                 kernelInfo); 

 mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource("okl/meshHaloExtract2D.okl",
               "meshHaloExtract2D",
               kernelInfo);   


 #elif TIME_DISC==LSIMEX
 
 // RESIDUAL UPDATE KERNELS
  printf("Compiling LSIMEX non-pml residual update kernel\n");
  mesh->residualUpdateKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
               "boltzmannLSIMEXResidualUpdate2D",
               kernelInfo);

  printf("Compiling pml residual update kernel\n");
  mesh->pmlResidualUpdateKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
               "boltzmannLSIMEXSplitPmlResidualUpdate2D",
               kernelInfo);

  
  printf("Compiling LSIMEX non-pml implicit update kernel\n");
  mesh->implicitUpdateKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
               "boltzmannLSIMEXImplicitUpdate2D",
               kernelInfo);

  printf("Compiling LSIMEX pml implicit update kernel\n");
  mesh->pmlImplicitUpdateKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
               "boltzmannLSIMEXSplitPmlImplicitUpdate2D",
               kernelInfo);
     
    


    #if CUBATURE_ENABLED
       printf("Compiling LSIMEX non-pml Implicit Iteration  kernel\n");

         mesh->implicitSolveKernel = 
         mesh->device.buildKernelFromSource("okl/boltzmannLSIMEXImplicitSolve2D.okl",
                 "boltzmannLSIMEXImplicitSolveCub2D",
                 kernelInfo); 
      
         printf("Compiling LSIMEX pml Implicit Iteration  kernel\n");
         mesh->pmlImplicitSolveKernel = 
         mesh->device.buildKernelFromSource("okl/boltzmannLSIMEXImplicitSolve2D.okl",
                 "boltzmannLSIMEXSplitPmlImplicitSolveCub2D",
                 kernelInfo); 
    #else
        mesh->implicitSolveKernel = 
         mesh->device.buildKernelFromSource("okl/boltzmannLSIMEXImplicitSolve2D.okl",
                 "boltzmannLSIMEXImplicitSolve2D",
                 kernelInfo); 
      
         printf("Compiling LSIMEX pml Implicit Iteration  kernel\n");
         mesh->pmlImplicitSolveKernel = 
         mesh->device.buildKernelFromSource("okl/boltzmannLSIMEXImplicitSolve2D.okl",
                 "boltzmannLSIMEXSplitPmlImplicitSolve2D",
                 kernelInfo); 
    #endif


    printf("Compiling LSIMEX volume kernel no stiff term!\n");
    mesh->volumeKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannVolume2D.okl",
                 "boltzmannVolumeCub2D",
                 kernelInfo);
      
    printf("Compiling LSIMEX pml volume kernel with cubature integration\n");
    mesh->pmlVolumeKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannVolume2D.okl",
             "boltzmannSplitPmlVolumeCub2D",
             kernelInfo);

    printf("Compiling surface kernel\n");
    mesh->surfaceKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannSurface2D.okl",
               "boltzmannSurface2D",
               kernelInfo);

    printf("Compiling pml surface kernel\n");
    mesh->pmlSurfaceKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannSurface2D.okl",
               "boltzmannSplitPmlSurface2D",
               kernelInfo);


   printf("Compiling LSIMEX non-pml update kernel\n");
  mesh->updateKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
               "boltzmannLSIMEXUpdate2D",
               kernelInfo);

  printf("Compiling LSIMEX pml update kernel\n");
  mesh->pmlUpdateKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
               "boltzmannLSIMEXSplitPmlUpdate2D",
               kernelInfo);



mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource("okl/meshHaloExtract2D.okl",
               "meshHaloExtract2D",
               kernelInfo);

#endif

  
}
