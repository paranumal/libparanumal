#include "boltzmann2D.h"

// NBN: toggle use of 2nd stream
#define USE_2_STREAMS
// #undef USE_2_STREAMS

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

  // // // Define problem parameters
  // dfloat Re = 1000.;
  // dfloat Ma = 0.1;

  // dfloat Lc = 1.0;
  // dfloat Uc = 1.0; 
  // mesh->RT  = Uc*Uc/(Ma*Ma);; 
  // mesh->sqrtRT = sqrt(mesh->RT);   


  // dfloat nu = Uc*Lc/Re; 
  // mesh->tauInv = mesh->RT/nu;

  // printf("Tauinv = %.4f\n", mesh->tauInv);

  // printf("starting initial conditions\n");
  // dfloat rho = 1, u = 1.0, v = 0; //u = 1.f/sqrt(2.f), v = 1.f/sqrt(2.f); 
  // dfloat sigma11 = 0, sigma12 = 0, sigma22 = 0;
  // //  dfloat ramp = 0.5*(1.f+tanh(10.f*(0-.5f)));
  // dfloat ramp, drampdt;
  // boltzmannRampFunction2D(0, &ramp, &drampdt);
  
  
     

  

  




  mesh->RT = 9.0;
  mesh->sqrtRT = sqrt(mesh->RT);  
  
  printf("starting initial conditions\n");
   dfloat rho = 1, u = 1, v = 0; //u = 1.f/sqrt(2.f), v = 1.f/sqrt(2.f); 
  dfloat sigma11 = 0, sigma12 = 0, sigma22 = 0;
  //  dfloat ramp = 0.5*(1.f+tanh(10.f*(0-.5f)));
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

  printf("starting parameters\n");

  // set BGK collision relaxation rate
  // nu = R*T*tau
  // 1/tau = RT/nu
  
  dfloat Re = 300; 
  dfloat Ma = 0.1;
  mesh->tauInv = mesh->sqrtRT * Re / Ma;

  dfloat nu = mesh->RT/mesh->tauInv; 

  





  
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

  dfloat cfl = .5; // depends on the stability region size (was .4)

  // dt ~ cfl (h/(N+1)^2)/(Lambda^2*fastest wave speed)
  // too small ???
  

  // AK: Set time step size
  #if TIME_DISC==LSERK
  printf("Time discretization method: LSERK with CFL: %.2f \n",cfl);

  dfloat magVelocity = sqrt(q2bar*q2bar+q3bar*q3bar)/(q1bar/mesh->sqrtRT);

  dfloat dtex = cfl*hmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*mesh->sqrtRT);

  dfloat dtim = cfl*4./(mesh->tauInv*magVelocity);

  dfloat dt = mymin(dtex,dtim);

  printf("dt = %.4e explicit-dt = %.4e , implicit-dt= %.4e \n", dt,dtex,dtim);

  #elif TIME_DISC==SAAB
  printf("Time discretization method: SAAB  with 1/3*CFL(%.2f)\n",cfl);
   // Adams-Bashforth stability region is approximated as 1/3 of LSERK scheme. 
  dfloat magVelocity = sqrt(q2bar*q2bar+q3bar*q3bar)/(q1bar/mesh->sqrtRT);

  dfloat dtex = 1./3. * cfl*hmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*mesh->sqrtRT);

  dfloat dtim = 1./3. * cfl*4./(mesh->tauInv*magVelocity);

  dfloat dt = mymin(dtex,dtim);

  printf("dt = %.4e explicit-dt = %.4e , implicit-dt= %.4e \n", dt,dtex,dtim);
  
  #elif TIME_DISC==LSIMEX
  printf("Time discretization method: Low Storage IMEX  with CFL: %.2f \n",cfl);
  dfloat magVelocity = sqrt(q2bar*q2bar+q3bar*q3bar)/(q1bar/mesh->sqrtRT);

  // dfloat dt = cfl*mymin(hmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*mesh->sqrtRT),
  //             4./(mesh->tauInv*magVelocity));

  dfloat dt = cfl*hmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*mesh->sqrtRT);
  #endif



  printf("hmin = %g\n", hmin);
  printf("hmax = %g\n", hmax);
  printf("cfl = %g\n", cfl);
  printf("dt = %g ", dt);
  printf("max wave speed = %g\n", sqrt(3.)*mesh->sqrtRT);
  
  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  //
  mesh->finalTime = 10.0;
  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  mesh->dt = mesh->finalTime/mesh->NtimeSteps;

  // errorStep
  mesh->errorStep = 1000;

  printf("dt = %g\n", mesh->dt);

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // use rank to choose DEVICE
  // sprintf(deviceConfig, "mode = CUDA, deviceID = %d", (rank+1)%3);
  sprintf(deviceConfig, "mode = OpenCL, deviceID = 1, platformID = 0");
  //    sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  // sprintf(deviceConfig, "mode = Serial");	 






  occa::kernelInfo kernelInfo;

  meshOccaSetup2D(mesh, deviceConfig,  kernelInfo);




  #if TIME_DISC==LSERK 
  // pml variables
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

  mesh->o_pmlNT =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlNT);
  mesh->o_rhspmlNT =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlNT);
  mesh->o_respmlNT =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->respmlNT);
  
  #elif TIME_DISC==SAAB
// Extra Storage for History, 
   mesh->o_rhsq2 =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);
   mesh->o_rhsq3 =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);

  // pml variables
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
 

  mesh->o_pmlNT =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlNT);
  mesh->o_rhspmlNT =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlNT);
  mesh->o_rhspmlNT2 =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlNT);
  mesh->o_rhspmlNT3 =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlNT);

  #elif TIME_DISC==LSIMEX
   mesh->o_qY =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
// pml variables
  mesh->o_qZ =    
     mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
// 
  mesh->o_qS =
        mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);

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



  mesh->o_pmlNT =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlNT);
  mesh->o_qSnt =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqx);
  mesh->o_qYnt =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlNT);
  mesh->o_qZnt =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->respmlNT);

  #endif  
  
  



  mesh->o_sigmax =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), mesh->sigmax);

  mesh->o_sigmay =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), mesh->sigmay);




 // All Functions Related to Adams-Bashforth
#if TIME_DISC==SAAB
    
         
  // Clasical AB coefficients
  dfloat ab1 = 23./14.*dt;
  dfloat ab2 = -4./3. *dt;
  dfloat ab3 = 5./12. *dt;
  kernelInfo.addDefine("p_ab1",  ab1 );
  kernelInfo.addDefine("p_ab2",  ab2  );
  kernelInfo.addDefine("p_ab3",  ab3 );
  //
  // Modefied AB Coefficients
  dfloat saab1 =  -(2.*exp(-mesh->dt*mesh->tauInv)  + 5.*mesh->dt*mesh->tauInv 
                  - 6.*pow(mesh->dt*mesh->tauInv,2) + 2.*pow(mesh->dt*mesh->tauInv,2)*exp(-mesh->dt*mesh->tauInv) 
                  - 3.*mesh->dt*mesh->tauInv*exp(-mesh->dt*mesh->tauInv) - 2.)
                   /(2.*pow(mesh->dt,2)*pow(mesh->tauInv,3));


  dfloat saab2 = -(3.*pow(mesh->dt*mesh->tauInv,2) - 4.*mesh->dt*mesh->tauInv - 2.*exp(-mesh->dt*mesh->tauInv) 
                 + 2.*mesh->dt*mesh->tauInv*exp(-mesh->dt*mesh->tauInv) + 2.)
                  /(pow(mesh->dt,2)*pow(mesh->tauInv,3));

  dfloat saab3 = (2.*pow(mesh->dt*mesh->tauInv,2) - 3.*mesh->dt*mesh->tauInv 
                - 2.*exp(-mesh->dt*mesh->tauInv) + mesh->dt*mesh->tauInv*exp(-mesh->dt*mesh->tauInv) + 2.)
                 /(2*pow(mesh->dt,2)*pow(mesh->tauInv,3));

   printf("%.5e  %.5e  %.5e\n", saab1,saab2,saab3 ); 


  kernelInfo.addDefine("p_saab1",  saab1 );
  kernelInfo.addDefine("p_saab2",  saab2  );
  kernelInfo.addDefine("p_saab3",  saab3 );

  dfloat expdt = exp(-mesh->tauInv*dt);
  kernelInfo.addDefine("p_expdt",  expdt);


#endif
    




  mesh->nonPmlNelements = nonPmlNelements;
  mesh->pmlNelements = pmlNelements;

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

  #if TIME_DISC==LSIMEX
  int G = pow(2, ceil( log(mesh->Np * NblockV / 32) /log(2) )) ;
  int S = 32; 

  kernelInfo.addDefine("p_G", G);
  kernelInfo.addDefine("p_S", S);

  int imex_iter_max = 50 ;  
  dfloat imex_tol   = 1e-7;  
  dfloat nodetol     = 1e-12; 
  kernelInfo.addDefine("p_LSIMEX_MAXITER", imex_iter_max);
  kernelInfo.addDefine("p_LSIMEX_TOL", imex_tol);
  kernelInfo.addDefine("p_NODETOL", nodetol);
  #endif
 






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

#elif TIME_DISC==SAAB

   #if CUBATURE_ENABLED
      printf("Compiling SAAB volume kernel with cubature integration\n");
      mesh->volumeKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannVolume2D.okl",
                   "boltzmannVolumeCub2D",
                   kernelInfo);

      printf("Compiling SAAB pml volume kernel with cubature integration\n");
      mesh->pmlVolumeKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannVolume2D.okl",
               "boltzmannSplitPmlVolumeCub2D",
               kernelInfo);

       printf("Compiling SAAB relaxation kernel with cubature integration\n");
       mesh->relaxationKernel =
       mesh->device.buildKernelFromSource("okl/boltzmannRelaxation2D.okl",
               "boltzmannSAABRelaxationCub2D",
               kernelInfo); 

      printf("Compiling SAAB pml relaxation kernel with cubature integration\n");
       mesh->pmlRelaxationKernel =
       mesh->device.buildKernelFromSource("okl/boltzmannRelaxation2D.okl",
               "boltzmannSAABSplitPmlRelaxationCub2D",
               kernelInfo); 


    #else
      printf("Compiling SAAB volume kernel\n");
      mesh->volumeKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannVolume2D.okl",
                 "boltzmannSAABVolume2D",
                 kernelInfo);

      printf("Compiling SAAB pml volume kernel\n");
      mesh->pmlVolumeKernel =
      mesh->device.buildKernelFromSource("okl/boltzmannVolume2D.okl",
               "boltzmannSAABSplitPmlVolume2D",
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



   //SAAB FIRST ORDER UPDATE
  printf("compiling non-pml 1st order update kernel\n");
  mesh->updateFirstOrderKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
               "boltzmannSAABUpdateFirst2D",
               kernelInfo);

    //SAAB SECOND ORDER UPDATE
  printf("compiling non-pml 2nd order update kernel\n");
  mesh->updateSecondOrderKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
               "boltzmannSAABUpdateSecond2D",
               kernelInfo);


     //SAAB STAGE UPDATE
  printf("compiling non-pml  update kernel\n");
  mesh->updateKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
               "boltzmannSAABUpdate2D",
               kernelInfo); 



      //SAAB FIRST ORDER UPDATE
  printf("Compiling SAAB pml 1st order update kernel\n");
  mesh->pmlUpdateFirstOrderKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
               "boltzmannSAABSplitPmlUpdateFirst2D",
               kernelInfo);

    //SAAB SECOND ORDER UPDATE
  printf("Compiling SAAB pml 2nd order update kernel\n");
  mesh->pmlUpdateSecondOrderKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
               "boltzmannSAABSplitPmlUpdateSecond2D",
               kernelInfo);


     //SAAB STAGE UPDATE
  printf("Compiling SAAB non-pml update kernel\n");
  mesh->pmlUpdateKernel =
    mesh->device.buildKernelFromSource("okl/boltzmannUpdate2D.okl",
               "boltzmannSAABSplitPmlUpdate2D",
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
         mesh->NRIterationKernel = 
         mesh->device.buildKernelFromSource("okl/boltzmannLSIMEXImplicitIteration2D.okl",
                 "boltzmannLSIMEXImplicitIterationCub2D",
                 kernelInfo); 

         printf("Compiling LSIMEX pml Implicit Iteration  kernel\n");
         mesh->pmlNRIterationKernel = 
         mesh->device.buildKernelFromSource("okl/boltzmannLSIMEXImplicitIteration2D.okl",
                 "boltzmannLSIMEXSplitPmlImplicitIterationCub2D",
                 kernelInfo); 


        printf("Compiling LSIMEX non-pml implicit volume kernel\n");
         mesh->implicitVolumeKernel = 
         mesh->device.buildKernelFromSource("okl/boltzmannLSIMEXImplicit2D.okl",
                 "boltzmannLSIMEXImplicitVolumeCub2D",
                 kernelInfo); 

        printf("Compiling LSIMEX pml implicit volume kernel\n");
         mesh->pmlImplicitVolumeKernel = 
         mesh->device.buildKernelFromSource("okl/boltzmannLSIMEXImplicit2D.okl",
                 "boltzmannLSIMEXSplitPmlImplicitVolumeCub2D",
                 kernelInfo);

    #else
        printf("Compiling LSIMEX non-pml Implicit Iteration  kernel\n");
         mesh->NRIterationKernel = 
         mesh->device.buildKernelFromSource("okl/boltzmannLSIMEXImplicitIteration2D.okl",
                 "boltzmannLSIMEXImplicitIteration2D",
                 kernelInfo); 

         printf("Compiling LSIMEX pml Implicit Iteration  kernel\n");
         mesh->pmlNRIterationKernel = 
         mesh->device.buildKernelFromSource("okl/boltzmannLSIMEXImplicitIteration2D.okl",
                 "boltzmannLSIMEXSplitPmlImplicitIteration2D",
                 kernelInfo); 


        printf("Compiling LSIMEX non-pml implicit volume kernel\n");
         mesh->implicitVolumeKernel = 
         mesh->device.buildKernelFromSource("okl/boltzmannLSIMEXImplicit2D.okl",
                 "boltzmannLSIMEXImplicitVolume2D",
                 kernelInfo); 

        printf("Compiling LSIMEX pml implicit volume kernel\n");
         mesh->pmlImplicitVolumeKernel = 
         mesh->device.buildKernelFromSource("okl/boltzmannLSIMEXImplicit2D.okl",
                 "boltzmannLSIMEXSplitPmlImplicitVolume2D",
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
