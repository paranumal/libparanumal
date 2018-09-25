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

#include "bns.h"

bns_t *bnsSetup(mesh_t *mesh, setupAide &options){
  
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

  // check = options.getArgs("FIXED TIME STEP", bns->fixed_dt);
  // if(!check) printf("WARNING setup file does not include FIXED TIME STEP\n");

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

  bns->readRestartFile = 0; 
  options.getArgs("RESTART FROM FILE", bns->readRestartFile);
  
  bns->writeRestartFile = 0; 
  options.getArgs("WRITE RESTART FILE", bns->writeRestartFile);
  
  if(options.compareArgs("PML INTEGRATION", "COLLOCATION"))
    bns->pmlcubature = 0;
  else
    bns->pmlcubature = 1; 

  
  // Set time discretization scheme:fully explicit or not
  bns->fexplicit = 0; 
  if(options.compareArgs("TIME INTEGRATOR", "LSERK") ) // fully explicit schemes
    bns->fexplicit = 1; 
  // Report / error time steps for fixed dt integration
  options.getArgs("TSTEPS FOR SOLUTION OUTPUT", bns->reportStep);
  options.getArgs("TSTEPS FOR ERROR COMPUTE",   bns->errorStep);

  // Output interval for variable dt integration
  options.getArgs("OUTPUT INTERVAL",   bns->outputInterval);
  
  // steps between force output
  bns->outputForceStep = 0;
  
  options.getArgs("TSTEPS FOR FORCE OUTPUT",   bns->outputForceStep);
  
  if(mesh->rank==0){
    printf("=============WRITING INPUT PARAMETERS===================\n");

    printf("VISCOSITY\t:\t%.2e\n", bns->nu);
    printf("SPEED OF SOUND\t:\t%.2e\n", bns->sqrtRT);
    printf("CFL NUMBER\t:\t%.2e\n", bns->cfl);
    printf("START TIME\t:\t%.2e\n", bns->startTime);
    printf("FINAL TIME\t:\t%.2e\n", bns->finalTime);
    // printf("FIXED DT\t:\t%d\n", bns->fixed_dt);
    printf("PML FORMULATION\t:\t%d\n", bns->pmlFlag);
    if(bns->pmlFlag){
      printf("PML PROFILE N\t:\t%d\n", bns->pmlOrder);
      printf("PML SIGMA X\t:\t%.2e\n", bns->sigmaXmax);
      printf("PML SIGMA Y\t:\t%.2e\n", bns->sigmaYmax);
      if(bns->dim==3)
        printf("PML SIGMA Z\t:\t%.2e\n", bns->sigmaZmax);
      printf("PML CUBATURE\t:\t%d\n", bns->pmlcubature);
    }
    printf("ERROR STEP\t:\t%d\n", bns->errorStep);
    printf("RESTART READ\t:\t%d\n", bns->readRestartFile);
    printf("RESTART WRITE\t:\t%d\n", bns->writeRestartFile);
  }

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
  
  // bns->RT       = Uref*Uref/(bns->Ma*bns->Ma);
  bns->RT     = bns->sqrtRT*bns->sqrtRT;  
  bns->tauInv = bns->RT/bns->nu;

  bns->Re = (fabs(u)>0) ? u/bns->nu   : 1.0; // just for postprocessing , assume unit characteristic length
  bns->Ma = (fabs(u)>0) ? u/bns->sqrtRT : 1.0; // just for postprocessing
  // Set penalty parameter for flux setting
  bns->Lambda2 = 0.5/(bns->sqrtRT);


  dfloat sigma11 = 0.f , sigma12 = 0.f, sigma13 = 0.f;
  dfloat sigma22 = 0.f , sigma23 = 0.f;
  dfloat sigma33 = 0.f;

  dfloat q1bar =0., q2bar =0., q3bar =0., q4bar =0., q5bar = 0.; 
  dfloat q6bar =0., q7bar =0., q8bar =0., q9bar =0., q10bar =0.; 
  // Set time step size
  if(bns->dim==2){
    if (mesh->rank==0) printf("MESH DIMENSION\t:\t%d\n", bns->dim);
    q1bar = rho;
    q2bar = rho*u/bns->sqrtRT;
    q3bar = rho*v/bns->sqrtRT;
    q4bar = (rho*u*v - sigma12)/bns->RT;
    q5bar = (rho*u*u - sigma11)/(sqrt(2.)*bns->RT);
    q6bar = (rho*v*v - sigma22)/(sqrt(2.)*bns->RT);    
  } else{
    if (mesh->rank==0) printf("MESH DIMENSION\t:\t%d\n", bns->dim);
    
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
      // not good for all element types (some have varying geofacs)
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


  printf("ghmin =  %lg\n", ghmin);

  
  printf("dtex = %.5e dtim = %.5e \n", bns->cfl*ghmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*bns->sqrtRT), bns->cfl*1.f/(bns->tauInv));

  // Set multiRate element groups/group  
  if(options.compareArgs("TIME INTEGRATOR", "MRSAAB") ){
    int maxLevels =0; options.getArgs("MAX MRAB LEVELS", maxLevels);

    if (mesh->rank==0) printf("MR MAX LEVELS\t:\t%d\n", maxLevels);

    if(bns->dim==2)
      bns->dt = meshMRABSetup2D(mesh,EtoDT,maxLevels, (bns->finalTime-bns->startTime));
    else
      bns->dt = meshMRABSetup3D(mesh,EtoDT,maxLevels, (bns->finalTime-bns->startTime));

    bns->NtimeSteps =  (bns->finalTime-bns->startTime)/(pow(2,mesh->MRABNlevels-1)*bns->dt);

    if (mesh->rank==0) printf("MR  LEVELS\t:\t%d\n", mesh->MRABNlevels);
  }
  else{
    // printf("MESH DIMENSION\t:\t%d\n", bns->dt);  
    //!!!!!!!!!!!!!! Fix time step to compute the error in postprecessing step  
    // MPI_Allreduce to get global minimum dt
    MPI_Allreduce(&dt, &(bns->dt), 1, MPI_DFLOAT, MPI_MIN, mesh->comm);
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

    dlong localElements =  mesh->Nelements;
    MPI_Allreduce(&localElements, &(bns->totalElements), 1, MPI_LONG, MPI_SUM, mesh->comm);

    // compute samples of q at interpolation nodes
    bns->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    //
    bns->NrkStages = 7; // 7 stages order(54) or 6 stages order(43) 
    
    bns->rkq      = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rkrhsq   = (dfloat*) calloc(bns->NrkStages*mesh->Nelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    bns->rkerr    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*bns->Nfields, sizeof(dfloat));
    
    dlong Ntotal  = mesh->Nelements*mesh->Np*bns->Nfields;
    bns->Nblock   = (Ntotal+blockSize-1)/blockSize;
    // printf("blockSize: %d %d \n", blockSize, bns->Nblock);
    bns->errtmp  =  (dfloat*) calloc(bns->Nblock, sizeof(dfloat)); 
  }

 
  dfloat Gy = 0.0; options.getArgs("BODYFORCE-Y",Gy);
  dfloat time = bns->startTime + 0.0;
  
  // Define Initial Mean Velocity
  dfloat fx, fy, fz, intfx, intfy, intfz;
  bnsBodyForce(time, &fx, &fy, &fz, &intfx, &intfy, &intfz);
  
  // INITIALIZE PROBLEM
  if(options.compareArgs("INITIAL CONDITION", "BROWN-MINION")){
        bnsBrownMinionQuad3D(bns);
  }
  else{

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
	  bns->q[id+1*mesh->Np] = q1bar*intfx/bns->sqrtRT;
	  bns->q[id+2*mesh->Np] = q1bar*intfy/bns->sqrtRT;
	  bns->q[id+3*mesh->Np] = q1bar*intfx*intfy/bns->sqrtRT;
	  bns->q[id+4*mesh->Np] = q1bar*intfx*intfx/(sqrt(2.)*bns->sqrtRT);
	  bns->q[id+5*mesh->Np] = q1bar*intfy*intfy/(sqrt(2.)*bns->sqrtRT);
	}
	if(bns->dim==3 && bns->elementType==QUADRILATERALS){
	  bns->q[id+0*mesh->Np] = q1bar*(2+exp(-40*((x-1)*(x-1) + y*y + z*z))); 
	  bns->q[id+1*mesh->Np] = 0;
	  bns->q[id+2*mesh->Np] = 0;
	  bns->q[id+3*mesh->Np] = 0;
	  
	  bns->q[id+4*mesh->Np] = 0;
	  bns->q[id+5*mesh->Np] = 0;
	  bns->q[id+6*mesh->Np] = 0;
	  
	  bns->q[id+7*mesh->Np] = 0;
	  bns->q[id+8*mesh->Np] = 0;
	  bns->q[id+9*mesh->Np] = 0;
	  
	}
	else if(bns->dim==3){
	  bns->q[id+0*mesh->Np] = q1bar; 
	  bns->q[id+1*mesh->Np] = q1bar*intfx/bns->sqrtRT;
	  bns->q[id+2*mesh->Np] = q1bar*intfy/bns->sqrtRT;
	  bns->q[id+3*mesh->Np] = q1bar*intfz/bns->sqrtRT;
	  
	  bns->q[id+4*mesh->Np] = q1bar*intfx*intfy/bns->sqrtRT;
	  bns->q[id+5*mesh->Np] = q1bar*intfx*intfz/bns->sqrtRT;
	  bns->q[id+6*mesh->Np] = q1bar*intfy*intfz/bns->sqrtRT;
	  
	  bns->q[id+7*mesh->Np] = q1bar*intfx*intfx/(sqrt(2.)*bns->sqrtRT);
	  bns->q[id+8*mesh->Np] = q1bar*intfy*intfy/(sqrt(2.)*bns->sqrtRT);
	  bns->q[id+9*mesh->Np] = q1bar*intfz*intfz/(sqrt(2.)*bns->sqrtRT);
	  
	}
	
      }
    }
  }

  
  // Write Problem Info 
  if(mesh->rank==0){
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

  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();
  
  if(bns->dim==3){
    if(bns->elementType != QUADRILATERALS)
      meshOccaSetup3D(mesh, options, kernelInfo);
    else
      meshOccaSetupQuad3D(mesh, options, kernelInfo);
  }
  else
    meshOccaSetup2D(mesh, options, kernelInfo);
  
  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";   

  // Setup MRAB PML
  if(options.compareArgs("TIME INTEGRATOR","MRSAAB")){

    printf("Preparing Pml for multirate rate\n"); 
    bnsMRABPmlSetup(bns, options);

    mesh->o_MRABelementIds = new occa::memory[mesh->MRABNlevels];
    mesh->o_MRABhaloIds    = new occa::memory[mesh->MRABNlevels];
    for (int lev=0;lev<mesh->MRABNlevels;lev++) {

      if (mesh->MRABNelements[lev]){
        mesh->o_MRABelementIds[lev] = mesh->device.malloc(mesh->MRABNelements[lev]*sizeof(dlong),mesh->MRABelementIds[lev]); 
      }
      if (mesh->MRABNhaloElements[lev]){
        mesh->o_MRABhaloIds[lev]    = mesh->device.malloc(mesh->MRABNhaloElements[lev]*sizeof(dlong), mesh->MRABhaloIds[lev]);
      }
    }
  } else{

    if(bns->pmlFlag){
     printf("Preparing Pml for single rate integrator\n");
      bnsPmlSetup(bns, options); 
    }

    if (mesh->nonPmlNelements)
      mesh->o_nonPmlElementIds = mesh->device.malloc(mesh->nonPmlNelements*sizeof(dlong), mesh->nonPmlElementIds);

  }
  

  bns->Nvort      = 3;   // hold wx, wy, wz
  // Set Iso-surfacing stuf here
  if(options.compareArgs("OUTPUT FILE FORMAT","ISO") && bns->dim==3){
    
    // Only one field is exported for iso-surface to reduce the file size
    bns->isoNfields  = 1;   //1 + (bns->dim) + (1 + bns->dim) ; // p, u.v,w, vort_x, vort_y, vort_z, wort_mag 
    bns->isoMaxNtris = 1.E7; 

    bns->procid = gethostid();

    //
    options.getArgs("ISOSURFACE FIELD ID", bns->isoField); 
    options.getArgs("ISOSURFACE COLOR ID", bns->isoColorField); 
    options.getArgs("ISOSURFACE LEVEL NUMBER", bns->isoNlevels);
    options.getArgs("ISOSURFACE CONTOUR MAX", bns->isoMaxVal); 
    options.getArgs("ISOSURFACE CONTOUR MIN", bns->isoMinVal);


    bns->isoMax    = (bns->dim + bns->isoNfields)*3*bns->isoMaxNtris;
    bns->isoNtris  = (int*) calloc(1, sizeof(int));
    bns->isoq      = (dfloat*) calloc(bns->isoMax, sizeof(dfloat)); 

  
    bns->o_isoq      = mesh->device.malloc(bns->isoMax*sizeof(dfloat), bns->isoq);
    bns->o_isoNtris  = mesh->device.malloc(1*sizeof(int), bns->isoNtris);


   

    // meshParallelGatherScatter(mesh, ogs, o_q);
    // Create all contour levels
    dfloat *isoLevels = (dfloat*) calloc(bns->isoNlevels, sizeof(dfloat));
    for(int l=0;l<bns->isoNlevels;++l)
      isoLevels[l] = bns->isoMinVal + (bns->isoMaxVal-bns->isoMinVal)*l/(dfloat)(bns->isoNlevels-1);



    // GROUP LEVELS of ISOCONTOURS

    int levelsInGroup = 0; 
    options.getArgs("ISOSURFACE GROUP NUMBER", levelsInGroup);

    if(levelsInGroup==0) {printf("Number of levels in each group can not be zero!!!\n");  exit(EXIT_FAILURE);} 
    if(levelsInGroup){

      // Number of groups for isosurfaces
      bns->isoGNgroups        = bns->isoNlevels/(levelsInGroup);  
      if(bns->isoNlevels%(levelsInGroup))
        bns->isoGNgroups++; 

      bns->isoGNlevels        = (int *) calloc(bns->isoGNgroups,sizeof(int));
      bns->isoGLvalues        = (dfloat **) calloc(bns->isoGNgroups,sizeof(dfloat*));

      for(int gr =0; gr<bns->isoGNgroups; gr++)
  {
    int nlevels = (gr+1)*levelsInGroup > bns->isoNlevels ? (bns->isoNlevels%levelsInGroup) : levelsInGroup;  
    bns->isoGNlevels[gr] = nlevels;  
    printf("Isosurface Group %d has %d levels\n", gr, bns->isoGNlevels[gr]);
  }

      // Allocate memory for levels in each group
      for (int gr =0;gr<bns->isoGNgroups;gr++)
        bns->isoGLvalues[gr] = (dfloat *) calloc(bns->isoGNlevels[gr],sizeof(dfloat));

      int sk = 0; 
      for (int gr =0;gr<bns->isoGNgroups;gr++){
        printf("Isosurface Group %d Values\n", gr);        
        for (int l=0;l<bns->isoGNlevels[gr];l++){
          bns->isoGLvalues[gr][l] = isoLevels[sk + l];
          printf("%.4f\t", bns->isoGLvalues[gr][l]);
        }
        printf("\n");
  sk += bns->isoGNlevels[gr]; 
      }

      // Create levels for each group
      bns->o_isoGLvalues     = (occa::memory *) malloc(bns->isoGNgroups*sizeof(occa::memory));
      for (int gr =0;gr<bns->isoGNgroups;gr++)
        bns->o_isoGLvalues[gr] = mesh->device.malloc(bns->isoGNlevels[gr]*sizeof(dfloat),bns->isoGLvalues[gr]);
    
    }

    // Interpolation operators form Np to PlotNp (equisapaced nodes of order >N generally)
    dfloat *plotInterp = (dfloat*) calloc(mesh->plotNp*mesh->Np, sizeof(dfloat));
    for(int n=0;n<mesh->plotNp;++n){
      for(int m=0;m<mesh->Np;++m){
        plotInterp[n+m*mesh->plotNp] = mesh->plotInterp[n*mesh->Np+m];
      }
    }
    bns->o_plotInterp = mesh->device.malloc(mesh->plotNp*mesh->Np*sizeof(dfloat), plotInterp);

    // EToV for local triangulation
    int *plotEToV = (int*) calloc(mesh->plotNp*mesh->Np, sizeof(int));
    for(int n=0;n<mesh->plotNelements;++n){
      for(int m=0;m<mesh->plotNverts;++m){
        plotEToV[n+m*mesh->plotNelements] = mesh->plotEToV[n*mesh->plotNverts+m];
      }
    }
    bns->o_plotEToV = mesh->device.malloc(mesh->plotNp*mesh->Np*sizeof(int), plotEToV);
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
  
    bns->o_rkq =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->rkq);

    bns->o_saveq =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->rkq);

    bns->o_rkrhsq =
      mesh->device.malloc(bns->NrkStages*mesh->Np*mesh->Nelements*bns->Nfields*sizeof(dfloat), bns->rkrhsq);
    bns->o_rkerr =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*bns->Nfields*sizeof(dfloat), bns->rkerr);
  
    // dlong Ntotal    = mesh->Nelements*mesh->Np*bns->Nfields;
    // printf("blockSize: %d  %d %d \n", Ntotal, blockSize, bns->Nblock);
    bns->o_errtmp = mesh->device.malloc(bns->Nblock*sizeof(dfloat), bns->errtmp);

    bns->o_rkA = mesh->device.malloc(bns->NrkStages*bns->NrkStages*sizeof(dfloat), bns->rkA);
    bns->o_rkE = mesh->device.malloc(bns->NrkStages*sizeof(dfloat), bns->rkE);

  }

  bns->Vort      = (dfloat*) calloc(bns->Nvort*mesh->Nelements*mesh->Np, sizeof(dfloat));
  bns->VortMag   = (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
  
  bns->o_Vort    = mesh->device.malloc(bns->Nvort*mesh->Nelements*mesh->Np*sizeof(dfloat), bns->Vort);
  bns->o_VortMag = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), bns->VortMag);

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  int maxCubNodes = mymax(maxNodes,mesh->cubNp);

  kernelInfo["defines/" "p_maxNodes"]= maxNodes;
  kernelInfo["defines/" "p_maxCubNodes"]= maxCubNodes;

  if(bns->elementType==QUADRILATERALS && mesh->dim==3){
    kernelInfo["defines/" "p_fainv"] = (dfloat) 0.0;
    kernelInfo["defines/" "p_invRadiusSq"] = (dfloat) 1./(mesh->sphereRadius*mesh->sphereRadius);
  }
  int NblockV = 128/mesh->Np; // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = 128/maxNodes; // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  int NblockCub = 128/mesh->cubNp; // works for CUDA
  kernelInfo["defines/" "p_NblockCub"]= NblockCub;

  // physics 
  kernelInfo["defines/" "p_Lambda2"]= 0.5f;
  kernelInfo["defines/" "p_sqrtRT"]= bns->sqrtRT;
  kernelInfo["defines/" "p_isqrtRT"]= (dfloat)(1./bns->sqrtRT);
  kernelInfo["defines/" "p_sqrt2"]= (dfloat)sqrt(2.);
  kernelInfo["defines/" "p_isq12"]= (dfloat)sqrt(1./12.);
  kernelInfo["defines/" "p_isq6"]= (dfloat)sqrt(1./6.);
  kernelInfo["defines/" "p_invsqrt2"]= (dfloat)sqrt(1./2.);
  kernelInfo["defines/" "p_tauInv"]= bns->tauInv;

  dfloat AX = 0, AY = 0, AZ = 0;

  if(options.getArgs("BODYFORCE-X", AX))
    if(AX)
      kernelInfo["defines/" "p_AX"]= AX/bns->sqrtRT;
  
  if(options.getArgs("BODYFORCE-Y", AY))
    if(AY)
      kernelInfo["defines/" "p_AY"]= AY/bns->sqrtRT;

  if(bns->dim==3){
    if(options.getArgs("BODYFORCE-Z", AZ))
      if(AZ)
        kernelInfo["defines/" "p_AZ"]= AZ/bns->sqrtRT;
  }

  kernelInfo["defines/" "p_q1bar"]= q1bar;
  kernelInfo["defines/" "p_q2bar"]= q2bar;
  kernelInfo["defines/" "p_q3bar"]= q3bar;
  kernelInfo["defines/" "p_q4bar"]= q4bar;
  kernelInfo["defines/" "p_q5bar"]= q5bar;
  kernelInfo["defines/" "p_q6bar"]= q6bar;
  
  if(bns->dim==3){
    kernelInfo["defines/" "p_q7bar"]= q7bar;
    kernelInfo["defines/" "p_q8bar"]= q8bar;
    kernelInfo["defines/" "p_q9bar"]= q9bar;
    kernelInfo["defines/" "p_q10bar"]= q10bar;
  }

  kernelInfo["defines/" "p_alpha0"]= (dfloat).01f;
  kernelInfo["defines/" "p_pmlAlpha"]= (dfloat)0.2f;
  kernelInfo["defines/" "p_blockSize"]= blockSize;
  kernelInfo["defines/" "p_NrkStages"]= bns->NrkStages;

  if(options.compareArgs("ABSORBING LAYER", "PML"))
    kernelInfo["defines/" "p_PML"]= (int) 1;
  else
    kernelInfo["defines/" "p_PML"]= (int) 0;



  if(bns->fexplicit) // full explicit or semi-analaytic
    kernelInfo["defines/" "p_SEMI_ANALYTIC"]= (int) 0;
  else
    kernelInfo["defines/" "p_SEMI_ANALYTIC"]= (int) 1;

  if(options.compareArgs("TIME INTEGRATOR", "MRSAAB"))
    kernelInfo["defines/" "p_MRSAAB"]= (int) 1;
  else
    kernelInfo["defines/" "p_MRSAAB"]= (int) 0;


  kernelInfo["defines/" "p_Nvort"]= bns->Nvort;

  if(bns->dim==3){
    kernelInfo["defines/" "p_isoNfields"]= bns->isoNfields;
    
    // Define Isosurface Area Tolerance
    kernelInfo["defines/" "p_triAreaTol"]= (dfloat) 1.0E-16;

    kernelInfo["defines/" "p_dim"]= bns->dim;
    kernelInfo["defines/" "p_plotNp"]= mesh->plotNp;
    kernelInfo["defines/" "p_plotNelements"]= mesh->plotNelements;
    
    int plotNthreads = mymax(mesh->Np, mymax(mesh->plotNp, mesh->plotNelements));
    kernelInfo["defines/" "p_plotNthreads"]= plotNthreads;
    
  } 

  // set kernel name suffix
  char *suffix, *suffixUpdate;
  
  if(bns->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(bns->elementType==QUADRILATERALS){
    if(bns->dim==2)
      suffix = strdup("Quad2D");
    else
      suffix = strdup("Quad3D");
  }
  if(bns->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(bns->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  if(bns->elementType==TRIANGLES || bns->elementType==QUADRILATERALS)
    if(bns->dim==2)
      suffixUpdate = strdup("2D");
    else
      suffixUpdate = strdup("3D");

  if(bns->elementType==TETRAHEDRA || bns->elementType==HEXAHEDRA)
    suffixUpdate = strdup("3D");
  

  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo["includes"] += (char*)boundaryHeaderFileName.c_str();

  char fileName[BUFSIZ], kernelName[BUFSIZ];
  for (int r=0;r<mesh->size;r++){

    if (r==mesh->rank) {

      // Volume kernels
      sprintf(fileName, DBNS "/okl/bnsVolume%s.okl", suffix);

      sprintf(kernelName, "bnsVolume%s", suffix);
      bns->volumeKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      if(bns->pmlFlag){
      // No that nonlinear terms are always integrated using cubature rules
      // this cubature shift is for sigma terms on pml formulation
        if(bns->pmlcubature){
          sprintf(kernelName, "bnsPmlVolumeCub%s", suffix);
          bns->pmlVolumeKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);
        }else{
          sprintf(kernelName, "bnsPmlVolume%s", suffix);
          bns->pmlVolumeKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);        
         }
      }
      
      // Relaxation kernels
      sprintf(fileName, DBNS "/okl/bnsRelaxation%s.okl", suffix);

      sprintf(kernelName, "bnsRelaxation%s", suffix);
      bns->relaxationKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      if(bns->pmlFlag){
        if(bns->pmlcubature){
          sprintf(kernelName, "bnsPmlRelaxationCub%s", suffix);        
          bns->pmlRelaxationKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);        
        }else{
          sprintf(kernelName, "bnsPmlRelaxation%s", suffix);        
          bns->pmlRelaxationKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);        
        }
      }

      
      // Surface kernels 
      sprintf(fileName, DBNS "/okl/bnsSurface%s.okl", suffix);

      if(options.compareArgs("TIME INTEGRATOR","MRSAAB")){
        sprintf(kernelName, "bnsMRSurface%s", suffix);
        bns->surfaceKernel = mesh->device.buildKernel(fileName,kernelName, kernelInfo);

        if(bns->pmlFlag){
          sprintf(kernelName, "bnsMRPmlSurface%s", suffix);
          bns->pmlSurfaceKernel = mesh->device.buildKernel(fileName,kernelName, kernelInfo);
        }
        }else{
          sprintf(kernelName, "bnsSurface%s", suffix);
          bns->surfaceKernel = mesh->device.buildKernel(fileName,kernelName, kernelInfo);

        if(bns->pmlFlag){
          sprintf(kernelName, "bnsPmlSurface%s", suffix);
          bns->pmlSurfaceKernel = mesh->device.buildKernel(fileName,kernelName, kernelInfo);
        }
      }

      
      sprintf(fileName, DBNS "/okl/bnsUpdate%s.okl", suffixUpdate);
      // Update Kernels
      if(options.compareArgs("TIME INTEGRATOR","LSERK")){
        sprintf(kernelName, "bnsLSERKUpdate%s", suffixUpdate);
        bns->updateKernel = mesh->device.buildKernel(fileName, kernelName,kernelInfo);

        if(bns->pmlFlag){
          sprintf(kernelName, "bnsLSERKPmlUpdate%s", suffixUpdate);
          bns->pmlUpdateKernel = mesh->device.buildKernel(fileName, kernelName,kernelInfo);
        }
      } else if(options.compareArgs("TIME INTEGRATOR","SARK")){
        sprintf(kernelName, "bnsSARKUpdateStage%s", suffixUpdate);
        bns->updateStageKernel = mesh->device.buildKernel(fileName,kernelName, kernelInfo);

        if(bns->pmlFlag){
          sprintf(kernelName, "bnsSARKPmlUpdateStage%s", suffixUpdate);
          bns->pmlUpdateStageKernel = mesh->device.buildKernel(fileName,kernelName, kernelInfo);
        }
  
        sprintf(kernelName, "bnsSARKUpdate%s", suffixUpdate);
        bns->updateKernel = mesh->device.buildKernel(fileName, kernelName,kernelInfo);

      if(bns->pmlFlag){
        sprintf(kernelName, "bnsSARKPmlUpdate%s", suffixUpdate);
        bns->pmlUpdateKernel = mesh->device.buildKernel(fileName, kernelName,kernelInfo);
      }
      sprintf(fileName, DBNS "/okl/bnsErrorEstimate.okl");
      sprintf(kernelName, "bnsErrorEstimate");
      bns->errorEstimateKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      } else if(options.compareArgs("TIME INTEGRATOR","MRSAAB")){
      
        sprintf(kernelName, "bnsMRSAABTraceUpdate%s", suffixUpdate);
        bns->traceUpdateKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "bnsMRSAABUpdate%s", suffixUpdate);
        bns->updateKernel = mesh->device.buildKernel(fileName, kernelName,kernelInfo);

        if(bns->pmlFlag){
          sprintf(kernelName, "bnsMRSAABPmlUpdate%s", suffixUpdate);
          bns->pmlUpdateKernel = mesh->device.buildKernel(fileName, kernelName,kernelInfo);
        }
      }

      sprintf(fileName, DBNS "/okl/bnsVorticity%s.okl",suffix);
      sprintf(kernelName, "bnsVorticity%s", suffix);
      bns->vorticityKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);


      // This needs to be unified
      mesh->haloExtractKernel =
      mesh->device.buildKernel(DHOLMES "/okl/meshHaloExtract3D.okl","meshHaloExtract3D",kernelInfo);


      if(bns->dim==3){

        mesh->gatherKernel = 
          mesh->device.buildKernel(DHOLMES "/okl/gather.okl","gather", kernelInfo);

        mesh->scatterKernel =
          mesh->device.buildKernel(DHOLMES "/okl/scatter.okl","scatter",kernelInfo);

        mesh->gatherScatterKernel =
          mesh->device.buildKernel(DHOLMES "/okl/gatherScatter.okl", "gatherScatter", kernelInfo);

        mesh->getKernel = 
          mesh->device.buildKernel(DHOLMES "/okl/get.okl", "get", kernelInfo);

        mesh->putKernel =
          mesh->device.buildKernel(DHOLMES "/okl/put.okl", "put",kernelInfo);

        bns->dotMultiplyKernel = mesh->device.buildKernel(DBNS "/okl/bnsDotMultiply.okl", "bnsDotMultiply", kernelInfo);

        // kernels from volume file
  if(bns->elementType!=QUADRILATERALS){
    sprintf(fileName, DBNS "/okl/bnsIsoSurface3D.okl");
    sprintf(kernelName, "bnsIsoSurface3D");
    
    bns->isoSurfaceKernel =
      mesh->device.buildKernel(fileName, kernelName, kernelInfo);
  }
      }
    }
    MPI_Barrier(mesh->comm);
  }



  // Setup GatherScatter
  if(bns->dim==3){
    int verbose = 1;
    dlong Ntotal = mesh->Np*mesh->Nelements;
    meshParallelGatherScatterSetup(mesh, Ntotal, mesh->globalIds, mesh->comm, verbose);
  }

  return bns; 
}




