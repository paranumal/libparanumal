#include "mppf.h"
#include "omp.h"
#include <unistd.h>

mppf_t *mppfSetup(mesh_t *mesh, setupAide options){
  
  mppf_t *mppf = (mppf_t *) calloc(1,sizeof(mppf_t));

  mppf->mesh    = mesh; 
  mppf->options = options; 

  options.getArgs("MESH DIMENSION", mppf->dim);
  options.getArgs("ELEMENT TYPE", mppf->elementType);

  // Flow Side
  mppf->NVfields = (mppf->dim==3) ? 3:2; //  Total Number of Velocity Fields
  mppf->NTfields = (mppf->dim==3) ? 4:3; // Total Velocity + Pressure

  // Phase Tracking Side
  mesh->Nfields = 1; // set mesh Nfields to 1, this may cause problems in meshOccaSetup 


  mppf->ARKswitch = 0; 

  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    mppf->extbdfA = (dfloat*) calloc(3, sizeof(dfloat));
    mppf->extbdfB = (dfloat*) calloc(3, sizeof(dfloat));
    mppf->extbdfC = (dfloat*) calloc(3, sizeof(dfloat));
    mppf->extC    = (dfloat*) calloc(3, sizeof(dfloat));
  }

  if (       options.compareArgs("TIME INTEGRATOR", "EXTBDF1")) {
    mppf->Nstages = 1;
    mppf->temporalOrder = 1;
    mppf->g0 = 1.0;
  } else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF2")) {
    mppf->Nstages = 2;
    mppf->temporalOrder = 2;
    mppf->g0 = 1.5;
  } else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF3")) {
    mppf->Nstages = 3;
    mppf->temporalOrder = 3;
    mppf->g0 = 11.f/6.f;
  }



  // Initialize fields
  dlong Nlocal = mesh->Np*mesh->Nelements;
  dlong Ntotal = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

  mppf->Ntotal = Ntotal;
  mppf->fieldOffset = Ntotal;
  mppf->Nblock = (Nlocal+blockSize-1)/blockSize;

  // flow side storage
  mppf->U       = (dfloat*) calloc(mppf->NVfields*  mppf->Nstages*Ntotal,sizeof(dfloat));
  mppf->P       = (dfloat*) calloc(                 mppf->Nstages*Ntotal,sizeof(dfloat));

  // interface side storage
  mppf->Phi     = (dfloat*) calloc(                 mppf->Nstages*Ntotal,sizeof(dfloat));
  mppf->Psi     = (dfloat*) calloc(                               Ntotal,sizeof(dfloat));
  mppf->Rho     = (dfloat*) calloc(                               Ntotal,sizeof(dfloat));
  mppf->Mu      = (dfloat*) calloc(                               Ntotal,sizeof(dfloat));
  mppf->GMu     = (dfloat*) calloc(                mppf->NVfields*Ntotal,sizeof(dfloat));
  
  // rhs storage
  mppf->rhsU    = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  mppf->rhsV    = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  mppf->rhsW    = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  mppf->rhsP    = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  mppf->rhsPhi  = (dfloat*) calloc(Ntotal,sizeof(dfloat));

  //additional field storage
  mppf->NU      = (dfloat*) calloc(mppf->NVfields*(mppf->Nstages+1)*Ntotal,sizeof(dfloat));
  mppf->LU      = (dfloat*) calloc(mppf->NVfields*(mppf->Nstages+1)*Ntotal,sizeof(dfloat));
  mppf->GP      = (dfloat*) calloc(mppf->NVfields*(mppf->Nstages+1)*Ntotal,sizeof(dfloat));
  mppf->NPhi    = (dfloat*) calloc(               (mppf->Nstages+1)*Ntotal,sizeof(dfloat));
  mppf->HPhi    = (dfloat*) calloc(               (mppf->Nstages+1)*Ntotal,sizeof(dfloat));
  
  // This needs to be changed too much storage!!!!!!!
  mppf->GU      = (dfloat*) calloc(mppf->NVfields*mppf->NVfields*Ntotal,sizeof(dfloat));

  mppf->GPhi   =  (dfloat*) calloc(mppf->NVfields*Ntotal, sizeof(dfloat));

  mppf->rkU     = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
  mppf->rkP     = (dfloat*) calloc(               Ntotal,sizeof(dfloat));
  mppf->rkPhi   = (dfloat*) calloc(               Ntotal,sizeof(dfloat));
  mppf->PI      = (dfloat*) calloc(               Ntotal,sizeof(dfloat));

  mppf->rkNU    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
  mppf->rkLU    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
  mppf->rkGP    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));

  //plotting fields
  mppf->Vort    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
  mppf->Div     = (dfloat*) calloc(               Nlocal,sizeof(dfloat));
  mppf->Ue      = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));

  //  Substepping will come  here // NOT IN USE CURRENTLY
  if(mppf->elementType==HEXAHEDRA)
    mppf->cU = (dfloat *) calloc(mppf->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
  else 
    mppf->cU = mppf->U;

  mppf->Nsubsteps = 0;
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF"))
    options.getArgs("SUBCYCLING STEPS",mppf->Nsubsteps);
    
    
  // if(mppf->Nsubsteps){
  //   mppf->Ud    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
  //   mppf->Ue    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
  //   mppf->resU  = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
  //   mppf->rhsUd = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));

  //   if(mppf->elementType==HEXAHEDRA)
  //     mppf->cUd = (dfloat *) calloc(mppf->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
  //   else 
  //     mppf->cUd = mppf->U;
  // }



  // No gravitational acceleration
  dfloat g[3]; g[0] = 0.0; g[1] = 0.0; g[2] = 0.0;  
  
  // Define mean-flow
  options.getArgs("UBAR", mppf->ubar);
  options.getArgs("VBAR", mppf->vbar);
  if (mppf->dim==3)
    options.getArgs("WBAR", mppf->wbar);
  options.getArgs("PBAR", mppf->pbar);

  // Define phase properties
  options.getArgs("VISCOSITY PHASE 1", mppf->mu1);
  options.getArgs("VISCOSITY PHASE 2", mppf->mu2);

  options.getArgs("DENSITY PHASE 1", mppf->rho1);
  options.getArgs("DENSITY PHASE 2", mppf->rho2);
  
  
  //Reynolds number defined on the phase with lower material properties
  int phase_check = 1; 
  phase_check = (mppf->mu1 > mppf->mu2) ? 0: 1; 
  if(phase_check == 0){ printf("Phase 1 viscosity has to be smaller than Phase 2 \n"); exit(EXIT_FAILURE);} 

  phase_check = (mppf->rho1 > mppf->rho2) ? 0: 1; 
  if(phase_check == 0){ printf("Phase 1 density has to be smaller than Phase 2 \n"); exit(EXIT_FAILURE);} 


  mppf->Re = mppf->ubar*mppf->rho1/mppf->mu1;

  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  if(mppf->dim==3)
    meshOccaSetup3D(mesh, options, kernelInfo);
  else
    meshOccaSetup2D(mesh, options, kernelInfo);

  
  occa::properties kernelInfoV    = kernelInfo;
  occa::properties kernelInfoP    = kernelInfo;
  occa::properties kernelInfoPhi  = kernelInfo;
  occa::properties kernelInfoPsi  = kernelInfo;

  // ADD-DEFINES
  kernelInfo["defines/" "p_pbar"]= mppf->pbar;
  kernelInfo["defines/" "p_ubar"]= mppf->ubar;
  kernelInfo["defines/" "p_vbar"]= mppf->vbar;
  kernelInfo["defines/" "p_wbar"]= mppf->wbar;

  // Multiphase
  kernelInfo["defines/" "p_mu1"] = mppf->mu1;
  kernelInfo["defines/" "p_mu2"] = mppf->mu2;
  kernelInfo["defines/" "p_rho1"] = mppf->rho1;
  kernelInfo["defines/" "p_rho2"] = mppf->rho2;

  mppf->rho0 = mymin(mppf->rho1, mppf->rho2);
  mppf->nu0  = 0.5*mymax(mppf->mu1, mppf->mu2)/ mymin(mppf->rho1, mppf->rho2); 

  kernelInfo["defines/" "p_invrho0"] = 1.0/mppf->rho0;
  kernelInfo["defines/" "p_rho0"] = mppf->rho0;
  kernelInfo["defines/" "p_nu0"] = mppf->nu0;

  kernelInfo["defines/" "p_NTfields"]= mppf->NTfields;
  kernelInfo["defines/" "p_NVfields"]= mppf->NVfields;
  kernelInfo["defines/" "p_NfacesNfp"]=  mesh->Nfaces*mesh->Nfp;
  kernelInfo["defines/" "p_Nstages"]=  mppf->Nstages;
  kernelInfo["defines/" "p_SUBCYCLING"]=  mppf->Nsubsteps;


  //add boundary data to kernel info
  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo["includes"] += (char*)boundaryHeaderFileName.c_str();

  mppf->o_U   = mesh->device.malloc(mppf->NVfields*mppf->Nstages*Ntotal*sizeof(dfloat), mppf->U);
  mppf->o_P   = mesh->device.malloc(               mppf->Nstages*Ntotal*sizeof(dfloat), mppf->P);
  mppf->o_Phi = mesh->device.malloc(               mppf->Nstages*Ntotal*sizeof(dfloat), mppf->Phi);
  mppf->o_Rho = mesh->device.malloc(                             Ntotal*sizeof(dfloat), mppf->Rho);
  mppf->o_Mu  = mesh->device.malloc(                             Ntotal*sizeof(dfloat), mppf->Mu);
  mppf->o_GMu = mesh->device.malloc(              mppf->NVfields*Ntotal*sizeof(dfloat), mppf->GMu);
  mppf->o_Ue  = mesh->device.malloc(              mppf->NVfields*Ntotal*sizeof(dfloat), mppf->Ue);
  
   for (int r=0;r<mesh->size;r++) {
    if (r==mesh->rank) {
      if (mppf->dim==2){ 
        mppf->setFlowFieldKernel =  mesh->device.buildKernel(DMPPF "/okl/mppfSetFlowField2D.okl", "mppfSetFlowField2D", kernelInfo);  
        mppf->setPhaseFieldKernel =  mesh->device.buildKernel(DMPPF "/okl/mppfSetPhaseField2D.okl", "mppfSetPhaseField2D", kernelInfo);  
        mppf->setMaterialPropertyKernel =  mesh->device.buildKernel(DMPPF "/okl/mppfSetMaterialProperty2D.okl", "mppfSetMaterialProperty2D", kernelInfo);
      }else{
        mppf->setFlowFieldKernel =  mesh->device.buildKernel(DMPPF "/okl/mppfSetFlowField3D.okl", "mppfSetFlowField3D", kernelInfo);  
        mppf->setPhaseFieldKernel =  mesh->device.buildKernel(DMPPF "/okl/mppfSetPhaseField3D.okl", "mppfSetPhaseField3D", kernelInfo);  
        mppf->setMaterialPropertyKernel =  mesh->device.buildKernel(DMPPF "/okl/mppfSetMaterialProperty3D.okl", "mppfSetMaterialProperty3D", kernelInfo);
      }
    }
    MPI_Barrier(mesh->comm);
  }

  mppf->startTime =0.0;
  options.getArgs("START TIME", mppf->startTime);
  mppf->setFlowFieldKernel(mesh->Nelements,
                          mppf->startTime,
                          mesh->o_x,
                          mesh->o_y,
                          mesh->o_z,
                          mppf->fieldOffset,
                          mppf->o_U,
                          mppf->o_P);
  
  mppf->o_U.copyTo(mppf->U);


  

// Set interface thickness and  time step size
  // Find max initial velocity and minum mesh thickness
  dfloat hmin = 1e9, hmax = 0, umax = 0;
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      dfloat hest = 2./(sJ*invJ);

      hmin = mymin(hmin, hest);
      hmax = mymax(hmax, hest);
    }

    for(int n=0;n<mesh->Np;++n){
      const dlong id = n + mesh->Np*e;
      dfloat t = 0;
      dfloat uxn = mppf->U[id+0*mppf->fieldOffset];
      dfloat uyn = mppf->U[id+1*mppf->fieldOffset];
      dfloat uzn = 0.0;
      if (mppf->dim==3) uzn = mppf->U[id+2*mppf->fieldOffset];
      //Squared maximum velocity
      dfloat numax;
      if (mppf->dim==2) numax = uxn*uxn + uyn*uyn;
      else             numax = uxn*uxn + uyn*uyn + uzn*uzn;
      umax = mymax(umax, numax);
    }
  }

  // Maximum Velocity
  umax = sqrt(umax);
  dfloat magVel = mymax(umax,1.0); // Correction for initial zero velocity

  options.getArgs("CFL", mppf->cfl);
  dfloat dt     = mppf->cfl* hmin/( (mesh->N+1.)*(mesh->N+1.) * magVel) ;

  // MPI_Allreduce to get global minimum dt and hmin
  MPI_Allreduce(&dt  , &(mppf->dti),  1, MPI_DFLOAT, MPI_MIN, mesh->comm);
  MPI_Allreduce(&hmin, &(mppf->hmin), 1, MPI_DFLOAT, MPI_MIN, mesh->comm);
  
  // characteristic length of interface thickness
  mppf->eta =  mppf->hmin; // Change this later, currently no mppfide !!!!!

  mppf->eta = 0.1;  // WARNING MUST COMMENT OUT !!!!!!!!!!!!!!!!!!!!!

  mppf->dt = mppf->dti;

  options.getArgs("FINAL TIME", mppf->finalTime);
  options.getArgs("START TIME", mppf->startTime);
  
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    mppf->NtimeSteps = (mppf->finalTime-mppf->startTime)/mppf->dt;

    if(mppf->Nsubsteps){
      mppf->dt         = mppf->Nsubsteps*mppf->dt;
      mppf->NtimeSteps = (mppf->finalTime-mppf->startTime)/mppf->dt;
      mppf->dt         = (mppf->finalTime-mppf->startTime)/mppf->NtimeSteps;
      mppf->sdt        = mppf->dt/mppf->Nsubsteps;
    } else{
      mppf->NtimeSteps = (mppf->finalTime-mppf->startTime)/mppf->dt;
      mppf->dt         = (mppf->finalTime-mppf->startTime)/mppf->NtimeSteps;
    }
  }

  mppf->dtMIN = 1E-2*mppf->dt; //minumum allowed timestep
  if (mppf->Nsubsteps && mesh->rank==0) 
    printf("dt: %.8f and sdt: %.8f ratio: %.8f \n", mppf->dt, mppf->sdt, mppf->dt/mppf->sdt);


  // // Set pahse field function on device
  mppf->setPhaseFieldKernel(mesh->Nelements,
                          mppf->startTime,
                          mppf->eta,
                          mesh->o_x,
                          mesh->o_y,
                          mesh->o_z,
                          mppf->o_Phi);


  // Smooth density and viscosity on device 
  mppf->setMaterialPropertyKernel(mesh->Nelements,
                                  mppf->o_Phi,
                                  mppf->o_Rho,
                                  mppf->o_Mu);
  
#if 1
  mppf->o_Rho.copyTo(mppf->Rho);
  mppf->o_Mu.copyTo(mppf->Mu);
  mppf->o_P.copyTo(mppf->P);
  mppf->o_Phi.copyTo(mppf->Phi);
  char fname[BUFSIZ];
  string outName;
  mppf->options.getArgs("OUTPUT FILE NAME", outName);
  sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(), mesh->rank, mppf->frame++);

  mppfPlotVTU(mppf, fname);
#endif


  mppf->outputStep = 0;
  options.getArgs("TSTEPS FOR SOLUTION OUTPUT", mppf->outputStep);
  // if (mesh->rank==0) printf("Nsteps = %d NerrStep= %d dt = %.8e\n", mppf->NtimeSteps,mppf->outputStep, mppf->dt);

  mppf->outputForceStep = 0;
  options.getArgs("TSTEPS FOR FORCE OUTPUT", mppf->outputForceStep);

  options.getArgs("MIXING ENERGY DENSITY", mppf->chL);
  options.getArgs("MOBILITY", mppf->chM);


  kernelInfo["defines/" "p_chL"]        = mppf->chL;   // mixing energy
  kernelInfo["defines/" "p_chM"]        = mppf->chM;   // mobility
  kernelInfo["defines/" "p_chInvLM"]    = 1.0/ (mppf->chL*mppf->chM);   // mobility

  // Some derived parameters
  mppf->idt     = 1.0/mppf->dt;
  mppf->eta2    = mppf->eta*mppf->eta; 
  mppf->inveta2 = 1.0/ mppf->eta2; 

  mppf->factorS = 1.5; // has to be >1.0 

  mppf->chS = mppf->factorS*mppf->eta2*sqrt(4.0*mppf->g0/ (mppf->chM*mppf->chL*mppf->dt));   
  mppf->chA  = -mppf->chS/(2.0*mppf->eta2) * (1.0 - sqrt(1 - 4.0*mppf->g0*mppf->eta2*mppf->eta2/(mppf->chM*mppf->chL*mppf->dt*mppf->chS*mppf->chS)));   
  
  // Hold S/n^2
  mppf->chSeta2 = mppf->chS/mppf->eta2; 
  
  
  // Helmholtz solve lambda's i.e. -laplace*psi + [alpha+ S/eta^2]*psi = -Q 
  mppf->lambdaPsi = mppf->chA + mppf->chS*mppf->inveta2;
  // Helmholtz solve lambda's i.e. -laplace*phi +[-alpha]*phi = -psi 
  mppf->lambdaPhi = -mppf->chA;

  mppf->phiOptions = options;
  mppf->phiOptions.setArgs("KRYLOV SOLVER",        options.getArgs("PHASE FIELD KRYLOV SOLVER"));
  mppf->phiOptions.setArgs("DISCRETIZATION",       options.getArgs("PHASE FIELD DISCRETIZATION"));
  mppf->phiOptions.setArgs("BASIS",                options.getArgs("PHASE FIELD BASIS"));
  mppf->phiOptions.setArgs("PRECONDITIONER",       options.getArgs("PHASE FIELD PRECONDITIONER"));
  mppf->phiOptions.setArgs("MULTIGRID COARSENING", options.getArgs("PHASE FIELD MULTIGRID COARSENING"));
  mppf->phiOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("PHASE FIELD MULTIGRID SMOOTHER"));
  mppf->phiOptions.setArgs("PARALMOND CYCLE",      options.getArgs("PHASE FIELD PARALMOND CYCLE"));
  mppf->phiOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("PHASE FIELD PARALMOND SMOOTHER"));
  mppf->phiOptions.setArgs("PARALMOND PARTITION",  options.getArgs("PHASE FIELD PARALMOND PARTITION"));

  mppf->vOptions = options;
  mppf->vOptions.setArgs("KRYLOV SOLVER",        options.getArgs("VELOCITY KRYLOV SOLVER"));
  mppf->vOptions.setArgs("DISCRETIZATION",       options.getArgs("VELOCITY DISCRETIZATION"));
  mppf->vOptions.setArgs("BASIS",                options.getArgs("VELOCITY BASIS"));
  mppf->vOptions.setArgs("PRECONDITIONER",       options.getArgs("VELOCITY PRECONDITIONER"));
  mppf->vOptions.setArgs("MULTIGRID COARSENING", options.getArgs("VELOCITY MULTIGRID COARSENING"));
  mppf->vOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("VELOCITY MULTIGRID SMOOTHER"));
  mppf->vOptions.setArgs("PARALMOND CYCLE",      options.getArgs("VELOCITY PARALMOND CYCLE"));
  mppf->vOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("VELOCITY PARALMOND SMOOTHER"));
  mppf->vOptions.setArgs("PARALMOND PARTITION",  options.getArgs("VELOCITY PARALMOND PARTITION"));

  mppf->pOptions = options;
  mppf->pOptions.setArgs("KRYLOV SOLVER",        options.getArgs("PRESSURE KRYLOV SOLVER"));
  mppf->pOptions.setArgs("DISCRETIZATION",       options.getArgs("PRESSURE DISCRETIZATION"));
  mppf->pOptions.setArgs("BASIS",                options.getArgs("PRESSURE BASIS"));
  mppf->pOptions.setArgs("PRECONDITIONER",       options.getArgs("PRESSURE PRECONDITIONER"));
  mppf->pOptions.setArgs("MULTIGRID COARSENING", options.getArgs("PRESSURE MULTIGRID COARSENING"));
  mppf->pOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("PRESSURE MULTIGRID SMOOTHER"));
  mppf->pOptions.setArgs("PARALMOND CYCLE",      options.getArgs("PRESSURE PARALMOND CYCLE"));
  mppf->pOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("PRESSURE PARALMOND SMOOTHER"));
  mppf->pOptions.setArgs("PARALMOND PARTITION",  options.getArgs("PRESSURE PARALMOND PARTITION"));

  if (mesh->rank==0) printf("==================ELLIPTIC SOLVE SETUP=========================\n");

  // SetUp Boundary Flags types for Elliptic Solve
  // bc = 1 -> wall
  // bc = 2 -> inflow
  // bc = 3 -> outflow
  // bc = 4 -> x-aligned slip
  // bc = 5 -> y-aligned slip
  // bc = 6 -> z-aligned slip
  int uBCType[7]   = {0,1,1,2,1,2,2}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int vBCType[7]   = {0,1,1,2,2,1,2}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int wBCType[7]   = {0,1,1,2,2,2,1}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int pBCType[7]   = {0,2,2,1,2,2,2}; // bc=3 => outflow => Dirichlet => pBCType[3] = 1, etc.
  // int phiBCType[7] = {0,2,2,2,2,2,2}; // All homogenous Neumann BCs for Phi and Psi solves 
  int phiBCType[7] = {0,2,2,2,2,2,2}; // All homogenous Neumann BCs for Phi and Psi solves 

  //Solver tolerances 
  mppf->presTOL = 1E-8;
  mppf->velTOL  = 1E-8;
  mppf->phiTOL  = 1E-8;

  // Use third Order Velocity Solve: full rank should converge for low orders
  if (mesh->rank==0) printf("================PHASE-FIELD SOLVE SETUP=========================\n");
 
  if (mesh->rank==0) printf("==================Phi Solve Setup=========================\n");
  // -laplace(Phi) + [-alpha ]*Phi = -Psi 
  mppf->phiSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  mppf->phiSolver->mesh = mesh;
  mppf->phiSolver->options = mppf->phiOptions;
  mppf->phiSolver->dim = mppf->dim;
  mppf->phiSolver->elementType = mppf->elementType;
  mppf->phiSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(mppf->phiSolver->BCType,phiBCType,7*sizeof(int));
  ellipticSolveSetup(mppf->phiSolver, mppf->lambdaPhi, kernelInfoPhi);

  if (mesh->rank==0) printf("==================Psi Solve Setup=========================\n");
  // -laplace(Psi) + [alpha+ S/eta^2]*Psi = -Q 
  mppf->psiSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  mppf->psiSolver->mesh = mesh;
  mppf->psiSolver->options = mppf->phiOptions; // Using the same options with Phi solver
  mppf->psiSolver->dim = mppf->dim;
  mppf->psiSolver->elementType = mppf->elementType;
  mppf->psiSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(mppf->phiSolver->BCType,phiBCType,7*sizeof(int)); // Using the same boundary flags
  ellipticSolveSetup(mppf->psiSolver, mppf->lambdaPsi, kernelInfoPsi);


  if (mesh->rank==0) printf("==================PRESSURE SOLVE SETUP=========================\n");
  mppf->pSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  mppf->pSolver->mesh = mesh;
  mppf->pSolver->options = mppf->pOptions;
  mppf->pSolver->dim = mppf->dim;
  mppf->pSolver->elementType = mppf->elementType;
  mppf->pSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(mppf->pSolver->BCType,pBCType,7*sizeof(int));
  ellipticSolveSetup(mppf->pSolver, 0.0, kernelInfoP);





  //make node-wise boundary flags
  mppf->VmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  mppf->PmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  for (int e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) mppf->VmapB[n+e*mesh->Np] = 1E9;
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[f+e*mesh->Nfaces];
      if (bc>0) {
        for (int n=0;n<mesh->Nfp;n++) {
          int fid = mesh->faceNodes[n+f*mesh->Nfp];
          mppf->VmapB[fid+e*mesh->Np] = mymin(bc,mppf->VmapB[fid+e*mesh->Np]);
          mppf->PmapB[fid+e*mesh->Np] = mymax(bc,mppf->PmapB[fid+e*mesh->Np]);
        }
      }
    }
  }
  gsParallelGatherScatter(mesh->hostGsh, mppf->VmapB, "int", "min"); 
  gsParallelGatherScatter(mesh->hostGsh, mppf->PmapB, "int", "max"); 

  for (int n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (mppf->VmapB[n] == 1E9) {
      mppf->VmapB[n] = 0.;
    }
  }
  mppf->o_VmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), mppf->VmapB);
  mppf->o_PmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), mppf->PmapB);
  
  // Kernel defines
  kernelInfo["defines/" "p_blockSize"]= blockSize;
  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  kernelInfo["defines/" "p_blockSize"]= blockSize;
  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

   // if(options.compareArgs("TIME INTEGRATOR", "EXTBDF"))
  kernelInfo["defines/" "p_EXTBDF"]= 1;
  // else
    // kernelInfo["defines/" "p_EXTBDF"]= 0;

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = 256/mesh->Np; // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = 256/maxNodes; // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;


  int maxNodesVolumeCub = mymax(mesh->cubNp,mesh->Np);  
  kernelInfo["defines/" "p_maxNodesVolumeCub"]= maxNodesVolumeCub;
  int cubNblockV = mymax(1,256/maxNodesVolumeCub);
  //
  int maxNodesSurfaceCub = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  kernelInfo["defines/" "p_maxNodesSurfaceCub"]=maxNodesSurfaceCub;
  int cubNblockS = mymax(256/maxNodesSurfaceCub,1); // works for CUDA
  //
  kernelInfo["defines/" "p_cubNblockV"]=cubNblockV;
  kernelInfo["defines/" "p_cubNblockS"]=cubNblockS;

  
  // IsoSurface related
  if(mppf->dim==3){
    if(options.compareArgs("OUTPUT FILE FORMAT", "ISO")){
      kernelInfo["defines/" "p_isoNfields"]= mppf->isoNfields;
      // Define Isosurface Area Tolerance
      kernelInfo["defines/" "p_triAreaTol"]= (dfloat) 1.0E-16;

      kernelInfo["defines/" "p_dim"]= mppf->dim;
      kernelInfo["defines/" "p_plotNp"]= mesh->plotNp;
      kernelInfo["defines/" "p_plotNelements"]= mesh->plotNelements;
      
      int plotNthreads = mymax(mesh->Np, mymax(mesh->plotNp, mesh->plotNelements));
      kernelInfo["defines/" "p_plotNthreads"]= plotNthreads;     
    }
 } 


 if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    dfloat rkC[4] = {1.0, 0.0, -1.0, -2.0};

    mppf->o_rkC  = mesh->device.malloc(4*sizeof(dfloat),rkC);
    mppf->o_extbdfA = mesh->device.malloc(3*sizeof(dfloat));
    mppf->o_extbdfB = mesh->device.malloc(3*sizeof(dfloat));
    mppf->o_extbdfC = mesh->device.malloc(3*sizeof(dfloat)); 

    mppf->o_extC = mesh->device.malloc(3*sizeof(dfloat)); 

    mppf->o_prkA = mppf->o_extbdfC;
    mppf->o_prkB = mppf->o_extbdfC;
  }

  // MEMORY ALLOCATION
  mppf->o_Psi       = mesh->device.malloc(Ntotal*sizeof(dfloat), mppf->Psi);
  mppf->o_lapPhi    = mesh->device.malloc(Ntotal*sizeof(dfloat), mppf->Psi);

  mppf->o_rhsU    = mesh->device.malloc(Ntotal*sizeof(dfloat), mppf->rhsU);
  mppf->o_rhsV    = mesh->device.malloc(Ntotal*sizeof(dfloat), mppf->rhsV);
  mppf->o_rhsW    = mesh->device.malloc(Ntotal*sizeof(dfloat), mppf->rhsW);
  mppf->o_rhsP    = mesh->device.malloc(Ntotal*sizeof(dfloat), mppf->rhsP);
  mppf->o_rhsPhi  = mesh->device.malloc(Ntotal*sizeof(dfloat), mppf->rhsPhi);

  mppf->o_GPhi  = mesh->device.malloc(                  mppf->NVfields*Ntotal*sizeof(dfloat), mppf->GPhi);
  mppf->o_NPhi  = mesh->device.malloc(               (mppf->Nstages+1)*Ntotal*sizeof(dfloat), mppf->NPhi);
  mppf->o_HPhi  = mesh->device.malloc(               (mppf->Nstages+1)*Ntotal*sizeof(dfloat), mppf->HPhi);
  
  mppf->o_NU    = mesh->device.malloc(mppf->NVfields*(mppf->Nstages+1)*Ntotal*sizeof(dfloat), mppf->NU);
  mppf->o_LU    = mesh->device.malloc(mppf->NVfields*(mppf->Nstages+1)*Ntotal*sizeof(dfloat), mppf->LU);
  mppf->o_GP    = mesh->device.malloc(mppf->NVfields*(mppf->Nstages+1)*Ntotal*sizeof(dfloat), mppf->GP);
  
  // This needs to be changed, too much storage
  mppf->o_GU    = mesh->device.malloc(mppf->NVfields*mppf->NVfields*Ntotal*sizeof(dfloat), mppf->GU);
  
  mppf->o_rkU   = mesh->device.malloc(mppf->NVfields*Ntotal*sizeof(dfloat), mppf->rkU);
  mppf->o_rkP   = mesh->device.malloc(              Ntotal*sizeof(dfloat), mppf->rkP);
  mppf->o_rkPhi = mesh->device.malloc(              Ntotal*sizeof(dfloat), mppf->rkPhi);
  mppf->o_PI    = mesh->device.malloc(              Ntotal*sizeof(dfloat), mppf->PI);
  
  mppf->o_rkNU  = mesh->device.malloc(mppf->NVfields*Ntotal*sizeof(dfloat), mppf->rkNU);
  mppf->o_rkLU  = mesh->device.malloc(mppf->NVfields*Ntotal*sizeof(dfloat), mppf->rkLU);
  mppf->o_rkGP  = mesh->device.malloc(mppf->NVfields*Ntotal*sizeof(dfloat), mppf->rkGP);

  //storage for helmholtz solves
  mppf->o_UH   = mesh->device.malloc(Ntotal*sizeof(dfloat));
  mppf->o_VH   = mesh->device.malloc(Ntotal*sizeof(dfloat));
  mppf->o_WH   = mesh->device.malloc(Ntotal*sizeof(dfloat));
  // mppf->o_PhiH = mesh->device.malloc(Ntotal*sizeof(dfloat));

  //plotting fields
  mppf->o_Vort = mesh->device.malloc(mppf->NVfields*Ntotal*sizeof(dfloat), mppf->Vort);
  mppf->o_Div  = mesh->device.malloc(              Nlocal*sizeof(dfloat), mppf->Div);

  if(mppf->elementType==HEXAHEDRA)
    mppf->o_cU = mesh->device.malloc(mppf->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), mppf->cU);
  else 
    mppf->o_cU = mppf->o_U;

  if(mesh->totalHaloPairs){//halo setup
    dlong vHaloBytes    = mesh->totalHaloPairs*mesh->Np*(mppf->NVfields)*sizeof(dfloat);
    dlong pHaloBytes    = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
    dlong phiHaloBytes  = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
    dlong vGatherBytes  = mppf->NVfields*mesh->ogs->NhaloGather*sizeof(dfloat);
    
    mppf->o_vHaloBuffer   = mesh->device.malloc(vHaloBytes);
    mppf->o_pHaloBuffer   = mesh->device.malloc(pHaloBytes);
    mppf->o_phiHaloBuffer = mesh->device.malloc(phiHaloBytes);

    occa::memory o_vSendBuffer,   o_vRecvBuffer; 
    occa::memory o_pSendBuffer,   o_pRecvBuffer;
    occa::memory o_phiSendBuffer, o_phiRecvBuffer;
    occa::memory o_gatherTmpPinned;

    mppf->vSendBuffer   = (dfloat*) occaHostMallocPinned(mesh->device, vHaloBytes, NULL, mppf->o_vSendBuffer);
    mppf->vRecvBuffer   = (dfloat*) occaHostMallocPinned(mesh->device, vHaloBytes, NULL, mppf->o_vRecvBuffer);

    mppf->pSendBuffer   = (dfloat*) occaHostMallocPinned(mesh->device, pHaloBytes, NULL, mppf->o_pSendBuffer);
    mppf->pRecvBuffer   = (dfloat*) occaHostMallocPinned(mesh->device, pHaloBytes, NULL, mppf->o_pRecvBuffer);

    mppf->phiSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, phiHaloBytes, NULL, mppf->o_phiSendBuffer);
    mppf->phiRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, phiHaloBytes, NULL, mppf->o_phiRecvBuffer);

    mppf->velocityHaloGatherTmp = (dfloat*) occaHostMallocPinned(mesh->device, vGatherBytes, NULL, mppf->o_gatherTmpPinned);
    
    mppf->o_velocityHaloGatherTmp = mesh->device.malloc(vGatherBytes,  mppf->velocityHaloGatherTmp);
  }


  // set kernel name suffix
  char *suffix;
  
  if(mppf->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(mppf->elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(mppf->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(mppf->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  for (int r=0;r<mesh->size;r++) {
    if (r==mesh->rank) {

      sprintf(fileName, DMPPF "/okl/mppfHaloExchange.okl");
      sprintf(kernelName, "mppfVelocityHaloExtract");
      mppf->velocityHaloExtractKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfVelocityHaloScatter");
      mppf->velocityHaloScatterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfPressureHaloExtract");
      mppf->pressureHaloExtractKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfPressureHaloScatter");
      mppf->pressureHaloScatterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfPhaseFieldHaloExtract");
      mppf->phaseFieldHaloExtractKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfPhaseFieldHaloScatter");
      mppf->phaseFieldHaloScatterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

       // ===========================================================================
      printf("Compiling Kernels\n");
      sprintf(fileName, DMPPF "/okl/mppfPhaseFieldAdvection%s.okl", suffix);
      sprintf(kernelName, "mppfPhaseFieldAdvectionCubatureVolume%s", suffix);
      mppf->phaseFieldAdvectionVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfPhaseFieldAdvectionCubatureSurface%s", suffix);
      mppf->phaseFieldAdvectionSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      // ===========================================================================//
      sprintf(fileName, DMPPF "/okl/mppfPhaseFieldDivGrad%s.okl", suffix);
      sprintf(kernelName, "mppfPhaseFieldDivGrad%s", suffix);
      mppf->phaseFieldDivGradKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // ===========================================================================//

      sprintf(fileName, DMPPF "/okl/mppfPhaseFieldRhs%s.okl", suffix);
      sprintf(kernelName, "mppfPhaseFieldRhsSolve1%s", suffix);
      mppf->phaseFieldRhsSolve1Kernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfPhaseFieldRhsSolve2%s", suffix);
      mppf->phaseFieldRhsSolve2Kernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(fileName, DMPPF "/okl/mppfPhaseFieldBC%s.okl", suffix);
      sprintf(kernelName, "mppfPhaseFieldIpdgBC%s", suffix);
      mppf->phaseFieldRhsIpdgBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      //===========================================================================//

      sprintf(fileName, DMPPF "/okl/mppfVorticity%s.okl", suffix);
      sprintf(kernelName, "mppfVorticity%s", suffix);
      mppf->vorticityKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      //===========================================================================//
      
      sprintf(fileName, DMPPF "/okl/mppfDivergence%s.okl", suffix);
      sprintf(kernelName, "mppfDivergenceVolume%s", suffix);
      mppf->divergenceVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfDivergenceSurface%s", suffix);
      mppf->divergenceSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);


      //===========================================================================//

      sprintf(fileName, DMPPF "/okl/mppfPhaseFieldGradient%s.okl", suffix);
      sprintf(kernelName, "mppfPhaseFieldGradientVolume%s", suffix);
      mppf->phaseFieldGradientVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfPhaseFieldGradientSurface%s", suffix);
      mppf->phaseFieldGradientSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // ===========================================================================

      sprintf(fileName, DMPPF "/okl/mppfAdvection%s.okl", suffix);
      sprintf(kernelName, "mppfAdvectionCubatureVolume%s", suffix);
      mppf->advectionCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfAdvectionCubatureSurface%s", suffix);
      mppf->advectionCubatureSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfAdvectionVolume%s", suffix);
      mppf->advectionVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfAdvectionSurface%s", suffix);
      mppf->advectionSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);


       // ===========================================================================

      sprintf(fileName, DMPPF "/okl/mppfAdvectionUpdate%s.okl", suffix);
      sprintf(kernelName, "mppfAdvectionUpdate%s", suffix);
      mppf->advectionUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // ===========================================================================

      sprintf(fileName, DMPPF "/okl/mppfPressureGradient%s.okl", suffix);
      sprintf(kernelName, "mppfPressureGradientVolume%s", suffix);
      mppf->pressureGradientVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfPressureGradientSurface%s", suffix);
      mppf->pressureGradientSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // ===========================================================================

      sprintf(fileName, DMPPF "/okl/mppfVelocityGradient%s.okl", suffix);
      sprintf(kernelName, "mppfVelocityGradientVolume%s", suffix);
      mppf->velocityGradientVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfVelocityGradientSurface%s", suffix);
      mppf->velocityGradientSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

       // ===========================================================================

      sprintf(fileName, DMPPF "/okl/mppfVelocityExtrapolate.okl");
      sprintf(kernelName, "mppfVelocityExtrapolate");
      mppf->velocityExtrapolateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

       // ===========================================================================

      sprintf(fileName, DMPPF "/okl/mppfPressureRhs%s.okl", suffix);
      sprintf(kernelName, "mppfPressureRhs%s", suffix);
      mppf->pressureRhsKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(fileName, DMPPF "/okl/mppfPressureBC%s.okl", suffix);
      sprintf(kernelName, "mppfPressureIpdgBC%s", suffix);
      mppf->pressureRhsIpdgBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfPressureBC%s", suffix);
      mppf->pressureRhsBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "mppfPressureAddBC%s", suffix);
      mppf->pressureAddBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      
    }
   MPI_Barrier(mesh->comm);
  }

 
if(mesh->rank==0){
    printf("=============WRITING INPUT PARAMETERS===================\n");

    printf("INTERFACE LENGTH\t:\t%.2e\n", mppf->eta);
    printf("INTERFACE THICKNESS\t:\t%.2e\n", mppf->hmin);
    printf("MINUM TIME STEP SIZE\t:\t%.2e\n", mppf->dt);
    printf("# TIME STEP\t:\t%d\n", mppf->NtimeSteps);
    printf("# SUBSTEPS\t\t:\t%d\n", mppf->Nsubsteps);
    printf("# chSeta2\t:\t%.4e\n", mppf->chSeta2);
    printf("# eta2\t\t:\t%.4e\n", mppf->eta2);



 printf("============================================================\n");
    printf("VISCOSITY PHASE 1\t:\t%.2e\n", mppf->mu1);
    printf("VISCOSITY PHASE 2\t:\t%.2e\n", mppf->mu2);
    printf("DENSITY PHASE 1\t\t:\t%.2e\n", mppf->rho1);
    printf("DENSITY PHASE 2\t\t:\t%.2e\n", mppf->rho2);
    printf("MIXING ENERGY\t\t:\t%.2e\n", mppf->chL);
    printf("MOBILITY\t\t:\t%.2e\n", mppf->chM);
    printf("CH HELMHOLTZ LAMBDA PSI\t:\t%.2e\n", mppf->lambdaPsi);
    printf("CH HELMHOLTZ LAMBDA PHI\t:\t%.2e\n", mppf->lambdaPhi);
 printf("============================================================\n");
    printf("# TIME STEPS\t\t:\t%d\n", mppf->NtimeSteps);
    printf("# OUTPUT STEP\t\t:\t%d\n", mppf->outputStep);
    printf("# FORCE STEP\t\t:\t%d\n", mppf->outputForceStep);
 printf("============================================================\n");

  }










return mppf;
}