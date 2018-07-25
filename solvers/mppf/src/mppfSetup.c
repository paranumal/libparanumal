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

  mppf->GU      = (dfloat*) calloc(mppf->NVfields*Ntotal*4,sizeof(dfloat));

  mppf->rkU     = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
  mppf->rkP     = (dfloat*) calloc(               Ntotal,sizeof(dfloat));
  mppf->PI      = (dfloat*) calloc(               Ntotal,sizeof(dfloat));

  mppf->rkNU    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
  mppf->rkLU    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
  mppf->rkGP    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));

  //plotting fields
  mppf->Vort    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
  mppf->Div     = (dfloat*) calloc(               Nlocal,sizeof(dfloat));

  //  Substepping will come  here // NOT IN USE CURRENTLY
  if(mppf->elementType==HEXAHEDRA)
    mppf->cU = (dfloat *) calloc(mppf->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
  else 
    mppf->cU = mppf->U;

  mppf->Nsubsteps = 0;
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF"))
    options.getArgs("SUBCYCLING STEPS",mppf->Nsubsteps);

  if(mppf->Nsubsteps){
    mppf->Ud    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
    mppf->Ue    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
    mppf->resU  = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
    mppf->rhsUd = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));

    if(mppf->elementType==HEXAHEDRA)
      mppf->cUd = (dfloat *) calloc(mppf->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
    else 
      mppf->cUd = mppf->U;
  }



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
  mppf->eta =  mppf->hmin; // Change this later, currently no inside !!!!!

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
  
#if 0

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
  // if (mesh->rank==0) printf("Nsteps = %d NerrStep= %d dt = %.8e\n", mppf->NtimeSteps,ins->outputStep, mppf->dt);

  mppf->outputForceStep = 0;
  options.getArgs("TSTEPS FOR FORCE OUTPUT", mppf->outputForceStep);

  options.getArgs("MIXING ENERGY DENSITY", mppf->chL);
  options.getArgs("MOBILITY", mppf->chM);



  // Some derived parameters
  mppf->idt  = 1.0/mppf->dt;
  mppf->eta2 = mppf->eta*mppf->eta; 
  // Define coefficients of Helmholtz solves in Chan-Hilliard equation
  mppf->chS = 1.5 * mppf->eta2 *sqrt(4.0*mppf->g0/ (mppf->chM*mppf->chL*mppf->dt));   
  mppf->chA = -mppf->chS/(2.0*mppf->eta2) * (1.0 + sqrt(1 - 4.0*mppf->g0*mppf->eta2*mppf->eta2/(mppf->chM*mppf->chL*mppf->dt*mppf->chS*mppf->chS)));   
  
  // Helmholtz solve lambda's i.e. -laplace*psi + [alpha+ S/eta^2]*psi = -Q 
  mppf->lambdaPsi = mppf->chA + mppf->chS/mppf->eta2;
  // Helmholtz solve lambda's i.e. -laplace*phi +[-alpha]*phi = -psi 
  mppf->lambdaPhi = -mppf->chA;









  // mppf->lambdaVel = mppf->g0 / (mppf->dt * mppf->nu);
 
if(mesh->rank==0){
    printf("=============WRITING INPUT PARAMETERS===================\n");

    printf("INTERFACE LENGTH\t:\t%.2e\n", mppf->eta);
    printf("INTERFACE THICKNESS\t:\t%.2e\n", mppf->hmin);
    printf("MINUM TIME STEP SIZE\t:\t%.2e\n", mppf->dt);
    printf("# SUBSTEPS\t\t:\t%d\n", mppf->Nsubsteps);
 printf("============================================================\n");
    printf("VISCOSITY PHASE 1\t:\t%.2e\n", mppf->mu1);
    printf("VISCOSITY PHASE 2\t:\t%.2e\n", mppf->mu2);
    printf("DENSITY PHASE 1\t\t:\t%.2e\n", mppf->rho1);
    printf("DENSITY PHASE 2\t\t:\t%.2e\n", mppf->rho2);
    printf("MIXING ENERGY\t\t:\t%.2e\n", mppf->chL);
    printf("MOBILITY\t\t:\t%.2e\n", mppf->chM);
    printf("CH HELMHOLTZ LAMBDA 1\t:\t%.2e\n", mppf->lambdaPsi);
    printf("CH HELMHOLTZ LAMBDA 2\t:\t%.2e\n", mppf->lambdaPhi);
 printf("============================================================\n");
    printf("# TIME STEPS\t\t:\t%d\n", mppf->NtimeSteps);
    printf("# OUTPUT STEP\t\t:\t%d\n", mppf->outputStep);
    printf("# FORCE STEP\t\t:\t%d\n", mppf->outputForceStep);
  }










return mppf;
}