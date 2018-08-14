#include "mppf.h"

void extbdfCoefficents(mppf_t *mppf, int order);

void mppfRun(mppf_t *mppf){

  mesh_t *mesh = mppf->mesh;
  
  occa::initTimer(mesh->device);
  occaTimerTic(mesh->device,"MPPF");

  for(int tstep=0;tstep<mppf->NtimeSteps;++tstep){
  // for(int tstep=0;tstep<100;++tstep){

    if(tstep<1){
      extbdfCoefficents(mppf,tstep+1);
      // Initialize Pressure Gradient GP for non-zero Pressure IC
      mppfPressureGradient(mppf, 0.0, mppf->o_P, mppf->o_GP);      
    } 
    else if(tstep<2 && mppf->temporalOrder>=2) 
      extbdfCoefficents(mppf,tstep+1);
    else if(tstep<3 && mppf->temporalOrder>=3) 
      extbdfCoefficents(mppf,tstep+1);
  
    dfloat time = mppf->startTime + tstep*mppf->dt;
    
    dfloat time_new = time + mppf->dt; 

    // // Interface Solver
    // // Compute Rhs
    mppfCahnHilliardRhs(mppf, time_new);
    // Solve for Psi  and Phi
    mppfCahnHilliardSolve(mppf, time_new);
    // Update Phi and Compute Gradient Phi
    mppfCahnHilliardUpdate(mppf, time_new);

    // Compute Nonlinear Term N(U) 
    mppfAdvection(mppf, time, mppf->o_U, mppf->o_NU);

    // Compute intermediate velocity, save to Uhat , take divergence and compute Pr Rhs
    mppfPressureRhs(mppf, time, mppf->o_Uhat);

#if 0
    mppf->setFlowFieldKernel(mesh->Nelements,
                              time_new, 
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              mppf->fieldOffset,
                              mppf->o_rkU,
                              mppf->o_rkPhi); // pressure is not exact anymore
#endif


     mppfPressureSolve(mppf,time, mppf->o_rkP);

    
    //cycle history for update
    for (int s=mppf->Nstages;s>1;s--) {
     
      mppf->o_P.copyFrom(mppf->o_P, mppf->Ntotal*sizeof(dfloat), 
                                  (s-1)*mppf->Ntotal*sizeof(dfloat), 
                                  (s-2)*mppf->Ntotal*sizeof(dfloat));

      mppf->o_GP.copyFrom(mppf->o_GP, mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
                                  (s-1)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
                                  (s-2)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat));
    }

    // update pressure
    mppf->o_P.copyFrom(mppf->o_rkP,mppf->Ntotal*sizeof(dfloat));  
    // // update pressure gradient
    mppfPressureGradient(mppf, time_new, mppf->o_P, mppf->o_GP);   

    
    mppfVelocityRhs(mppf, time_new, mppf->o_rhsU, mppf->o_rhsV, mppf->o_rhsW);
    mppfVelocitySolve(mppf, time_new, mppf->o_rhsU, mppf->o_rhsV, mppf->o_rhsW, mppf->o_rkU);

     //cycle history 
    for (int s=mppf->Nstages;s>1;s--) {
      mppf->o_U.copyFrom(mppf->o_U, mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
                                  (s-1)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
                                  (s-2)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat));

      mppf->o_NU.copyFrom(mppf->o_NU, mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
                                  (s-1)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
                                  (s-2)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat));
    }

   // copy updated fields
   mppf->o_U.copyFrom(mppf->o_rkU,   mppf->NVfields*mppf->Ntotal*sizeof(dfloat));  
    
    
    occaTimerTic(mesh->device,"Report");

    if(mppf->outputStep){
      if(((tstep+1)%(mppf->outputStep))==0){
        if (mppf->dim==2 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d Psi - %3d Phi - %3d \n", tstep+1, mppf->NiterU, mppf->NiterV, mppf->NiterP, mppf->NiterPsi, mppf->NiterPhi);
        if (mppf->dim==3 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d Psi - %3d Phi - %3d \n", tstep+1, mppf->NiterU, mppf->NiterV, mppf->NiterW, mppf->NiterP, mppf->NiterPsi, mppf->NiterPhi);
        mppfReport(mppf, time+mppf->dt, tstep+1);

        // // Write a restart file
        // if(mppf->writeRestartFile){
        //   if(mesh->rank==0) printf("\nWriting Binary Restart File....");
        //     mppfRestartWrite(mppf, mppf->options, time+mppf->dt);
        //   if(mesh->rank==0) printf("done\n");
        // }
      }
    }

    // if (mppf->dim==2 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d Psi - %3d Phi - %3d \n", tstep+1, mppf->NiterU, mppf->NiterV, mppf->NiterP, mppf->NiterPsi, mppf->NiterPhi);
    // if (mppf->dim==3 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d Psi - %3d Phi - %3d \n", tstep+1, mppf->NiterU, mppf->NiterV, mppf->NiterW, mppf->NiterP, mppf->NiterPsi, mppf->NiterPhi); 
    
    occaTimerToc(mesh->device,"Report");
  }
  occaTimerToc(mesh->device,"MPPF");


  dfloat finalTime = mppf->NtimeSteps*mppf->dt;
  printf("\n");

  if(mppf->outputStep) mppfReport(mppf, finalTime,mppf->NtimeSteps);
  
  // if(mesh->rank==0) occa::printTimer();
}


void extbdfCoefficents(mppf_t *mppf, int order) {

  if(order==1) {
     //advection, first order in time, increment
    mppf->g0 =  1.0f; 
    dfloat extbdfB[3] = {1.0f, 0.0f, 0.0f};
    dfloat extbdfA[3] = {1.0f, 0.0f, 0.0f};
    dfloat extbdfC[3] = {1.0f, 0.0f, 0.0f};
    
    memcpy(mppf->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(mppf->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(mppf->extbdfC, extbdfC, 3*sizeof(dfloat));

    mppf->o_extbdfB.copyFrom(extbdfB);
    mppf->o_extbdfA.copyFrom(extbdfA);
    mppf->o_extbdfC.copyFrom(extbdfC);

    mppf->ExplicitOrder = 1; 

    // Define coefficients of Helmholtz solves in Chan-Hilliard equation
    mppf->chS = mppf->factorS*mppf->eta2*sqrt(4.0*mppf->g0/ (mppf->chM*mppf->chL*mppf->dt));   
    mppf->chA  = -mppf->chS/(2.0*mppf->eta2) * (1.0 - sqrt(1 - 4.0*mppf->g0*mppf->eta2*mppf->eta2/(mppf->chM*mppf->chL*mppf->dt*mppf->chS*mppf->chS)));   

    mppf->chSeta2 = mppf->chS/mppf->eta2; 
 
    // Helmholtz solve lambda's i.e. -laplace*psi + [alpha+ S/eta^2]*psi = -Q 
    mppf->lambdaPsi = mppf->chA + mppf->chSeta2;
    // Helmholtz solve lambda's i.e. -laplace*phi +[-alpha]*phi = -psi 
    mppf->lambdaPhi = -mppf->chA;
    // // 
    // mppf->lambda = mppf->g0 / (mppf->dt * mppf->nu);
    // mppf->ig0 = 1.0/mppf->g0; 

    printf("# chSeta2\t:\t%.4e\n", mppf->chSeta2);
    printf("# chAlpha\t\t:\t%.4e\n", mppf->chA);
  } else if(order==2) {
    //advection, second order in time, increment
    mppf->g0 =  1.5f;
    dfloat extbdfB[3] = {2.0f,-0.5f, 0.0f};
    dfloat extbdfA[3] = {2.0f,-1.0f, 0.0f};
    dfloat extbdfC[3] = {1.0f, 0.0f, 0.0f};

    memcpy(mppf->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(mppf->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(mppf->extbdfC, extbdfC, 3*sizeof(dfloat));

    mppf->o_extbdfB.copyFrom(extbdfB);
    mppf->o_extbdfA.copyFrom(extbdfA);
    mppf->o_extbdfC.copyFrom(extbdfC);

    mppf->ExplicitOrder=2;

    // Define coefficients of Helmholtz solves in Chan-Hilliard equation
    mppf->chS = mppf->factorS*mppf->eta2*sqrt(4.0*mppf->g0/ (mppf->chM*mppf->chL*mppf->dt));   
    mppf->chA  = -mppf->chS/(2.0*mppf->eta2) * (1.0 - sqrt(1 - 4.0*mppf->g0*mppf->eta2*mppf->eta2/(mppf->chM*mppf->chL*mppf->dt*mppf->chS*mppf->chS)));   

    mppf->chSeta2 = mppf->chS/mppf->eta2; 
 
    // Helmholtz solve lambda's i.e. -laplace*psi + [alpha+ S/eta^2]*psi = -Q 
    mppf->lambdaPsi = mppf->chA + mppf->chSeta2;
    // Helmholtz solve lambda's i.e. -laplace*phi +[-alpha]*phi = -psi 
    mppf->lambdaPhi = -mppf->chA;
    // mppf->lambda = mppf->g0 / (mppf->dt * mppf->nu);

    printf("# chSeta2\t:\t%.4e\n", mppf->chSeta2);
    printf("# chAlpha\t\t:\t%.4e\n", mppf->chA);
    // mppf->ig0 = 1.0/mppf->g0; 
  } else if(order==3) {
    //advection, third order in time, increment
    mppf->g0 =  11.f/6.f;
    dfloat extbdfB[3] = {3.0f,-1.5f, 1.0f/3.0f};
    dfloat extbdfA[3] = {3.0f,-3.0f, 1.0f};
    dfloat extbdfC[3] = {2.0f,-1.0f, 0.0f};
    
    memcpy(mppf->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(mppf->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(mppf->extbdfC, extbdfC, 3*sizeof(dfloat));

    mppf->o_extbdfB.copyFrom(extbdfB);
    mppf->o_extbdfA.copyFrom(extbdfA);
    mppf->o_extbdfC.copyFrom(extbdfC);

    mppf->ExplicitOrder=3;

    // Define coefficients of Helmholtz solves in Chan-Hilliard equation
    mppf->chS = mppf->factorS*mppf->eta2*sqrt(4.0*mppf->g0/ (mppf->chM*mppf->chL*mppf->dt));   
    mppf->chA  = -mppf->chS/(2.0*mppf->eta2) * (1.0 - sqrt(1 - 4.0*mppf->g0*mppf->eta2*mppf->eta2/(mppf->chM*mppf->chL*mppf->dt*mppf->chS*mppf->chS)));   

    mppf->chSeta2 = mppf->chS/mppf->eta2; 
 
    // Helmholtz solve lambda's i.e. -laplace*psi + [alpha+ S/eta^2]*psi = -Q 
    mppf->lambdaPsi = mppf->chA + mppf->chSeta2;
    // Helmholtz solve lambda's i.e. -laplace*phi +[-alpha]*phi = -psi 
    mppf->lambdaPhi = -mppf->chA;

   printf("# chSeta2\t:\t%.4e\n", mppf->chSeta2);
   printf("# chAlpha\t\t:\t%.4e\n", mppf->chA);
  }
}
