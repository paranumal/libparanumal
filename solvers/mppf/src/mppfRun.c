#include "mppf.h"

void extbdfCoefficents(mppf_t *mppf, int order);

void mppfRun(mppf_t *mppf){

  mesh_t *mesh = mppf->mesh;
  
  occa::initTimer(mesh->device);
  occaTimerTic(mesh->device,"MPPF");

  // char fname[BUFSIZ];
  // string outName;
  // mppf->options.getArgs("OUTPUT FILE NAME", outName);

  // mppf->o_Phi.copyTo(mppf->Phi);
  // mppf->o_U.copyTo(mppf->U);
  // mppf->o_P.copyTo(mppf->P);
  // sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(), mesh->rank, mppf->frame++);

  // mppfPlotVTU(mppf, fname);

  // Write Initial Data
  if(mppf->outputStep) mppfReport(mppf, mppf->startTime, 0);

  for(int tstep=0;tstep<mppf->NtimeSteps;++tstep){
  // for(int tstep=0;tstep<100;++tstep){

    if(tstep<1) 
      extbdfCoefficents(mppf,tstep+1);
    else if(tstep<2 && mppf->temporalOrder>=2) 
      extbdfCoefficents(mppf,tstep+1);
    else if(tstep<3 && mppf->temporalOrder>=3) 
      extbdfCoefficents(mppf,tstep+1);
  
    dfloat time = mppf->startTime + tstep*mppf->dt;
    
    dfloat time_new = time + mppf->dt; 

    // printf("Calling CF Rhs Function\n");
    mppfPhaseFieldRhs(mppf, time_new);

    mppfCahnHilliardSolve(mppf, time, mppf->o_rkPhi);


    mppf->setFlowFieldKernel(mesh->Nelements,
                              time_new, 
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              mppf->fieldOffset,
                              mppf->o_rkU,
                              mppf->o_rkP);

     //cycle history
    for (int s=mppf->Nstages;s>1;s--) {
      mppf->o_U.copyFrom(mppf->o_U, mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
                                  (s-1)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
                                  (s-2)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat));
      mppf->o_P.copyFrom(mppf->o_P, mppf->Ntotal*sizeof(dfloat), 
                                  (s-1)*mppf->Ntotal*sizeof(dfloat), 
                                  (s-2)*mppf->Ntotal*sizeof(dfloat));

      mppf->o_Phi.copyFrom(mppf->o_Phi, mppf->Ntotal*sizeof(dfloat), 
                                  (s-1)*mppf->Ntotal*sizeof(dfloat), 
                                  (s-2)*mppf->Ntotal*sizeof(dfloat));
    }


     //copy updated fields
    mppf->o_Phi.copyFrom(mppf->o_rkPhi, mppf->Ntotal*sizeof(dfloat));
    mppf->o_U.copyFrom(  mppf->o_rkU,   mppf->NVfields*mppf->Ntotal*sizeof(dfloat));  
    mppf->o_P.copyFrom(  mppf->o_rkP,                  mppf->Ntotal*sizeof(dfloat));  



     //cycle rhs history
    for (int s=mppf->Nstages;s>1;s--) {
      // mppf->o_NU.copyFrom(mppf->o_NU, mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
      //                             (s-1)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
      //                             (s-2)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat));
      // mppf->o_GP.copyFrom(mppf->o_GP, mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
      //                             (s-1)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
      //                             (s-2)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat));

       mppf->o_NPhi.copyFrom(mppf->o_NPhi, mppf->Ntotal*sizeof(dfloat), 
                                  (s-1)*mppf->Ntotal*sizeof(dfloat), 
                                  (s-2)*mppf->Ntotal*sizeof(dfloat));
    }


    

    // if(mppf->Nsubsteps) {
    //   mppfSubCycle(mppf, time, mppf->Nstages, mppf->o_U, mppf->o_NU);
    // } else {
    //   mppfAdvection(mppf, time, mppf->o_U, mppf->o_NU);
    // } 
    // mppfGradient (mppf, time, mppf->o_P, mppf->o_GP);

    // mppfVelocityRhs  (mppf, time+mppf->dt, mppf->Nstages, mppf->o_rhsU, mppf->o_rhsV, mppf->o_rhsW);
    // mppfVelocitySolve(mppf, time+mppf->dt, mppf->Nstages, mppf->o_rhsU, mppf->o_rhsV, mppf->o_rhsW, mppf->o_rkU);

    // mppfPressureRhs  (mppf, time+mppf->dt, mppf->Nstages);
    // mppfPressureSolve(mppf, time+mppf->dt, mppf->Nstages); 

    // mppfPressureUpdate(mppf, time+mppf->dt, mppf->Nstages, mppf->o_rkP);
    // mppfGradient(mppf, time+mppf->dt, mppf->o_rkP, mppf->o_rkGP);

    // //cycle history
    // for (int s=mppf->Nstages;s>1;s--) {
    //   mppf->o_U.copyFrom(mppf->o_U, mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
    //                               (s-1)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
    //                               (s-2)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat));
    //   mppf->o_P.copyFrom(mppf->o_P, mppf->Ntotal*sizeof(dfloat), 
    //                               (s-1)*mppf->Ntotal*sizeof(dfloat), 
    //                               (s-2)*mppf->Ntotal*sizeof(dfloat));
    // }

    // //copy updated pressure
    // mppf->o_P.copyFrom(mppf->o_rkP, mppf->Ntotal*sizeof(dfloat)); 

    // //update velocity
    // mppfVelocityUpdate(mppf, time+mppf->dt, mppf->Nstages, mppf->o_rkGP, mppf->o_rkU);

    // //copy updated pressure
    // mppf->o_U.copyFrom(mppf->o_rkU, mppf->NVfields*mppf->Ntotal*sizeof(dfloat)); 

    // //cycle rhs history
    // for (int s=mppf->Nstages;s>1;s--) {
    //   mppf->o_NU.copyFrom(mppf->o_NU, mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
    //                               (s-1)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
    //                               (s-2)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat));
    //   mppf->o_GP.copyFrom(mppf->o_GP, mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
    //                               (s-1)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat), 
    //                               (s-2)*mppf->Ntotal*mppf->NVfields*sizeof(dfloat));
    // }

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

    // mppf->lambda = mppf->g0 / (mppf->dt * mppf->nu);
    // mppf->ig0 = 1.0/mppf->g0; 
  }
}
