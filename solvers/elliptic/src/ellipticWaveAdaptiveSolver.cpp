

#if USE_ADAPTIVE==1
  int tstep = 0;
  int allStep = 0;
  
  // START ADAPTIVE TIME-STEPPER
  dfloat dtMIN = 1E-9; //minumum allowed timestep
  dfloat ATOL = 1E-6;  //absolute error tolerance
  dfloat RTOL = 1E-6;  //relative error tolerance
  dfloat safe = 0.8;   //safety factor
    
  //error control parameters
  dfloat errBeta = 0.05;
  dfloat factor1 = 0.2;
  dfloat factor2 = 10.0;

  // ???
  dfloat exp1 = 0.2 - 0.75*errBeta;
  dfloat invfactor1 = 1.0/factor1;
  dfloat invfactor2 = 1.0/factor2;
//  dfloat facold = 1E-4;
    dfloat facold = 1E-4;
  dfloat sqrtinvNall = 1.0/sqrt(2.*(dfloat)Nall);

  dfloat errEstimate = 0;
    
  t = 0;

  while(t<finalTime){

    LIBP_ABORT("Time step became too small at time step = " << tstep,
               dt<dtMIN);
    LIBP_ABORT("Solution became unstable at time step = " << tstep,
               std::isnan(dt));
    
    //check for final timestep
    if (t+dt > finalTime){
      dt = finalTime-t;
    }
    
    // STEP
    {
      invDt = 1./dt;
      invGammaDt = 1./(gam*dt);
      lambdaSolve = 1./(gam*gam*dt*dt);
        
      int iter = 0;
        
      dfloat scD = 1. + invGamma*alpha(2,1);
      dfloat scP = scD*invGamma*invDt;
      waveStepInitializeKernel(mesh.Nelements, scD, scP, lambdaSolve,
                               o_WJ, o_MM, o_DL, o_PL, o_DhatL, o_PhatL, o_DrhsL);
        
      // LOOP OVER IMPLICIT STAGES
      for(int stage=2;stage<=Nstages;++stage){

        lambda = lambdaSolve;
        
        // solve for D
        if(esc)
           esc->setLocalRHS(o_DrhsL);

        int iterD = Solve(DLinearSolver, o_DtildeL, o_DrhsL, tol, maxIter, verbose, stoppingCriteria);
          
        iter += iterD;
          
        // transform DtildeL to DhatL, compute DrhsL for next stage (if appropriate)
        waveStageFinalizeKernel(mesh.Nelements,
                                dt,
                                invGammaDt,
                                invGamma,
                                gam,
                                lambdaSolve,
                                Nstages,
                                stage,
                                o_alpha,
                                o_alphatilde,
                                o_gammatilde,
                                o_WJ,
                                o_MM,
                                o_DtildeL, 
                                o_DhatL,
                                o_PhatL,
                                o_DrhsL); // remember 1-index
      }
      //      printf("====> step=%d, sum(iterD)=%d, ave(iterD)=%3.2f\n", tstep, iter, iter/(double)(Nstages-1));
      printf("====> time=%g, dt=%g, step=%d, sum(iterD)=%d, ave(iterD)=%3.2f\n", t, dt, tstep, iter, iter/(double)(Nstages-1));
        
      // finalize: 
//      waveCombineKernel(Nall, Nstages, o_beta, o_PhatL, o_scratch1);
      waveCombineKernel(Nall, Nstages, dt, o_beta, o_betaAlpha, o_PhatL, o_DhatL, o_scratch1);
      lambda = 0;
      Operator(o_scratch1, o_scratch2);
      lambda = lambdaSolve;
      waveStepFinalizeKernel(mesh.Nelements, dt, Nstages, o_beta, o_invWJ, o_invMM, o_scratch2, o_DhatL, o_PhatL, o_newDL, o_newPL);

      platform.device.finish();
      
#if 1
      // finalize with embedded  coefficients
      // TW NEEDS TO BE UPDATED
      waveCombineKernel(Nall, Nstages, o_betahat, o_betahatAlpha, o_PhatL, p_DhatL, o_scratch1);
      lambda = 0;
      Operator(o_scratch1, o_scratch2);
      lambda = lambdaSolve;
      waveStepFinalizeKernel(mesh.Nelements, dt, Nstages, o_betahat, o_invWJ, o_invMM, o_scratch2, o_DhatL, o_PhatL, o_embDL, o_embPL);
#endif
      
      // estimate error
      waveErrorEstimateKernel(Nall, ATOL, RTOL, o_DL, o_PL, o_newDL, o_newPL, o_embDL, o_embPL, o_errL);

#if 1
      errEstimate = platform.linAlg().sum(Nall, o_errL, mesh.comm);
      errEstimate = sqrt(errEstimate)*sqrtinvNall;
#endif
      printf("ESDIRK: errEStimate=%f\n", errEstimate);
    } // END STEP

      // build controller
    dfloat fac1 = pow(errEstimate,exp1);
    dfloat fac = fac1/pow(facold,errBeta);
      
    fac = std::max(invfactor2, std::min(invfactor1,fac/safe));
    dfloat dtnew = dt/fac;
    dtnew = dt;
    
    if (errEstimate<1.0) { //dt is accepted
      std::cout << "ESDIRK: Accepting step" << std::endl;

      // accept rkq
      o_DL.copyFrom(o_newDL);
      o_PL.copyFrom(o_newPL);

      t += dt;
      
      constexpr dfloat errMax = 1.0e-4;  // hard coded factor ?
      facold = std::max(errEstimate,errMax);
      
      tstep++;
    } else {
      dtnew = dt/(std::max(invfactor1,fac1/safe));
      printf("ESDIRK: old dt = %g, dtnew = %g\n", dt, dtnew);
    }
    if(dtnew>dt){
      printf("ESDIRK: INCREASING FROM old dt = %g, dtnew = %g\n", dt, dtnew);
    }
    dt = dtnew;
    allStep++;
  }
   
    
  // END NON-ADAPTIVE TIME-STEPPER
#endif
    
  
