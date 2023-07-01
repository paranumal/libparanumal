/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "wave.hpp"
#include "timer.hpp"
#include "ellipticPrecon.hpp"

void wave_t::Solve(deviceMemory<dfloat> &o_rDL,
                   deviceMemory<dfloat> &o_rPL,
                   deviceMemory<dfloat> &o_rFL){

  std::cout << "iostep = " << iostep << std::endl;

  linAlgMatrix_t<dfloat> filtD(1,Nstages);
  deviceMemory<dfloat> o_filtD = platform.malloc<dfloat>(Nstages);

  deviceMemory<dfloat> o_SurfaceForceL  = platform.malloc<dfloat>(mesh.Np*mesh.Nelements);
  
  dlong cumIter = 0;
  
  for(int tstep=0;tstep<Nsteps;++tstep){ // do adaptive later
    int iter = 0;
    
    dfloat t = tstep*dt;
    
    timePoint_t starts = GlobalPlatformTime(platform);
    
    // PhatL = PL, DhatL = DL, DrhsL = lambda*WJ*(scD*DL + scP*PL)
    dfloat scD = 1. + invGamma*alpha(2,1);
    dfloat scP = scD*invGamma*invDt;
    waveStepInitializeKernel(mesh.Nelements, scD, scP, lambdaSolve,
                             o_WJ, o_MM, o_rDL, o_rPL, o_DhatL, o_PhatL, o_DrhsL);

    // LOOP OVER IMPLICIT STAGES
    for(int stage=2;stage<=Nstages;++stage){

      elliptic.lambda = lambdaSolve;

      // entering surface source
      std::cout << "Entering surface source kernel " << std::endl;
      waveSurfaceSourceKernel(mesh.Nelements,
                              mesh.o_vmapM,
                              elliptic.tau,
                              mesh.o_vgeo,
                              mesh.o_sgeo,
                              o_EToPatch,
                              t,
                              dt,
                              stage-1,
                              Nstages,
                              o_alphatilde,
                              gammatilde(1,stage)/dt,
                              o_esdirkC,
                              xsource,
                              ysource,
                              zsource,
                              fsource,
                              mesh.o_x,
                              mesh.o_y,
                              mesh.o_z,
                              mesh.o_D,
                              mesh.o_LIFT,
                              mesh.o_MM,
                              o_SurfaceForceL);
      
    
      platform.linAlg().
         axpy(mesh.Np*mesh.Nelements, (dfloat)-1., o_SurfaceForceL, (dfloat)1., o_DrhsL);
      
      // record local RHS
      if(esc)
         esc->setLocalRHS(o_DrhsL);
      
      // gather rhs to globalDofs if c0
      if(disc_c0){
        elliptic.ogsMasked.Gather(o_Drhs,   o_DrhsL,   1, ogs::Add, ogs::Trans);
        elliptic.ogsMasked.Gather(o_Dtilde, o_DtildeL, 1, ogs::Add, ogs::NoTrans);
      }else{
        o_DtildeL.copyTo(o_Dtilde);
        o_DrhsL.copyTo(o_Drhs);
      }
      
      
      // changed to tol
      int iterD =
         elliptic.Solve(linearSolver, o_Dtilde, o_Drhs, tol, maxIter, verbose, stoppingCriteria);
      cumIter += iterD;

      if(disc_c0){

        // scatter x to LocalDofs if c0
        elliptic.ogsMasked.Scatter(o_DtildeL, o_Dtilde, 1, ogs::NoTrans);
        
      }else{
        o_Dtilde.copyTo(o_DtildeL);
      }
      
      iter += iterD;

      dfloat scF = 0;

      if(settings.compareSetting("SOLVER MODE", "WAVEHOLTZ")){
        if(stage<Nstages){
          for(int j=1;j<=stage+1;++j){
            dfloat tj = t + dt*esdirkC(1,j);
            scF += alpha(stage+1,j)*cos(omega*tj);
          }
        }
      }else{
        scF = 0; // no forcing if TIMEDOMAIN
      }

      // transform DtildeL to DhatL, compute DrhsL for next stage (if appropriate)
      waveStageFinalizeKernel(mesh.Nelements,
                              dt,
                              invGammaDt,
                              invGamma,
                              gamma,
                              lambdaSolve,
                              Nstages,
                              stage,
                              scF,
                              o_alpha,
                              o_alphatilde,
                              o_gammatilde,
                              o_WJ,
                              o_MM,
                              o_FL,
                              o_DtildeL, 
                              o_DhatL,
                              o_PhatL,
                              o_DrhsL); // remember 1-index
    }

    // KERNEL 4: FINALIZE
    // P = Phat(:,1) +     dt*(Dhat(:,1:Nstages)*beta'); 
    // D = Dhat(:,1) + dt*LAP*(Phat(:,1:Nstages)*beta');
    
    // a. Phat(:,1:Nstages)*beta' => o_scratchL
    waveCombineKernel(Nall, Nstages, dt, o_beta, o_betaAlpha, o_PhatL, o_DhatL, o_scratch1L);

    // b. L*(Phat(:,1:Nstages)*beta') => o_rDL
    elliptic.lambda = 0;
    if(disc_c0){
      // gather up RHS 
      elliptic.ogsMasked.Gather(o_scratch1, o_scratch1L, 1, ogs::Add, ogs::NoTrans);
      elliptic.Operator(o_scratch1, o_scratch2);
      elliptic.ogsMasked.Scatter(o_scratch2L, o_scratch2, 1, ogs::NoTrans);
    }else{
      std::cout << "waveSolve calling  elliptic Operator" << std::endl;
      elliptic.Operator(o_scratch1L, o_scratch2L);
      std::cout << "waveSolve finished elliptic Operator" << std::endl;
    }
    elliptic.lambda = lambdaSolve;

#if 1
    // finalizing force
    std::cout << "Entering finalize source kernel " << std::endl;
    waveSurfaceSourceKernel(mesh.Nelements,
                            mesh.o_vmapM,
                            elliptic.tau,
                            mesh.o_vgeo,
                            mesh.o_sgeo,
                            o_EToPatch,
                            t,
                            dt,
                            Nstages,  // last stage+1 (beta coefficients)
                            Nstages,
                            o_alphatilde,  /* use these coefficients */
                            (dfloat)1./dt,  /* use these coefficients */
                            o_esdirkC,
                            xsource,
                            ysource,
                            zsource,
                            fsource,
                            mesh.o_x,
                            mesh.o_y,
                            mesh.o_z,
                            mesh.o_D,
                            mesh.o_LIFT,
                            mesh.o_MM,
                            o_SurfaceForceL);
    
    platform.linAlg().
       axpy(mesh.Np*mesh.Nelements, dt, o_SurfaceForceL, (dfloat)1., o_scratch2L); 
#endif
    
    // c. finalize
    // P = Phat(:,1) +     dt*(Dhat(:,1:Nstages)*beta'); 
    // D = Dhat(:,1) + dt*LAP*(Phat(:,1:Nstages)*beta'); (LAP = -MM\L)
    // THIS INCLUDE MASS MATRIX INVERSE - (IF C0 USES INVERSE OF GLOBALIZED LUMPED MASS MATRIX)

    dfloat scF = 0;
    dfloat filtP = 0;
    filtD = (dfloat)0.;

    if(settings.compareSetting("SOLVER MODE", "WAVEHOLTZ")){
      for(int i=1;i<=Nstages;++i){
        dfloat ti = t + dt*esdirkC(1,i);
        scF += beta(1,i)*cos(omega*ti);
      }
      for(int i=1;i<=Nstages;++i){
        dfloat filti = 2.*(cos(omega*(t+dt*esdirkC(1,i)))-0.25)/finalTime;
        filtP += beta(1,i)*filti;
        for(int j=1;j<=i;++j){
          filtD(1,j) += beta(1,i)*filti*alpha(i,j);
        }
      }
    }else{
      scF = 0; // no forcing if TIMEDOMAIN
    }

    o_filtD.copyFrom(filtD.data);
    
    waveStepFinalizeKernel(mesh.Nelements, dt, Nstages, scF, o_beta, filtP, o_filtD,
                           o_invWJ, o_invMM, o_scratch2L,
                           o_DhatL, o_PhatL, o_FL, o_rDL, o_rPL, o_filtPL);  
    
    timePoint_t ends = GlobalPlatformTime(platform);

    printf("====> time=%g, dt=%g, step=%d, sum(iterD)=%d, ave(iterD)=%3.2f\n", t+dt, dt, tstep, iter, iter/(double)(Nstages-1));
    
    double elapsedTime = ElapsedTime(starts, ends);
    
    if ((mesh.rank==0) && verbose){
      printf("%d, " hlongFormat ", %g, %d, %g, %g; global: N, dofs, elapsed, iterations, time per node, nodes*iterations/time %s\n",
             mesh.N,
             NglobalDofs,
             elapsedTime,
             iter,
             elapsedTime/(NglobalDofs),
             NglobalDofs*((dfloat)iter/elapsedTime),
             (char*) elliptic.settings.getSetting("PRECONDITIONER").c_str());
    }

    if(1)
    if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {
      static int slice=0;
      if(tstep==0 || ((tstep+1)%iostep) == 0){
        // copy data back to host
        o_rPL.copyTo(PL);
        o_rDL.copyTo(DL);
        
        // output field files
        std::string name;
        settings.getSetting("OUTPUT FILE NAME", name);
        char fname[BUFSIZ];
        sprintf(fname, "DP_%04d_%04d.vtu",  mesh.rank, slice);
        PlotFields(DL, PL, fname);

        ++slice;

      }
    }
  }

  std::cout << "time=" << Nsteps*dt << ", Cumulative iterations: " << cumIter << std::endl;
}
