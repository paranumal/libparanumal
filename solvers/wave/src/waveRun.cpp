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

void wave_t::Run(){

  // set up initial conditions for D and P
  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);
  
  dfloat dt = 0.1;
  dfloat t = 0;
  int Nsteps = ceil(finalTime/dt);
  dt = finalTime/Nsteps;

  dfloat invGamma = 1./gamma;
  dfloat invDt = 1./dt;
  dfloat invGammaDt = 1./(gamma*dt);

  lambdaSolve = 1./(gamma*gamma*dt*dt);

  int NBCTypes = 7;
  memory<int> BCType(NBCTypes);
  // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  BCType[0] = 0;
  BCType[1] = 1;
  BCType[2] = 1;
  BCType[3] = 2;
  BCType[4] = 1;
  BCType[5] = 2;
  BCType[6] = 2;


  
  elliptic.Setup(platform, mesh, ellipticSettings, lambdaSolve, NBCTypes, BCType);
  
  int maxIter = 5000;
  int verbose = settings.compareSetting("VERBOSE", "TRUE") ? 1 : 0;
  dfloat tol = (sizeof(dfloat)==sizeof(double)) ? 1.0e-10 : 1.0e-5; // TW !!!
  elliptic.settings.getSetting("ITERATIVE CONVERGENCE TOLERANCE", tol);
  
  dlong Ndofs = elliptic.Ndofs;
  dlong Nhalo = elliptic.Nhalo;

  // rebuild precon for this lambda
  elliptic.lambda = lambdaSolve;
  if     (elliptic.settings.compareSetting("PRECONDITIONER", "JACOBI"))
     elliptic.precon.Setup<JacobiPrecon>(elliptic);
  else if(elliptic.settings.compareSetting("PRECONDITIONER", "MASSMATRIX"))
     elliptic.precon.Setup<MassMatrixPrecon>(elliptic);
  else if(elliptic.settings.compareSetting("PRECONDITIONER", "PARALMOND"))
     elliptic.precon.Setup<ParAlmondPrecon>(elliptic);
  else if(elliptic.settings.compareSetting("PRECONDITIONER", "MULTIGRID"))
     elliptic.precon.Setup<MultiGridPrecon>(elliptic);
  else if(elliptic.settings.compareSetting("PRECONDITIONER", "SEMFEM"))
     elliptic.precon.Setup<SEMFEMPrecon>(elliptic);
  else if(elliptic.settings.compareSetting("PRECONDITIONER", "OAS"))
     elliptic.precon.Setup<OASPrecon>(elliptic);
  else if(elliptic.settings.compareSetting("PRECONDITIONER", "NONE"))
     elliptic.precon.Setup<IdentityPrecon>(Ndofs);


  
  // build linear solvers for P and D, but typically only use one
  if (elliptic.settings.compareSetting("LINEAR SOLVER","NBPCG")){
    linearSolver.Setup<LinearSolver::nbpcg<dfloat> >(Ndofs, Nhalo, platform, elliptic.settings, comm);
  } else if (elliptic.settings.compareSetting("LINEAR SOLVER","NBFPCG")){
    linearSolver.Setup<LinearSolver::nbfpcg<dfloat> >(Ndofs, Nhalo, platform, elliptic.settings, comm);
  } else if (elliptic.settings.compareSetting("LINEAR SOLVER","PCG")){
    linearSolver.Setup<LinearSolver::pcg<dfloat> >(Ndofs, Nhalo, platform, elliptic.settings, comm);
  } else if (elliptic.settings.compareSetting("LINEAR SOLVER","PGMRES")){
    linearSolver.Setup<LinearSolver::pgmres<dfloat> >(Ndofs, Nhalo, platform, elliptic.settings, comm);
  } else if (elliptic.settings.compareSetting("LINEAR SOLVER","PMINRES")){
    linearSolver.Setup<LinearSolver::pminres<dfloat> >(Ndofs, Nhalo, platform, elliptic.settings, comm);
  }

  if (elliptic.settings.compareSetting("INITIAL GUESS STRATEGY", "LAST")) {
    linearSolver.SetupInitialGuess<InitialGuess::Last<dfloat> >(Ndofs, platform, elliptic.settings, comm);
  } else if (elliptic.settings.compareSetting("INITIAL GUESS STRATEGY", "ZERO")) {
    linearSolver.SetupInitialGuess<InitialGuess::Zero<dfloat> >(Ndofs, platform, elliptic.settings, comm);
  } else if (elliptic.settings.compareSetting("INITIAL GUESS STRATEGY", "CLASSIC")) {
    linearSolver.SetupInitialGuess<InitialGuess::ClassicProjection<dfloat> >(Ndofs, platform, elliptic.settings, comm);
  } else if (elliptic.settings.compareSetting("INITIAL GUESS STRATEGY", "QR")) {
    linearSolver.SetupInitialGuess<InitialGuess::RollingQRProjection<dfloat> >(Ndofs, platform, elliptic.settings, comm);
  } else if (elliptic.settings.compareSetting("INITIAL GUESS STRATEGY", "EXTRAP")) {
    linearSolver.SetupInitialGuess<InitialGuess::Extrap<dfloat> >(Ndofs, platform, elliptic.settings, comm);
  }

  //setup linear solver
  disc_c0 = elliptic.settings.compareSetting("DISCRETIZATION", "CONTINUOUS") ? 1 : 0;
  
  Nall = mesh.Np*(mesh.Nelements+mesh.totalHaloPairs);

  // set initial conditions
  waveInitialConditionsKernel(Nall, t, mesh.o_x, mesh.o_y, mesh.o_z, o_DL, o_PL);

  stoppingCriteria_t<dfloat> *stoppingCriteria = NULL;
  ellipticStoppingCriteria<dfloat> *esc = NULL;
  
  if(elliptic.settings.compareSetting("STOPPING CRITERIA", "ERRORESTIMATE")){
    printf("SETTING UP ESC\n");
    esc = new ellipticStoppingCriteria<dfloat>(&elliptic, NULL);
    esc->reset();
    stoppingCriteria = esc;
  }
  else
  {
    stoppingCriteria = new stoppingCriteria_t<dfloat>();
  }

  // WILL MOVE THIS TO Solve

  for(int tstep=0;tstep<Nsteps;++tstep){ // do adaptive later
    int iter = 0;
    
    t = tstep*dt;

    timePoint_t starts = GlobalPlatformTime(platform);

    // PhatL = PL, DhatL = DL, DrhsL = lambda*WJ*(scD*DL + scP*PL)
    dfloat scD = 1. + invGamma*alpha(2,1);
    dfloat scP = scD*invGamma*invDt;
    waveStepInitializeKernel(mesh.Nelements, scD, scP, lambdaSolve,
                             o_WJ, o_MM, o_DL, o_PL, o_DhatL, o_PhatL, o_DrhsL);

    // LOOP OVER IMPLICIT STAGES
    for(int stage=2;stage<=Nstages;++stage){

      elliptic.lambda = lambdaSolve;

      // record local RHS
      if(esc)
         esc->setLocalRHS(o_DrhsL);

      // gather rhs to globalDofs if c0
      if(disc_c0){
        elliptic.ogsMasked.Gather(o_Drhs,   o_DrhsL,   1, ogs::Add, ogs::Trans);
        
        // TW not sure about this (working assumption is that DtildeL is already implicitly continuous)
        elliptic.ogsMasked.Gather(o_Dtilde, o_DtildeL, 1, ogs::Add, ogs::NoTrans);
      }else{
        o_DtildeL.copyTo(o_Dtilde);
        o_DrhsL.copyTo(o_Drhs);
      }
      
      int iterD = elliptic.Solve(linearSolver, o_Dtilde, o_Drhs, tol, maxIter, verbose, stoppingCriteria);

      //add the boundary data to the masked nodes
      if(disc_c0){

        // scatter x to LocalDofs if c0
        elliptic.ogsMasked.Scatter(o_DtildeL, o_Dtilde, 1, ogs::NoTrans);
        
      }else{
        o_Dtilde.copyTo(o_DtildeL);
      }
      
      iter += iterD;
      
      // transform DtildeL to DhatL, compute DrhsL for next stage (if appropriate)
      waveStageFinalizeKernel(mesh.Nelements,
                              dt,
                              invGammaDt,
                              invGamma,
                              gamma,
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

    // KERNEL 4: FINALIZE
    // P = Phat(:,1) +     dt*(Dhat(:,1:Nstages)*beta'); 
    // D = Dhat(:,1) + dt*LAP*(Phat(:,1:Nstages)*beta');
    
    // a. Phat(:,1:Nstages)*beta' => o_scratchL
    waveCombineKernel(Nall, Nstages, dt, o_beta, o_betaAlpha, o_PhatL, o_DhatL, o_scratch1L);

    // b. L*(Phat(:,1:Nstages)*beta') => o_DL
    elliptic.lambda = 0;
    if(disc_c0){
      // gather up RHS 
      elliptic.ogsMasked.Gather(o_scratch1, o_scratch1L, 1, ogs::Add, ogs::NoTrans);
      Operator(o_scratch1, o_scratch2);
      elliptic.ogsMasked.Scatter(o_scratch2L, o_scratch2, 1, ogs::NoTrans);
    }else{
      Operator(o_scratch1L, o_scratch2L);
    }
    elliptic.lambda = lambdaSolve;
    
    // c. finalize
    // P = Phat(:,1) +     dt*(Dhat(:,1:Nstages)*beta'); 
    // D = Dhat(:,1) + dt*LAP*(Phat(:,1:Nstages)*beta'); (LAP = -MM\L)
    // THIS INCLUDE MASS MATRIX INVERSE - (IF C0 USES INVERSE OF GLOBALIZED LUMPED MASS MATRIX)
    waveStepFinalizeKernel(mesh.Nelements, dt, Nstages, o_beta, o_invWJ, o_invMM, o_scratch2L, o_DhatL, o_PhatL, o_DL, o_PL); 

    timePoint_t ends = GlobalPlatformTime(platform);

    printf("====> time=%g, dt=%g, step=%d, sum(iterD)=%d, ave(iterD)=%3.2f\n", t, dt, tstep, iter, iter/(double)(Nstages-1));
    
    double elapsedTime = ElapsedTime(starts, ends);
    
    if ((mesh.rank==0) && verbose){
      printf("%d, " hlongFormat ", %g, %d, %g, %g; global: N, dofs, elapsed, iterations, time per node, nodes*iterations/time %s\n",
             mesh.N,
             NglobalDofs,
             elapsedTime,
             iter,
             elapsedTime/(NglobalDofs),
             NglobalDofs*((dfloat)iter/elapsedTime),
             (char*) settings.getSetting("PRECONDITIONER").c_str());
    }

    int iostep = 1;
    if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {
      static int slice=0;
      if((tstep%iostep) == 0){
        ++slice;
        
        // copy data back to host
        o_PL.copyTo(PL);
        o_DL.copyTo(DL);

#if 0
        {
          std::cout << "DL === ";
          for(int n=0;n<100;++n){
            std::cout << DL[n] << ",";
          }
          std::cout << std::endl;
        }
        
        {
          dfloat val = platform.linAlg().norm2(Nall, o_DL, mesh.comm);
          dfloat minval = platform.linAlg().min(Nall, o_DL, mesh.comm);
          dfloat maxval = platform.linAlg().max(Nall, o_DL, mesh.comm);
          std::cout << "|DL|=" << val << "range=[" << minval << "," << maxval << "]" << std::endl;
        }

        {
          dfloat val = platform.linAlg().norm2(Nall, o_PL, mesh.comm);
          std::cout << "|PL|=" << val << std::endl;
        }
#endif
        
      // output field files
        std::string name;
        settings.getSetting("OUTPUT FILE NAME", name);
        char fname[BUFSIZ];
        sprintf(fname, "DP_%04d_%04d.vtu",  mesh.rank, slice);
        PlotFields(DL, PL, fname);
      }
    }
  }

#if 0
    // output norm of final solution
    {
      //compute q.M*q
      dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
      deviceMemory<dfloat> o_MDL = platform.reserve<dfloat>(Nentries);
      mesh.MassMatrixApply(o_DL, o_MDL);
      
      dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_DL, o_MDL, mesh.comm));
      
      if(mesh.rank==0)
         printf("Solution norm = %17.15lg\n", norm2);
    }
#endif

  
  
  
}
