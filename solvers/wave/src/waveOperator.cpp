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

void wave_t::Operator(deviceMemory<dfloat> &o_QL,
                      deviceMemory<dfloat> &o_AQL){

  // THIS IMPACTS -
  // o_DhatL, o_PhatL, o_DrhsL,  o_Dtilde, o_DtildeL, o_Drhs 
  
//  std::cout << "STARTING: wave solve" << std::endl;

  linAlgMatrix_t<dfloat> filtD(1,Nstages);
  deviceMemory<dfloat> o_filtD = platform.malloc<dfloat>(Nstages);
  
  // zero velocity divergence
  platform.linAlg().set(Nall, (dfloat)0., o_DL);
  
  // zero filtered solution accumulator
  platform.linAlg().set(Nall, (dfloat)0., o_AQL);

  // copy IC
  o_QL.copyTo(o_PL);
  
  for(int tstep=0;tstep<Nsteps;++tstep){ // do adaptive later
    int iter = 0;
    
    dfloat t = tstep*dt;
    
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
        elliptic.ogsMasked.Gather(o_Dtilde, o_DtildeL, 1, ogs::Add, ogs::NoTrans);
      }else{
        o_DtildeL.copyTo(o_Dtilde);
        o_DrhsL.copyTo(o_Drhs);
      }

      int modMaxIter = maxIter; // -4;
      int iterD = elliptic.Solve(linearSolver, o_Dtilde, o_Drhs, tol, modMaxIter, verbose, stoppingCriteria);

      //add the boundary data to the masked nodes
      if(disc_c0){

        // scatter x to LocalDofs if c0
        elliptic.ogsMasked.Scatter(o_DtildeL, o_Dtilde, 1, ogs::NoTrans);
        
      }else{
        o_Dtilde.copyTo(o_DtildeL);
      }
      
      iter += iterD;

      // unforced
      dfloat scF = 0;
      
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


    waveCombineKernel(Nall, Nstages, dt, o_beta, o_betaAlpha, o_PhatL, o_DhatL, o_scratch1L);

    elliptic.lambda = 0;
    if(disc_c0){
      // gather up RHS 
      elliptic.ogsMasked.Gather(o_scratch1, o_scratch1L, 1, ogs::Add, ogs::NoTrans);
      elliptic.Operator(o_scratch1, o_scratch2);
      elliptic.ogsMasked.Scatter(o_scratch2L, o_scratch2, 1, ogs::NoTrans);
    }else{
      elliptic.Operator(o_scratch1L, o_scratch2L);
    }
    elliptic.lambda = lambdaSolve;

    dfloat filtP = 0;
    filtD = (dfloat)0.;
    for(int i=1;i<=Nstages;++i){
      dfloat filti = 2.*(cos(omega*(t+dt*esdirkC(1,i)))-0.25)/finalTime;
      filtP += beta(i)*filti;
      for(int j=1;j<=i;++j){
        filtD(1,j) += beta(i)*filti*alpha(i,j);
      }
    }
    o_filtD.copyFrom(filtD.data);
    
    // unforced
    dfloat scF = 0;
    waveStepFinalizeKernel(mesh.Nelements, dt, Nstages, scF, o_beta, filtP, o_filtD,
                           o_invWJ, o_invMM,
                           o_scratch2L, o_DhatL, o_PhatL, o_FL, o_DL, o_PL, o_AQL); 

    // printf("=*=*=*=> time=%g, dt=%g, step=%d, sum(iterD)=%d, ave(iterD)=%3.2f\n", t+dt, dt, tstep, iter, iter/(double)(Nstages-1));
  }

  // AQL = (I - OP)*QL;
  platform.linAlg().axpy(Nall, (dfloat)1., o_QL, (dfloat)-1., o_AQL);
  
  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {
    static int slice=0;
//    if(tstep>0 && (tstep%iostep) == 0){
    {
      ++slice;
      
      // copy data back to host
      o_AQL.copyTo(PL);
      o_DL.copyTo(DL);
      
      // output field files
      std::string name;
      settings.getSetting("OUTPUT FILE NAME", name);
      char fname[BUFSIZ];
      sprintf(fname, "OP_%04d_%04d.vtu",  mesh.rank, slice);
      PlotFields(DL, PL, fname);
    }
  }
}
