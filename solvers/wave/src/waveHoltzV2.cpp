
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

void wave_t::waveHoltzV2(deviceMemory<dfloat> &o_qL){
  
#if 1

  const dlong Nhalo = elliptic.Nhalo;
  
  precon_t waveHoltzPrecon;
  waveHoltzPrecon.Setup<IdentityPrecon>(Nall);

  linearSolver_t<dfloat> waveHoltzLinearSolver;

  // configure linear solver to use pgmres
  waveHoltzLinearSolver.Setup<LinearSolver::pgmres<dfloat> >
     (Nall, Nhalo, platform, elliptic.settings, comm);

  // configure linear solver to use zero initial guess
  waveHoltzLinearSolver.SetupInitialGuess<InitialGuess::Last<dfloat> >
     (Nall, platform, elliptic.settings, comm);
  
  // compute RHS (this will work for IPDG for now)
  deviceMemory<dfloat> o_bL = platform.malloc<dfloat>(Nall);
  
  platform.linAlg().set(Nall, (dfloat)0, o_qL);
  platform.linAlg().set(Nall, (dfloat)0, o_DL);
  platform.linAlg().set(Nall, (dfloat)0, o_PL);
  platform.linAlg().set(Nall, (dfloat)0, o_filtPL);

  // solve homogeneous initial conditions with forcing to obtain RHS
//  Nsteps *=20;
  SolveV2(o_DL, o_PL, o_FL);
//  Nsteps /=20;
  std::cout << "**************** DONE INITIAL FORCING PHASE *******************" << std::endl;
  
  // this becomes the RHS 
  o_filtPL.copyTo(o_bL);
  
  stoppingCriteria_t<dfloat> *waveHoltzStoppingCriteria = new stoppingCriteria_t<dfloat>();
  int iterD =
     waveHoltzLinearSolver.Solve(*this,
                                 waveHoltzPrecon,
                                 o_qL,
                                 o_bL,
                                 tol,
                                 maxIter,
                                 verbose,
                                 waveHoltzStoppingCriteria);
  
  std::cout << " WOWSA " << std::endl;
  
//  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {
  {
    // copy data back to host
    o_qL.copyTo(PL);
    o_qL.copyTo(DL);
    
    // output field files
    std::string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "SOLN_%04d.vtu",  mesh.rank);
    PlotFields(DL, PL, fname);
  }
  
#endif

#if 0  
  // start with zeroed P, L, and filtered pressure
  platform.linAlg().set(Nall, (dfloat)0, o_DL);
  platform.linAlg().set(Nall, (dfloat)0, o_PL);

  // arbitrarily try 100 outer iterations
  for(int it=1;it<=100;++it){

    std::cout << "STARTING OUTER ITERATION: " << it << std::endl;
    
    platform.linAlg().set(Nall, (dfloat)0, o_DL);
    platform.linAlg().set(Nall, (dfloat)0, o_filtPL);

    // solve to T = 2pi/omega
    Solve(o_DL, o_PL, o_FL);

    // PL <= filtPL
    o_filtPL.copyTo(o_PL);
    
    if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {
      static int slice=0;
      ++slice;
      
      // copy data back to host
      o_PL.copyTo(PL);
      o_DL.copyTo(DL);
      
      // output field files
      std::string name;
      settings.getSetting("OUTPUT FILE NAME", name);
      char fname[BUFSIZ];
      sprintf(fname, "ITER_DP_%04d_%04d.vtu",  mesh.rank, slice);
      PlotFields(DL, PL, fname);
    }
  }
#endif
}
