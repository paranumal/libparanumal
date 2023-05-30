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
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  // harmonic forcing data
  settings.getSetting("OMEGA", omega);
  sigma = std::max(36., omega*omega);
  int NouterSteps = 10;
  finalTime = NouterSteps*(2.*M_PI/omega);
  dfloat t = startTime;
  
  // round time step
  settings.getSetting("TIME STEP", dt);
  
  Nsteps = std::max(NouterSteps*10., ceil(finalTime/dt));
  dt = finalTime/Nsteps;
  std::cout << "dt=" << dt << std::endl;
  
  iostep = 1;
  settings.getSetting("OUTPUT STEP", iostep);
  
  invGamma = 1./gamma;
  invDt = 1./dt;
  invGammaDt = 1./(gamma*dt);

  lambdaSolve = 1./(gamma*gamma*dt*dt);
  
  maxIter = 5000;
  verbose = settings.compareSetting("VERBOSE", "TRUE") ? 1 : 0;
  tol = (sizeof(dfloat)==sizeof(double)) ? 1.0e-10 : 1.0e-5; // TW !!!

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
  
  // build linear solvers for elliptic
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

  // build initial guess strategy
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

  // set initial conditions
  waveInitialConditionsKernel(Nall, t, mesh.o_x, mesh.o_y, mesh.o_z, o_DL, o_PL);

  // set up some monochromatic forcing (https://arxiv.org/pdf/1910.10148.pdf)
  // harmonic spatial forcing
  waveForcingKernel(Nall, t, sigma, omega, mesh.o_x, mesh.o_y, mesh.o_z, o_FL); // will use cos(omega*t)*FL
  
  if(elliptic.settings.compareSetting("STOPPING CRITERIA", "ERRORESTIMATE")){
    esc = new ellipticStoppingCriteria<dfloat>(&elliptic, NULL);
    esc->reset();
    stoppingCriteria = esc;
  }
  else
  {
    stoppingCriteria = new stoppingCriteria_t<dfloat>();
  }

  // integrate between startTime and endTime
  Solve(o_DL, o_PL, o_FL);

  // try WaveHoltz (not sure if this is quite right yet)
  waveHoltz(o_DL, o_PL, o_FL);
  
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
