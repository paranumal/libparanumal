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

void wave_t::Run(){

  // set up initial conditions for D and P
  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);
  
  dfloat dt = 0.1;
  dfloat t = 0;
  int Nsteps = ceil(finalTime/dt);
  dt = finalTime/Nsteps;

  dfloat invDt = 1./dt;
  dfloat invGammaDt = 1./(gamma*dt);
  dfloat lambdaSolve = 1./(gamma*gamma*dt*dt);

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

  
  
}
