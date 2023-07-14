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

  dfloat t = startTime;
  omega = 0;
  
  if(settings.compareSetting("SOLVER MODE", "WAVEHOLTZ")){
    
    // harmonic forcing data
    settings.getSetting("OMEGA", omega);
    sigma = std::max((dfloat)36., omega*omega);
    int NouterSteps = 15; // was 30
    finalTime = NouterSteps*(2.*M_PI/omega);
    
    // should try using more accurate pressure accumulator
    Nsteps = NouterSteps*100; // THIS PART IS CRITICAL. AT OMEGA=40.5 WE NEED 10

  }else{
    // round time step
    settings.getSetting("TIME STEP", dt);
    Nsteps = finalTime/dt;
  }
  
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

  std::cout << "Ndofs = " << Ndofs << ", Nall = " << Nall << std::endl;
  
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

#if 0
  // set up some monochromatic forcing (https://arxiv.org/pdf/1910.10148.pdf)
  // harmonic spatial forcing
  waveForcingKernel(Nall, t, sigma, omega, mesh.o_x, mesh.o_y, mesh.o_z, o_FL); // will use cos(omega*t)*FL
#endif
  
  if(elliptic.settings.compareSetting("STOPPING CRITERIA", "ERRORESTIMATE")){
    esc = new ellipticStoppingCriteria<dfloat>(&elliptic, NULL);
    esc->reset();
    stoppingCriteria = esc;
  }
  else
  {
    stoppingCriteria = new stoppingCriteria_t<dfloat>();
  }


  // choose which model to run
  if(settings.compareSetting("SOLVER MODE", "WAVEHOLTZ")){
    deviceMemory<dfloat> o_qL = platform.malloc<dfloat>(Nall);
    waveHoltz(o_qL);
  }else{

    timePoint_t starts = GlobalPlatformTime(platform);
    
    Solve(o_DL, o_PL, o_FL);

    timePoint_t ends = GlobalPlatformTime(platform);

    platform.device.finish();
    
    double elapsedTime = ElapsedTime(starts, ends);
    std::cout << "elapsedTime = " << std::scientific << elapsedTime << std::endl;
    
  
    deviceMemory<dfloat> o_exactDL = platform.malloc<dfloat>(Nall);
    deviceMemory<dfloat> o_exactPL = platform.malloc<dfloat>(Nall);
    
    waveInitialConditionsKernel(Nall, finalTime, mesh.o_x, mesh.o_y, mesh.o_z, o_exactDL, o_exactPL);

    o_PL.copyTo(PL);
    o_DL.copyTo(DL);       
    
    memory<dfloat> exactDL(Nall);
    memory<dfloat> exactPL(Nall);

    o_exactDL.copyTo(exactDL);
    o_exactPL.copyTo(exactPL);

    dfloat errDLMax = 0;
    dfloat errPLMax = 0;
    for(int n=0;n<Nall;++n){
      dfloat errDLn = fabs(DL[n] - exactDL[n]);
      dfloat errPLn = fabs(PL[n] - exactPL[n]);

      errDLMax = std::max(errDLMax, errDLn);
      errPLMax = std::max(errPLMax, errPLn);
    }
    std::cout << std::setprecision (6);
    std::cout << std::scientific;
    std::cout << "errDLMax = " << errDLMax << " errPLMax = " << errPLMax << std::endl;

    dfloat errPL2 = 0, errDL2 = 0;
    dfloat normExactPL2 = 0, normExactDL2 = 0;
    // change for other elements
    if(mesh.elementType==Mesh::TRIANGLES ||
       mesh.elementType==Mesh::TETRAHEDRA){
      for(dlong e=0;e<mesh.Nelements;++e){
        for(int n=0;n<mesh.Np;++n){
          dlong idn = e*mesh.Np + n;
          dfloat errPL2n = 0, normExactPL2n = 0;
          dfloat errDL2n = 0, normExactDL2n = 0;
          for(int m=0;m<mesh.Np;++m){
            dlong idm = e*mesh.Np + m;
            dfloat MMnm = mesh.MM[n*mesh.Np+m];
            errPL2n += MMnm*(PL[idm]-exactPL[idm]);
            errDL2n += MMnm*(DL[idm]-exactDL[idm]);
            normExactPL2n += MMnm*exactPL[idm];
            normExactDL2n += MMnm*exactDL[idm];
          }
          errPL2 += WJ[e]*(PL[idn]-exactPL[idn])*errPL2n;
          errDL2 += WJ[e]*(DL[idn]-exactDL[idn])*errDL2n;
          normExactPL2 += WJ[e]*(exactPL[idn])*normExactPL2n;
          normExactDL2 += WJ[e]*(exactDL[idn])*normExactDL2n;
        }
      }
    }

    if(mesh.elementType==Mesh::QUADRILATERALS ||
       mesh.elementType==Mesh::HEXAHEDRA){
      for(dlong e=0;e<mesh.Nelements;++e){
        for(int n=0;n<mesh.Np;++n){
          dlong idn = e*mesh.Np + n;
          dfloat WJn = WJ[idn];
          errPL2 += WJn*pow(PL[idn]-exactPL[idn],2);
          errDL2 += WJn*pow(DL[idn]-exactDL[idn],2);
          
          normExactPL2 += WJn*pow(exactPL[idn],2);
          normExactDL2 += WJn*pow(exactDL[idn],2);
        }
      }
    }
    
    dfloat relErrDL2 = sqrt(errDL2/normExactDL2);
    dfloat relErrPL2 = sqrt(errPL2/normExactPL2);
    
    std::cout << mesh.N << ", " << mesh.Nelements << ", " << dt << ", " << relErrDL2 << ", " <<  relErrPL2 << ", " << elapsedTime << "; %% N, Nelements, dt, relErrDL2, relErrPL2, elapsedTime" << std::endl;
    
    if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {
      // copy data back to host
      // output field files
      std::string name;
      settings.getSetting("OUTPUT FILE NAME", name);
      char fname[BUFSIZ];
      sprintf(fname, "SOLN_DP_%04d.vtu",  mesh.rank);
      PlotFields(DL, PL, fname);
    }

    {
      std::string name;
      settings.getSetting("OUTPUT FILE NAME", name);
      char fname[BUFSIZ];
      // write to binary file
      sprintf(fname, "SOLN_DP_%04d.bin",  mesh.rank);
      FILE *fp = fopen(fname, "w");
      fprintf(fp, "# libParanumal binary format: blocks of Np*Nel doubles for (x,y,DL,PL,..). First line is: dim,Np,Nel,Nfields,sizeof(dfloat)\n");
      fprintf(fp, "%d %d %d %d %d (x,y,PL,DL)\n", mesh.dim, mesh.Np, mesh.Nelements, 4, sizeof(dfloat));
      dlong NpNel = mesh.Np*mesh.Nelements;
      fwrite(mesh.x.ptr(), NpNel, sizeof(dfloat), fp);
      fwrite(mesh.y.ptr(), NpNel, sizeof(dfloat), fp);
      fwrite(DL.ptr(), NpNel, sizeof(dfloat), fp);
      fwrite(PL.ptr(), NpNel, sizeof(dfloat), fp);
      fclose(fp);
    }
  }
  
}