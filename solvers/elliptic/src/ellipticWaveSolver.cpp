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

#include "elliptic.hpp"
#include "timeStepper.hpp"
#include "ellipticPrecon.hpp"

#include "timer.hpp"

template <typename T>
void printMatrixLocal(linAlgMatrix_t<T> &A, const char *str){

#if 1
  std::cout << "matrix: " <<  std::string(str) << "[" << std::endl;
  for(int r=1;r<=A.rows();++r){
    for(int c=1;c<=A.cols();++c){
      //      std::cout << A(r,c) << " ";
      printf("% 5.4e ", A(r,c));
      
    }
    std::cout << std::endl;
  }
  std::cout << std::endl << "]" << std::endl;
#endif
}


void elliptic_t::WaveSolver(){

  //setup linear algebra module
  platform.linAlg().InitKernels({"set", "sum", "norm2", "min", "max"});
  
  // FOR THIS - WE ASSUME IPDG  and LAMBDA=1
  
  Nfields = 1;

  // initialize time stepper
  linAlgMatrix_t<dfloat> BTABLE;
  int Nstages = 0;
  int embedded = 0;
 
  // libp::TimeStepper::butcherTables("ESDIRK5(3)6L[2]SA", Nstages, embedded, BTABLE);
//  libp::TimeStepper::butcherTables("ESDIRK5(4)7L[2]SA2", Nstages, embedded, BTABLE);
  libp::TimeStepper::butcherTables("ESDIRK6(5)9L[2]SA", Nstages, embedded, BTABLE);
//    libp::TimeStepper::butcherTables("Kvaerno(7,4,5)-ESDIRK", Nstages, embedded, BTABLE);
//    libp::TimeStepper::butcherTables("ESDIRK3(2)4L[2]SA", Nstages, embedded, BTABLE);
  
//  libp::TimeStepper::butcherTables("TRBDF2-ESDIRK", Nstages, embedded, BTABLE);
 

  // extract alphas, betas
  linAlgMatrix_t<dfloat> alpha(Nstages, Nstages), beta(1,Nstages), betahat(1,Nstages);
  for(int n=1;n<=Nstages;++n){
    for(int m=1;m<=Nstages;++m){
      alpha(n,m) = BTABLE(n,m+1);
    }
    beta(1,n) = BTABLE(Nstages+1,n+1);

    if(embedded)
       betahat(1,n) = BTABLE(Nstages+2,n+1);

  }

  dfloat gam = alpha(2,2);
  dfloat invGamma = 1./gam;
  
  std::cout << "gamma = " << gam << std::endl;

  linAlgMatrix_t<dfloat> alphatilde(Nstages, Nstages);
  linAlgMatrix_t<dfloat> gammatilde(1,Nstages);
  linAlgMatrix_t<dfloat> betaAlpha(1,Nstages);
  linAlgMatrix_t<dfloat> betahatAlpha(1,Nstages);
  alphatilde = 0.;
  gammatilde = 0.;
  
  for(int i=0;i<=Nstages-1;++i){
//    alphatilde(i+1,1) = 1.;
    
    for(int j=1;j<=i;++j){
      alphatilde(i+1,j) += alpha(i+1,j)/gam;
    }
    for(int j=1;j<=i;++j){
      for(int k=1;k<=j;++k){
        alphatilde(i+1,k) += alpha(i+1,j)*alpha(j,k)/(gam*gam);
      }
    }
    gammatilde(1,i+1) = 1./gam;
    for(int j=1;j<=i;++j){
      gammatilde(1,i+1) += alpha(i+1,j)/(gam*gam);
    }
  }

  betaAlpha = 0.;
  betahatAlpha = 0.;
  for(int i=1;i<=Nstages;++i){
    for(int j=1;j<=i;++j){
      betaAlpha(1,j) += beta(1,i)*alpha(i,j);
      betahatAlpha(1,j) += betahat(1,i)*alpha(i,j);
    }
  }

#if 0
  printMatrixLocal(BTABLE, "BUTCHER TABLE");
  printMatrixLocal(alpha, "ALPHA BLOCK");
  printMatrixLocal(beta, "BETA BLOCK");
  printMatrixLocal(gammatilde, "GAMMA TILDE BLOCK");
  printMatrixLocal(alphatilde, "ALPHA TILDE BLOCK");
  printMatrixLocal(betaAlpha,  "BETAxALPHA  BLOCK");
#endif
  
  deviceMemory<dfloat> o_alphatilde = platform.malloc<dfloat>(Nstages*Nstages, alphatilde.data);
  deviceMemory<dfloat> o_gammatilde = platform.malloc<dfloat>(Nstages, gammatilde.data);
  deviceMemory<dfloat> o_betaAlpha  = platform.malloc<dfloat>(Nstages, betaAlpha.data);
  deviceMemory<dfloat> o_betahatAlpha  = platform.malloc<dfloat>(Nstages, betahatAlpha.data);
  
  deviceMemory<dfloat> o_alpha = platform.malloc<dfloat>(Nstages*Nstages, alpha.data);
  deviceMemory<dfloat> o_beta  = platform.malloc<dfloat>(Nstages, beta.data);
  deviceMemory<dfloat> o_betahat  = platform.malloc<dfloat>(Nstages, betahat.data);
  
  // set up initial conditions for D and P
  dfloat finalTime = 10;
  //  dfloat dt = 0.1;
  dfloat dt = 0.1;
  dfloat t = 0;
  int Nsteps = ceil(finalTime/dt);
  dt = finalTime/Nsteps;

  dfloat invDt = 1./dt;
  dfloat invGammaDt = 1./(gam*dt);
  dfloat lambdaSolve = 1./(gam*gam*dt*dt);
  
  int maxIter = 5000;
  int verbose = settings.compareSetting("VERBOSE", "TRUE") ? 1 : 0;
  dfloat tol = (sizeof(dfloat)==sizeof(double)) ? 1.0e-10 : 1.0e-5; // TW !!!
  settings.getSetting("ITERATIVE CONVERGENCE TOLERANCE", tol);
  
  // rebuild precon for this lambda
  lambda = lambdaSolve;
  if     (settings.compareSetting("PRECONDITIONER", "JACOBI"))
     precon.Setup<JacobiPrecon>(*this);
  else if(settings.compareSetting("PRECONDITIONER", "MASSMATRIX"))
     precon.Setup<MassMatrixPrecon>(*this);
  else if(settings.compareSetting("PRECONDITIONER", "PARALMOND"))
     precon.Setup<ParAlmondPrecon>(*this);
  else if(settings.compareSetting("PRECONDITIONER", "MULTIGRID"))
    precon.Setup<MultiGridPrecon>(*this);
  else if(settings.compareSetting("PRECONDITIONER", "SEMFEM"))
     precon.Setup<SEMFEMPrecon>(*this);
  else if(settings.compareSetting("PRECONDITIONER", "OAS"))
     precon.Setup<OASPrecon>(*this);
  else if(settings.compareSetting("PRECONDITIONER", "NONE"))
    precon.Setup<IdentityPrecon>(Ndofs);
  
  
  // build linear solvers for P and D, but typically only use one
  linearSolver_t<dfloat> PLinearSolver, DLinearSolver;
  if (settings.compareSetting("LINEAR SOLVER","NBPCG")){
    DLinearSolver.Setup<LinearSolver::nbpcg<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","NBFPCG")){
    DLinearSolver.Setup<LinearSolver::nbfpcg<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
//  } else if (settings.compareSetting("LINEAR SOLVER","RPCG")){
//    DLinearSolver.Setup<LinearSolver::rpcg<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","PCG")){
    DLinearSolver.Setup<LinearSolver::pcg<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","PGMRES")){
    DLinearSolver.Setup<LinearSolver::pgmres<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  } else if (settings.compareSetting("LINEAR SOLVER","PMINRES")){
    DLinearSolver.Setup<LinearSolver::pminres<dfloat> >(Ndofs, Nhalo, platform, settings, comm);
  }

  if (settings.compareSetting("INITIAL GUESS STRATEGY", "LAST")) {
    DLinearSolver.SetupInitialGuess<InitialGuess::Last<dfloat> >(Ndofs, platform, settings, comm);
  } else if (settings.compareSetting("INITIAL GUESS STRATEGY", "ZERO")) {
    DLinearSolver.SetupInitialGuess<InitialGuess::Zero<dfloat> >(Ndofs, platform, settings, comm);
  } else if (settings.compareSetting("INITIAL GUESS STRATEGY", "CLASSIC")) {
    DLinearSolver.SetupInitialGuess<InitialGuess::ClassicProjection<dfloat> >(Ndofs, platform, settings, comm);
  } else if (settings.compareSetting("INITIAL GUESS STRATEGY", "QR")) {
    DLinearSolver.SetupInitialGuess<InitialGuess::RollingQRProjection<dfloat> >(Ndofs, platform, settings, comm);
  } else if (settings.compareSetting("INITIAL GUESS STRATEGY", "EXTRAP")) {
    DLinearSolver.SetupInitialGuess<InitialGuess::Extrap<dfloat> >(Ndofs, platform, settings, comm);
  }

  //setup linear solver
  hlong NglobalDofs;
  if (disc_c0){
    NglobalDofs = ogsMasked.NgatherGlobal*Nfields;
  } else {
    NglobalDofs = mesh.NelementsGlobal*mesh.Np*Nfields;
  }

  dlong Nall = mesh.Np*(mesh.Nelements+mesh.totalHaloPairs);

  
  properties_t kernelInfo = mesh.props; //copy base occa properties

  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  //add standard boundary functions
  std::string boundaryHeaderFileName;
  if (mesh.dim==2)
    boundaryHeaderFileName = std::string(DELLIPTIC "/data/ellipticBoundary2D.h");
  else if (mesh.dim==3)
    boundaryHeaderFileName = std::string(DELLIPTIC "/data/ellipticBoundary3D.h");
  kernelInfo["includes"] += boundaryHeaderFileName;

  // set kernel name suffix
  std::string suffix;
  if(mesh.elementType==Mesh::TRIANGLES) {
    suffix = "Tri2D";
  } else if(mesh.elementType==Mesh::QUADRILATERALS) {
    suffix = "Quad2D";
  } else if(mesh.elementType==Mesh::TETRAHEDRA) {
    suffix = "Tet3D";
  } else { //mesh.elementType==Mesh::HEXAHEDRA)
    suffix = "Hex3D";
  }

  std::string oklFilePrefix = DELLIPTIC "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  // WE WILL ASSUME ZERO DIRICHLET BOUNDARIES
  
  fileName   = oklFilePrefix + "waveKernels"  + oklFileSuffix;

  
  kernelName = "waveStageUpdate";
  kernel_t waveStageUpdateKernel =
     platform.buildKernel(fileName, kernelName, kernelInfo);


  kernelName = "waveCombine";
  kernel_t waveCombineKernel =
     platform.buildKernel(fileName, kernelName, kernelInfo);


  kernelName = "waveErrorEstimate";
  kernel_t waveErrorEstimateKernel
     = platform.buildKernel(fileName, kernelName, kernelInfo);
  
  fileName   = oklFilePrefix + "waveKernels" + suffix + oklFileSuffix;

  kernelName = "waveStepInitialize" + suffix;
  kernel_t waveStepInitializeKernel =
     platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "waveStageFinalize" + suffix;
  kernel_t waveStageFinalizeKernel =
     platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "waveStepFinalize" + suffix;
  kernel_t waveStepFinalizeKernel =
     platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "waveInitialConditions" + suffix;
  kernel_t waveInitialConditionsKernel
     = platform.buildKernel(fileName, kernelName, kernelInfo);

  //create occa buffers
  memory<dfloat> DL(Nall), PL(Nall), DrhsL(Nall), PhatL(Nall), DhatL(Nall*Nstages);

  deviceMemory<dfloat> o_DL    = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_PL    = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_newDL    = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_newPL    = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_embDL    = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_embPL    = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_DtildeL  = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_DrhsL = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_DhatL = platform.malloc<dfloat>(Nall*Nstages);
  deviceMemory<dfloat> o_PhatL = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_scratch1L = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_scratch2L = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_errL    = platform.malloc<dfloat>(Nall);

  deviceMemory<dfloat> o_Dtilde   = platform.malloc<dfloat>(NglobalDofs);
  deviceMemory<dfloat> o_Drhs     = platform.malloc<dfloat>(NglobalDofs);
  deviceMemory<dfloat> o_scratch1 = platform.malloc<dfloat>(NglobalDofs);
  deviceMemory<dfloat> o_scratch2 = platform.malloc<dfloat>(NglobalDofs);
  
  memory<dfloat> invMM, V, MM;
  
  if(mesh.elementType==Mesh::TRIANGLES) {
    mesh.VandermondeTri2D(mesh.N, mesh.r, mesh.s, V);
    mesh.invMassMatrixTri2D(mesh.Np, V, invMM);
    mesh.MassMatrixTri2D(mesh.Np, V, MM);
  } else if(mesh.elementType==Mesh::QUADRILATERALS) {
    invMM.malloc(1);
    MM.malloc(1);
  } else if(mesh.elementType==Mesh::TETRAHEDRA) {
    mesh.VandermondeTet3D(mesh.N, mesh.r, mesh.s, mesh.t, V);
    mesh.invMassMatrixTet3D(mesh.Np, V, invMM);
    mesh.MassMatrixTet3D(mesh.Np, V, MM);
  } else { 
    printf("HEXES NOT IMPLEMENTED YET\n");
    exit(-1);
  }
  
  deviceMemory<dfloat> o_invMM = platform.malloc<dfloat>(mesh.Np*mesh.Np, invMM);
  deviceMemory<dfloat> o_MM    = platform.malloc<dfloat>(mesh.Np*mesh.Np, MM);

  // triangle specific
  memory<dfloat> invWJ;
  memory<dfloat> WJ(mesh.Nelements);
  deviceMemory<dfloat> o_invWJ;
  deviceMemory<dfloat> o_WJ;
  if(mesh.elementType==Mesh::TRIANGLES ||
     mesh.elementType==Mesh::TETRAHEDRA) {

    WJ.malloc(mesh.Nelements);
    invWJ.malloc(mesh.Nelements);
    for(int e=0;e<mesh.Nelements;++e){
      invWJ[e] = 1./mesh.wJ[e];
      WJ[e] = mesh.wJ[e];
    }

    o_invWJ = platform.malloc<dfloat>(mesh.Nelements, invWJ);
    o_WJ    = platform.malloc<dfloat>(mesh.Nelements, WJ);  

  }else{

    // use ungathered weights in WJ on device
    WJ.malloc(mesh.Np*mesh.Nelements);
    for(int n=0;n<mesh.Np*mesh.Nelements;++n){
      WJ[n] = mesh.wJ[n];
    }
    o_WJ    = platform.malloc<dfloat>(mesh.Nelements*mesh.Np, WJ);
    
    // gather weights for C0
    if(disc_c0){
      ogsMasked.GatherScatter(WJ, 1, ogs::Add, ogs::Sym);
    }
    
    // use globalized mass for C0
    invWJ.malloc(mesh.Np*mesh.Nelements);
    for(int n=0;n<mesh.Np*mesh.Nelements;++n){
      invWJ[n] = (WJ[n]) ? 1./WJ[n]:0;
    }
    
    o_invWJ = platform.malloc<dfloat>(mesh.Nelements*mesh.Np, invWJ);
  }
  
  waveInitialConditionsKernel(Nall, t, mesh.o_x, mesh.o_y, mesh.o_z, o_DL, o_PL);

  stoppingCriteria_t<dfloat> *stoppingCriteria = NULL;
  ellipticStoppingCriteria<dfloat> *esc = NULL;

  if(settings.compareSetting("STOPPING CRITERIA", "ERRORESTIMATE")){
    printf("SETTING UP ESC\n");
    esc = new ellipticStoppingCriteria<dfloat>(this, NULL);
    esc->reset();
    stoppingCriteria = esc;
  }
  else
  {
    stoppingCriteria = new stoppingCriteria_t<dfloat>();
  }

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

      lambda = lambdaSolve;

      // record local RHS
      if(esc)
         esc->setLocalRHS(o_DrhsL);

      // gather rhs to globalDofs if c0
      if(disc_c0){
        ogsMasked.Gather(o_Drhs,   o_DrhsL,   1, ogs::Add, ogs::Trans);

#if 0        
        {
          dfloat val = platform.linAlg().norm2(NglobalDofs, o_Drhs, mesh.comm);
          std::cout << "|Drhs|=" << val << std::endl;
        }
#endif
        
        // TW not sure about this (working assumption is that DtildeL is already implicitly continuous)
        ogsMasked.Gather(o_Dtilde, o_DtildeL, 1, ogs::Add, ogs::NoTrans);
      }else{
        o_DtildeL.copyTo(o_Dtilde);
        o_DrhsL.copyTo(o_Drhs);
      }
      
      int iterD = Solve(DLinearSolver, o_Dtilde, o_Drhs, tol, maxIter, verbose, stoppingCriteria);

      //add the boundary data to the masked nodes
      if(disc_c0){

#if 0
        {
          dfloat val = platform.linAlg().norm2(NglobalDofs, o_Dtilde, mesh.comm);
          std::cout << "|Dtilde|=" << val << std::endl;
        }
#endif
        
        // scatter x to LocalDofs if c0
        ogsMasked.Scatter(o_DtildeL, o_Dtilde, 1, ogs::NoTrans);

#if 0
        {
          dfloat val = platform.linAlg().norm2(Nall, o_DtildeL, mesh.comm);
          std::cout << "|DtildeL|=" << val << std::endl;
        }
#endif
        
      }else{
        o_Dtilde.copyTo(o_DtildeL);
      }
      
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

    // KERNEL 4: FINALIZE
    // P = Phat(:,1) +     dt*(Dhat(:,1:Nstages)*beta'); 
    // D = Dhat(:,1) + dt*LAP*(Phat(:,1:Nstages)*beta');
    
    // a. Phat(:,1:Nstages)*beta' => o_scratchL
    waveCombineKernel(Nall, Nstages, dt, o_beta, o_betaAlpha, o_PhatL, o_DhatL, o_scratch1L);

    // b. L*(Phat(:,1:Nstages)*beta') => o_DL
    lambda = 0;
    if(disc_c0){
      // gather up RHS 
      ogsMasked.Gather(o_scratch1, o_scratch1L, 1, ogs::Add, ogs::NoTrans);
      Operator(o_scratch1, o_scratch2);
      ogsMasked.Scatter(o_scratch2L, o_scratch2, 1, ogs::NoTrans);
    }else{
      Operator(o_scratch1L, o_scratch2L);
    }
    lambda = lambdaSolve;
    
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
        sprintf(fname, "P_%04d_%04d.vtu",  mesh.rank, slice);
        PlotFields(PL, fname);

        // output field files
        sprintf(fname, "D_%04d_%04d.vtu",  mesh.rank, slice);
        PlotFields(DL, fname);

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
