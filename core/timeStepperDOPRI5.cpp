/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "core.hpp"
#include "timeStepper.hpp"

namespace TimeStepper {

dopri5::dopri5(dlong Nelements, dlong NhaloElements,
               int Np, int Nfields, solver_t& _solver):
  timeStepper_t(Nelements, NhaloElements, Np, Nfields, _solver) {

  Nrk = 7;

  o_rhsq   = device.malloc(N*sizeof(dfloat));
  o_rkq    = device.malloc((N+Nhalo)*sizeof(dfloat));
  o_rkrhsq = device.malloc(Nrk*N*sizeof(dfloat));
  o_rkerr  = device.malloc(N*sizeof(dfloat));

  o_saveq  = device.malloc(N*sizeof(dfloat));

  Nblock = (N+BLOCKSIZE-1)/BLOCKSIZE;
  errtmp = (dfloat*) occaHostMallocPinned(device, Nblock*sizeof(dfloat),
                                       NULL, o_errtmp, h_errtmp);

  hlong Nlocal = N;
  hlong Ntotal;
  MPI_Allreduce(&Nlocal, &Ntotal, 1, MPI_HLONG, MPI_SUM, comm);

  occa::properties kernelInfo = props; //copy base occa properties from solver

  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)BLOCKSIZE;

  rkUpdateKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                    "timeStepperDOPRI5.okl",
                                    "dopri5RkUpdate",
                                    kernelInfo, comm);

  rkStageKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                    "timeStepperDOPRI5.okl",
                                    "dopri5RkStage",
                                    kernelInfo, comm);

  rkErrorEstimateKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                    "timeStepperDOPRI5.okl",
                                    "dopri5ErrorEstimate",
                                    kernelInfo, comm);

  // Dormand Prince -order (4) 5 with PID timestep control
  dfloat _rkC[7] = {0.0, 0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0};
  dfloat _rkA[7*7] ={             0.0,             0.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                                 0.2,             0.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                            3.0/40.0,        9.0/40.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                           44.0/45.0,      -56.0/15.0,       32.0/9.0,          0.0,             0.0,       0.0, 0.0,
                      19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0,             0.0,       0.0, 0.0,
                       9017.0/3168.0,     -355.0/33.0, 46732.0/5247.0,   49.0/176.0, -5103.0/18656.0,       0.0, 0.0,
                          35.0/384.0,             0.0,   500.0/1113.0,  125.0/192.0,  -2187.0/6784.0, 11.0/84.0, 0.0 };
  dfloat _rkE[7] = {71.0/57600.0,  0.0, -71.0/16695.0, 71.0/1920.0, -17253.0/339200.0, 22.0/525.0, -1.0/40.0 };

  rkC = (dfloat*) calloc(Nrk, sizeof(dfloat));
  rkE = (dfloat*) calloc(Nrk, sizeof(dfloat));
  rkA = (dfloat*) calloc(Nrk*Nrk, sizeof(dfloat));

  memcpy(rkC, _rkC, Nrk*sizeof(dfloat));
  memcpy(rkE, _rkE, Nrk*sizeof(dfloat));
  memcpy(rkA, _rkA, Nrk*Nrk*sizeof(dfloat));

  o_rkA = device.malloc(Nrk*Nrk*sizeof(dfloat), rkA);
  o_rkE = device.malloc(Nrk*sizeof(dfloat), rkE);

  dtMIN = 1E-9; //minumum allowed timestep
  ATOL = 1E-6;  //absolute error tolerance
  RTOL = 1E-6;  //relative error tolerance
  safe = 0.8;   //safety factor

  //error control parameters
  beta = 0.05;
  factor1 = 0.2;
  factor2 = 10.0;

  exp1 = 0.2 - 0.75*beta;
  invfactor1 = 1.0/factor1;
  invfactor2 = 1.0/factor2;
  facold = 1E-4;
  sqrtinvNtotal = 1.0/sqrt(Ntotal);
}

void dopri5::Run(occa::memory &o_q, dfloat start, dfloat end) {

  dfloat time = start;

  int rank;
  MPI_Comm_rank(comm, &rank);

  solver.Report(time,0);

  dfloat outputInterval;
  settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat outputTime = time + outputInterval;

  int tstep=0, allStep=0;

  while (time < end) {

    if (dt<dtMIN){
      stringstream ss;
      ss << "Time step became too small at time step = " << tstep;
      LIBP_ABORT(ss.str());
    }
    if (std::isnan(dt)) {
      stringstream ss;
      ss << "Solution became unstable at time step = " << tstep;
      LIBP_ABORT(ss.str());
    }

    //check for final timestep
    if (time+dt > end){
      dt = end-time;
    }

    Step(o_q, time, dt);

    // compute Dopri estimator
    dfloat err = Estimater(o_q);

    // build controller
    dfloat fac1 = pow(err,exp1);
    dfloat fac = fac1/pow(facold,beta);

    fac = mymax(invfactor2, mymin(invfactor1,fac/safe));
    dfloat dtnew = dt/fac;

    if (err<1.0) { //dt is accepted

      // check for output during this step and do a mini-step
      if (time<outputTime && time+dt>=outputTime) {
        dfloat savedt = dt;

        // save rkq
        Backup(o_rkq);

        // change dt to match output
        dt = outputTime-time;

        // if(!rank)
        //   printf("Taking output mini step: %g\n", dt);

        // time step to output
        Step(o_q, time, dt);

        // shift for output
        o_rkq.copyTo(o_q);

        // output  (print from rkq)
        // if (!rank) printf("\n");
        solver.Report(outputTime,tstep);

        // restore time step
        dt = savedt;

        // increment next output time
        outputTime += outputInterval;

        // accept saved rkq
        Restore(o_rkq);
        AcceptStep(o_q, o_rkq);
      } else {
        // accept rkq
        AcceptStep(o_q, o_rkq);
      }

      time += dt;
      while (time>outputTime) outputTime+= outputInterval; //catch up next output in case dt>outputInterval

      facold = mymax(err,1E-4); // hard coded factor ?

      // if (!rank)
      //   printf("\r time = %g (%d), dt = %g accepted                      ", time, allStep,  dt);

      tstep++;
    } else {
      dtnew = dt/(mymax(invfactor1,fac1/safe));

      // if (!rank)
      //   printf("\r time = %g (%d), dt = %g rejected, trying %g", time, allStep, dt, dtnew);
      if (!rank)
        printf("Repeating timestep %d. dt was %g, trying %g.\n", tstep, dt, dtnew);
    }
    dt = dtnew;
    allStep++;
  }

  if (!rank)
    printf("%d accepted steps and %d total steps\n", tstep, allStep);
}

void dopri5::Backup(occa::memory &o_Q) {
  o_saveq.copyFrom(o_Q, N*sizeof(dfloat));
}

void dopri5::Restore(occa::memory &o_Q) {
  o_saveq.copyTo(o_Q, N*sizeof(dfloat));
}

void dopri5::AcceptStep(occa::memory &o_q, occa::memory &o_rq) {
  o_q.copyFrom(o_rq, N*sizeof(dfloat));
}

void dopri5::Step(occa::memory &o_q, dfloat time, dfloat _dt) {

  //RK step
  for(int rk=0;rk<Nrk;++rk){

    // t_rk = t + C_rk*_dt
    dfloat currentTime = time + rkC[rk]*_dt;

    //compute RK stage
    // rkq = q + _dt sum_{i=0}^{rk-1} a_{rk,i}*rhsq_i
    rkStageKernel(N,
                  rk,
                  _dt,
                  o_rkA,
                  o_q,
                  o_rkrhsq,
                  o_rkq);

    //evaluate ODE rhs = f(q,t)
    solver.rhsf(o_rkq, o_rhsq, currentTime);

    // update solution using Runge-Kutta
    // rkrhsq_rk = rhsq
    // if rk==6
    //   q = q + _dt*sum_{i=0}^{rk} rkA_{rk,i}*rkrhs_i
    //   rkerr = _dt*sum_{i=0}^{rk} rkE_{rk,i}*rkrhs_i
    rkUpdateKernel(N,
                   rk,
                   _dt,
                   o_rkA,
                   o_rkE,
                   o_q,
                   o_rhsq,
                   o_rkrhsq,
                   o_rkq,
                   o_rkerr);
  }
}

dfloat dopri5::Estimater(occa::memory& o_q){

  //Error estimation
  //E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
  //      DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
  rkErrorEstimateKernel(N,
                        ATOL,
                        RTOL,
                        o_q,
                        o_rkq,
                        o_rkerr,
                        o_errtmp);

  o_errtmp.copyTo(errtmp);
  dfloat localerr = 0;
  dfloat err = 0;
  for(dlong n=0;n<Nblock;++n){
    localerr += errtmp[n];
  }
  MPI_Allreduce(&localerr, &err, 1, MPI_DFLOAT, MPI_SUM, comm);

  err = sqrt(err)*sqrtinvNtotal;

  return err;
}

dopri5::~dopri5() {
  if (o_rkq.size()) o_rkq.free();
  if (o_rkrhsq.size()) o_rkrhsq.free();
  if (o_rkerr.size()) o_rkerr.free();
  if (o_errtmp.size()) o_errtmp.free();
  if (o_rkA.size()) o_rkA.free();
  if (o_rkE.size()) o_rkE.free();

  if (rkC) free(rkC);
  if (rkA) free(rkA);
  if (rkE) free(rkE);

  rkUpdateKernel.free();
  rkStageKernel.free();
  rkErrorEstimateKernel.free();
}

/**************************************************/
/* PML version                                    */
/**************************************************/

dopri5_pml::dopri5_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
                      int Np, int Nfields, int Npmlfields, solver_t& _solver):
  dopri5(Nelements, NhaloElements, Np, Nfields, _solver),
  Npml(Npmlfields*Np*NpmlElements) {

  if (Npml) {
    dfloat *pmlq = (dfloat *) calloc(Npml,sizeof(dfloat));
    o_pmlq   = device.malloc(Npml*sizeof(dfloat), pmlq);
    free(pmlq);

    o_rhspmlq   = device.malloc(Npml*sizeof(dfloat));
    o_rkpmlq    = device.malloc(Npml*sizeof(dfloat));
    o_rkrhspmlq = device.malloc(Nrk*Npml*sizeof(dfloat));

    o_savepmlq  = device.malloc(Npml*sizeof(dfloat));

    occa::properties kernelInfo = props; //copy base occa properties from solver

    //add defines
    kernelInfo["defines/" "p_blockSize"] = (int)BLOCKSIZE;

    rkPmlUpdateKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                      "timeStepperDOPRI5.okl",
                                      "dopri5RkPmlUpdate",
                                      kernelInfo, comm);
  }
}

void dopri5_pml::Backup(occa::memory &o_Q) {
  o_saveq.copyFrom(o_Q, N*sizeof(dfloat));
  if (Npml)
    o_savepmlq.copyFrom(o_rkpmlq, Npml*sizeof(dfloat));
}

void dopri5_pml::Restore(occa::memory &o_Q) {
  o_saveq.copyTo(o_Q, N*sizeof(dfloat));
  if (Npml)
    o_savepmlq.copyTo(o_rkpmlq, Npml*sizeof(dfloat));
}

void dopri5_pml::AcceptStep(occa::memory &o_q, occa::memory &o_rq) {
  o_q.copyFrom(o_rq, N*sizeof(dfloat));
  if (Npml)
    o_pmlq.copyFrom(o_rkpmlq, Npml*sizeof(dfloat));
}

void dopri5_pml::Step(occa::memory &o_q, dfloat time, dfloat _dt) {

  //RK step
  for(int rk=0;rk<Nrk;++rk){

    // t_rk = t + C_rk*_dt
    dfloat currentTime = time + rkC[rk]*_dt;

    //compute RK stage
    // rkq = q + _dt sum_{i=0}^{rk-1} a_{rk,i}*rhsq_i
    rkStageKernel(N,
                  rk,
                  _dt,
                  o_rkA,
                  o_q,
                  o_rkrhsq,
                  o_rkq);
    if (Npml)
      rkStageKernel(Npml,
                    rk,
                    _dt,
                    o_rkA,
                    o_pmlq,
                    o_rkrhspmlq,
                    o_rkpmlq);

    //evaluate ODE rhs = f(q,t)
    solver.rhsf_pml(o_rkq, o_rkpmlq, o_rhsq, o_rhspmlq, currentTime);

    // update solution using Runge-Kutta
    // rkrhsq_rk = rhsq
    // if rk==6
    //   q = q + _dt*sum_{i=0}^{rk} rkA_{rk,i}*rkrhs_i
    //   rkerr = _dt*sum_{i=0}^{rk} rkE_{rk,i}*rkrhs_i
    rkUpdateKernel(N,
                   rk,
                   _dt,
                   o_rkA,
                   o_rkE,
                   o_q,
                   o_rhsq,
                   o_rkrhsq,
                   o_rkq,
                   o_rkerr);
    if (Npml)
      rkPmlUpdateKernel(Npml,
                     rk,
                     _dt,
                     o_rkA,
                     o_pmlq,
                     o_rhspmlq,
                     o_rkrhspmlq,
                     o_rkpmlq);
  }
}

dopri5_pml::~dopri5_pml() {
  if (o_pmlq.size()) o_pmlq.free();
  if (o_rkpmlq.size()) o_rkpmlq.free();
  if (o_rhspmlq.size()) o_rhspmlq.free();
  if (o_rkrhspmlq.size()) o_rkrhspmlq.free();
  if (o_savepmlq.size()) o_savepmlq.free();
}

} //namespace TimeStepper
