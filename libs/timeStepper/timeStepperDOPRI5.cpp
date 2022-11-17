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

#include "core.hpp"
#include "timeStepper.hpp"

namespace libp {

namespace TimeStepper {

constexpr int blocksize = 256;

dopri5::dopri5(dlong Nelements, dlong NhaloElements,
               int Np, int Nfields,
               platform_t& _platform, comm_t _comm):
  dopri5(Nelements, 0, NhaloElements,
         Np, Nfields, 0,
         _platform, _comm) {}

dopri5::dopri5(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
               int Np, int Nfields, int Npmlfields,
               platform_t& _platform, comm_t _comm):
  timeStepperBase_t(Nelements, NpmlElements, NhaloElements,
                    Np, Nfields, Npmlfields,
                    _platform, _comm) {

  Nrk = 7;

  hlong Ntotal = N;
  comm.Allreduce(Ntotal);

  properties_t kernelInfo = platform.props(); //copy base occa properties from solver

  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)blocksize;

  rkUpdateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperDOPRI5.okl",
                                    "dopri5RkUpdate",
                                    kernelInfo);

  rkPmlUpdateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                      "timeStepperDOPRI5.okl",
                                      "dopri5RkPmlUpdate",
                                      kernelInfo);

  rkStageKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperDOPRI5.okl",
                                    "dopri5RkStage",
                                    kernelInfo);

  rkErrorEstimateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperDOPRI5.okl",
                                    "dopri5ErrorEstimate",
                                    kernelInfo);

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

  rkC.malloc(Nrk);
  rkE.malloc(Nrk);
  rkA.malloc(Nrk*Nrk);

  rkC.copyFrom(_rkC);
  rkE.copyFrom(_rkE);
  rkA.copyFrom(_rkA);

  o_rkA = platform.malloc<dfloat>(rkA);
  o_rkE = platform.malloc<dfloat>(rkE);

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

void dopri5::Run(solver_t& solver,
                 deviceMemory<dfloat> o_q,
                 std::optional<deviceMemory<dfloat>> o_pmlq,
                 dfloat start, dfloat end) {

  /*Pre-reserve memory pool space to avoid some unnecessary re-sizing*/
  platform.reserve<dfloat>((N+Nhalo) + (Nrk+2)*N + (Nrk+2)*Npml
                           + 7 * platform.memPoolAlignment<dfloat>());

  deviceMemory<dfloat> o_rkq    = platform.reserve<dfloat>(N+Nhalo);
  deviceMemory<dfloat> o_rkpmlq = platform.reserve<dfloat>(Npml);
  deviceMemory<dfloat> o_rkerr  = platform.reserve<dfloat>(N);

  dfloat time = start;

  solver.Report(time,0);

  dfloat outputInterval=0.0;
  solver.settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat outputTime = time + outputInterval;

  int tstep=0, allStep=0;

  while (time < end) {

    LIBP_ABORT("Time step became too small at time step = " << tstep,
               dt<dtMIN);
    LIBP_ABORT("Solution became unstable at time step = " << tstep,
               std::isnan(dt));

    //check for final timestep
    if (time+dt > end){
      dt = end-time;
    }

    Step(solver, o_q, o_pmlq,
         o_rkq, o_rkpmlq, o_rkerr,
         time, dt);

    // compute Dopri estimator
    dfloat err = Estimater(o_q, o_rkq, o_rkerr);

    // build controller
    dfloat fac1 = pow(err,exp1);
    dfloat fac = fac1/pow(facold,beta);

    fac = std::max(invfactor2, std::min(invfactor1,fac/safe));
    dfloat dtnew = dt/fac;

    if (err<1.0) { //dt is accepted

      // check for output during this step and do a mini-step
      if (time<outputTime && time+dt>=outputTime) {
        dfloat savedt = dt;

        // save rkq
        deviceMemory<dfloat> o_saveq = platform.reserve<dfloat>(N);
        deviceMemory<dfloat> o_savepmlq  = platform.reserve<dfloat>(Npml);
        o_saveq.copyFrom(o_rkq, N, properties_t("async", true));
        if (o_pmlq.has_value()) {
          o_savepmlq.copyFrom(o_rkpmlq, Npml, properties_t("async", true));
        }

        // change dt to match output
        dt = outputTime-time;

        // time step to output
        Step(solver, o_q, o_pmlq,
             o_rkq, o_rkpmlq, o_rkerr,
             time, dt);

        // shift for output
        o_rkq.copyTo(o_q, properties_t("async", true));
        if (o_pmlq.has_value()) {
          o_rkpmlq.copyTo(o_pmlq.value(), properties_t("async", true));
        }

        // output  (print from rkq)
        solver.Report(outputTime,tstep);

        // restore time step
        dt = savedt;

        // increment next output time
        outputTime += outputInterval;

        // accept saved rkq
        o_saveq.copyTo(o_q, N, properties_t("async", true));
        if (o_pmlq.has_value()) {
          o_savepmlq.copyTo(o_pmlq.value(), Npml, properties_t("async", true));
        }
      } else {
        // accept rkq
        o_q.copyFrom(o_rkq, N, properties_t("async", true));
        if (o_pmlq.has_value()) {
          o_pmlq.value().copyFrom(o_rkpmlq, Npml, properties_t("async", true));
        }
      }

      time += dt;
      while (time>outputTime) outputTime+= outputInterval; //catch up next output in case dt>outputInterval

      constexpr dfloat errMax = 1.0e-4;  // hard coded factor ?
      facold = std::max(err,errMax);

      tstep++;
    } else {
      dtnew = dt/(std::max(invfactor1,fac1/safe));
    }
    dt = dtnew;
    allStep++;
  }
}


void dopri5::Step(solver_t& solver,
                  deviceMemory<dfloat> o_q,
                  std::optional<deviceMemory<dfloat>> o_pmlq,
                  deviceMemory<dfloat> o_rkq,
                  deviceMemory<dfloat> o_rkpmlq,
                  deviceMemory<dfloat> o_rkerr,
                  dfloat time, dfloat _dt) {

  deviceMemory<dfloat> o_rhsq   = platform.reserve<dfloat>(N);
  deviceMemory<dfloat> o_rkrhsq = platform.reserve<dfloat>(Nrk*N);
  deviceMemory<dfloat> o_rhspmlq   = platform.reserve<dfloat>(Npml);
  deviceMemory<dfloat> o_rkrhspmlq = platform.reserve<dfloat>(Nrk*Npml);

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
    if (o_pmlq.has_value()) {
      rkStageKernel(Npml,
                    rk,
                    _dt,
                    o_rkA,
                    o_pmlq.value(),
                    o_rkrhspmlq,
                    o_rkpmlq);
    }

    //evaluate ODE rhs = f(q,t)
    if (o_pmlq.has_value()) {
      solver.rhsf_pml(o_rkq, o_rkpmlq, o_rhsq, o_rhspmlq, currentTime);
    } else {
      solver.rhsf(o_rkq, o_rhsq, currentTime);
    }

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
    if (o_pmlq.has_value()) {
      rkPmlUpdateKernel(Npml,
                       rk,
                       _dt,
                       o_rkA,
                       o_pmlq.value(),
                       o_rhspmlq,
                       o_rkrhspmlq,
                       o_rkpmlq);
    }
  }
}

dfloat dopri5::Estimater(deviceMemory<dfloat>& o_q,
                         deviceMemory<dfloat>& o_rkq,
                         deviceMemory<dfloat>& o_rkerr){

  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  pinnedMemory<dfloat> h_errtmp = platform.hostReserve<dfloat>(1);
  deviceMemory<dfloat> o_errtmp = platform.reserve<dfloat>(blocksize);

  //Error estimation
  //E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
  //      DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
  rkErrorEstimateKernel(Nblock,
                        N,
                        ATOL,
                        RTOL,
                        o_q,
                        o_rkq,
                        o_rkerr,
                        o_errtmp);

  dfloat err;
  if (Nblock>0) {
    h_errtmp.copyFrom(o_errtmp, 1, 0, properties_t("async", true));
    platform.finish();
    err = h_errtmp[0];
  } else {
    err = 0.0;
  }

  comm.Allreduce(err);

  err = sqrt(err)*sqrtinvNtotal;

  return err;
}

} //namespace TimeStepper

} //namespace libp
