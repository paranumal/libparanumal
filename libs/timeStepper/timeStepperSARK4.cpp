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

#include <math.h>
#include "solver.hpp"
#include "timeStepper.hpp"
#include <complex>

namespace libp {

namespace TimeStepper {

using std::complex;

static void UpdateCoefficients(const int Nfields,
                               const memory<dfloat> &lambda,
                               const dfloat dt,
                               pinnedMemory<dfloat> &h_rkX,
                               pinnedMemory<dfloat> &h_rkA,
                               pinnedMemory<dfloat> &h_rkE);

sark4::sark4(dlong _Nelements, dlong _NhaloElements,
             int _Np, int _Nfields,
             memory<dfloat> _lambda,
             platform_t& _platform, comm_t _comm):
  timeStepperBase_t(_Nelements, _NhaloElements, _Np, _Nfields,
                    _platform, _comm),
  Np(_Np),
  Nfields(_Nfields),
  Nelements(_Nelements),
  NhaloElements(_NhaloElements) {

  lambda.malloc(Nfields);
  lambda.copyFrom(_lambda);

  //Nrk = 5;
  order = 4;
  embeddedOrder = 3;

  dlong Nlocal = Nelements*Np*Nfields;
  dlong Ntotal = (Nelements+NhaloElements)*Np*Nfields;

  o_rkq    = platform.malloc<dfloat>(Ntotal);
  o_rhsq   = platform.malloc<dfloat>(Nlocal);
  o_rkrhsq = platform.malloc<dfloat>(Nlocal*Nrk);
  o_rkerr  = platform.malloc<dfloat>(Nlocal);

  o_saveq  = platform.malloc<dfloat>(Nlocal);

  const int blocksize=256;

  Nblock = (N+blocksize-1)/blocksize;
  h_errtmp = platform.hostMalloc<dfloat>(Nblock);
  o_errtmp = platform.malloc<dfloat>(Nblock);

  hlong gNtotal = Nlocal;
  comm.Allreduce(gNtotal);

  //copy base occa properties from platform
  properties_t kernelInfo = platform.props();

  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)blocksize;
  kernelInfo["defines/" "p_Nrk"]     = (int)Nrk;
  kernelInfo["defines/" "p_Np"]      = (int)Np;
  kernelInfo["defines/" "p_Nfields"] = (int)Nfields;

  rkUpdateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperSARK.okl",
                                    "sarkRkUpdate",
                                    kernelInfo);

  rkStageKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperSARK.okl",
                                    "sarkRkStage",
                                    kernelInfo);

  rkErrorEstimateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperSARK.okl",
                                    "sarkErrorEstimate",
                                    kernelInfo);

  // Semi-Analytic Runge Kutta - order (3) 4 with PID timestep control
  dfloat _rkC[Nrk] = {0.0, 0.5, 0.5, 1.0, 1.0};
  rkC.malloc(Nrk);
  rkC.copyFrom(_rkC);

  h_rkX = platform.hostMalloc<dfloat>(Nfields*Nrk);
  h_rkA = platform.hostMalloc<dfloat>(Nfields*Nrk*Nrk);
  h_rkE = platform.hostMalloc<dfloat>(Nfields*Nrk);

  o_rkX = platform.malloc<dfloat>(Nfields*Nrk);
  o_rkA = platform.malloc<dfloat>(Nfields*Nrk*Nrk);
  o_rkE = platform.malloc<dfloat>(Nfields*Nrk);

  dtMIN = 1E-9; //minumum allowed timestep
  ATOL = 1E-5;  //absolute error tolerance
  RTOL = 1E-4;  //relative error tolerance
  safe = 0.8;   //safety factor

  //error control parameters
  beta = 0.05;
  factor1 = 0.2;
  factor2 = 10.0;

  exp1 = 1.0/(embeddedOrder+1) - 0.75*beta;
  invfactor1 = 1.0/factor1;
  invfactor2 = 1.0/factor2;
  facold = 1E-4;
  sqrtinvNtotal = 1.0/sqrt(gNtotal);
}

void sark4::Run(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat start, dfloat end) {

  dfloat time = start;

  solver.Report(time,0);

  dfloat outputInterval=0.0;
  solver.settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat outputTime = time + outputInterval;

  int tstep=0, allStep=0;

  //Compute Butcher Tableau
  UpdateCoefficients(Nfields, lambda, dt, h_rkX, h_rkA, h_rkE);

  // move data to platform
  h_rkX.copyTo(o_rkX, properties_t("async", true));
  h_rkA.copyTo(o_rkA, properties_t("async", true));
  h_rkE.copyTo(o_rkE, properties_t("async", true));

  while (time < end) {

    LIBP_ABORT("Time step became too small at time step = " << tstep,
               dt<dtMIN);
    LIBP_ABORT("Solution became unstable at time step = " << tstep,
               std::isnan(dt));

    //check for final timestep
    if (time+dt > end){
      dt = end-time;
    }

    Step(solver, o_q, time, dt);

    // compute Dopri estimator
    dfloat err = Estimater(o_q);

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
        Backup(o_rkq);

        // change dt to match output
        dt = outputTime-time;

        //Compute Butcher Tableau
        UpdateCoefficients(Nfields, lambda, dt, h_rkX, h_rkA, h_rkE);

        // move data to platform
        h_rkX.copyTo(o_rkX, properties_t("async", true));
        h_rkA.copyTo(o_rkA, properties_t("async", true));
        h_rkE.copyTo(o_rkE, properties_t("async", true));

        // time step to output
        Step(solver, o_q, time, dt);

        // shift for output
        o_rkq.copyTo(o_q, properties_t("async", true));

        // output  (print from rkq)
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

      constexpr dfloat errMax = 1.0e-4;  // hard coded factor ?
      facold = std::max(err,errMax);

      tstep++;
    } else {
      dtnew = dt/(std::max(invfactor1,fac1/safe));
    }
    dt = dtnew;

    //Compute Butcher Tableau
    UpdateCoefficients(Nfields, lambda, dt, h_rkX, h_rkA, h_rkE);

    // move data to platform
    h_rkX.copyTo(o_rkX, properties_t("async", true));
    h_rkA.copyTo(o_rkA, properties_t("async", true));
    h_rkE.copyTo(o_rkE, properties_t("async", true));

    allStep++;
  }
}

void sark4::Backup(deviceMemory<dfloat> &o_q) {
  o_saveq.copyFrom(o_q, N, properties_t("async", true));
}

void sark4::Restore(deviceMemory<dfloat> &o_q) {
  o_saveq.copyTo(o_q, N, properties_t("async", true));
}

void sark4::AcceptStep(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rq) {
  o_q.copyFrom(o_rq, N, properties_t("async", true));
}

void sark4::Step(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat time, dfloat _dt) {

  //RK step
  for(int rk=0;rk<Nrk;++rk){

    // t_rk = t + C_rk*_dt
    dfloat currentTime = time + rkC[rk]*_dt;

    //compute RK stage
    // rkq = x_{rk}*q + _dt sum_{i=0}^{rk-1} a_{rk,i}*rhsq_i
    rkStageKernel(Nelements,
                  rk,
                  _dt,
                  o_rkX,
                  o_rkA,
                  o_q,
                  o_rkrhsq,
                  o_rkq);

    //evaluate ODE rhs = f(q,t)
    solver.rhsf(o_rkq, o_rhsq, currentTime);

    // update solution using Runge-Kutta
    // rkrhsq_rk = rhsq
    // if rk==6
    //   q = rkX_{rk}*q + _dt*sum_{i=0}^{rk} rkA_{rk,i}*rkrhs_i
    //   rkerr = _dt*sum_{i=0}^{rk} rkE_{rk,i}*rkrhs_i
    rkUpdateKernel(Nelements,
                   rk,
                   _dt,
                   o_rkX,
                   o_rkA,
                   o_rkE,
                   o_q,
                   o_rhsq,
                   o_rkrhsq,
                   o_rkq,
                   o_rkerr);
  }
}

dfloat sark4::Estimater(deviceMemory<dfloat>& o_q){

  //Error estimation
  //E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
  //      DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
  rkErrorEstimateKernel(Nelements*Np*Nfields,
                        ATOL,
                        RTOL,
                        o_q,
                        o_rkq,
                        o_rkerr,
                        o_errtmp);

  h_errtmp.copyFrom(o_errtmp);
  dfloat err = 0;
  for(dlong n=0;n<Nblock;++n){
    err += h_errtmp[n];
  }
  comm.Allreduce(err);

  err = sqrt(err)*sqrtinvNtotal;

  return err;
}

static void UpdateCoefficients(const int Nfields,
                               const memory<dfloat> &lambda,
                               const dfloat dt,
                               pinnedMemory<dfloat> &h_rkX,
                               pinnedMemory<dfloat> &h_rkA,
                               pinnedMemory<dfloat> &h_rkE) {

  constexpr int Nrk = 5;
  const int Nr = 32;
  complex<dfloat> R[Nr];

  //contour integral integration points
  for(int ind=1; ind <= Nr; ++ind){
    const dfloat theta = (dfloat) (ind - 0.5) / (dfloat) Nr;
    complex<dfloat> z(0., M_PI* theta);
    R[ind-1] = exp(z);
  }

  for (int n=0;n<Nfields;n++) {

    dfloat* rkX = h_rkX.ptr() + n * Nrk;
    dfloat* rkA = h_rkA.ptr() + n * Nrk * Nrk;
    dfloat* rkE = h_rkE.ptr() + n * Nrk;

    if (lambda[n]==0) { //Zero exponential term, usual RK coefficients

      rkX[0] = 1.0; rkX[1] = 1.0; rkX[2] = 1.0; rkX[3] = 1.0; rkX[4] = 1.0;

      rkA[0+0*Nrk] =     0.0; rkA[1+0*Nrk] =     0.0; rkA[2+0*Nrk] =     0.0; rkA[3+0*Nrk] =     0.0; rkA[4+0*Nrk] = 0.0;
      rkA[0+1*Nrk] =     0.5; rkA[1+1*Nrk] =     0.0; rkA[2+1*Nrk] =     0.0; rkA[3+1*Nrk] =     0.0; rkA[4+1*Nrk] = 0.0;
      rkA[0+2*Nrk] =     0.0; rkA[1+2*Nrk] =     0.5; rkA[2+2*Nrk] =     0.0; rkA[3+2*Nrk] =     0.0; rkA[4+2*Nrk] = 0.0;
      rkA[0+3*Nrk] =     0.0; rkA[1+3*Nrk] =     0.0; rkA[2+3*Nrk] =     1.0; rkA[3+3*Nrk] =     0.0; rkA[4+3*Nrk] = 0.0;
      rkA[0+4*Nrk] = 1.0/6.0; rkA[1+4*Nrk] = 1.0/3.0; rkA[2+4*Nrk] = 1.0/3.0; rkA[3+4*Nrk] = 1.0/6.0; rkA[4+4*Nrk] = 0.0;

      rkE[0] = 0.0; rkE[1] = 0.0; rkE[2] = 0.0; rkE[3] = -1.0/6.0; rkE[4] = 1.0/6.0;

    } else {

      //Compute RK coefficients via contour integral
      dfloat alpha = lambda[n]*dt;

      // Initialize complex variables for contour integral
      complex<double> ca21(0.,0.);
      complex<double> ca31(0.,0.);
      complex<double> ca32(0.,0.);
      complex<double> ca41(0.,0.);
      complex<double> ca43(0.,0.);
      complex<double> ca51(0.,0.);
      complex<double> ca52(0.,0.);
      complex<double> ca53(0.,0.);
      complex<double> ca54(0.,0.);

      for(int i = 0; i<Nr; ++i ){
        complex<double> lr = alpha  + R[i];

        ca21 +=  (exp(lr/2.) - 1.)/lr;

        ca31 +=  (4. + exp(lr/2.)*(-4. + lr) + lr)/pow(lr,2);
        ca32 +=  (4.*exp(lr/2.) -2.*lr -4.)/pow(lr,2);

        ca41 +=  ((exp(lr) + 1.)*lr + 2. - 2.*exp(lr))/pow(lr,2);
        ca43 +=  (2.*exp(lr) - 2.*lr - 2.)/pow(lr,2);

        ca51 +=  (exp(lr)*pow(lr,2) + (- 3.*exp(lr) - 1.)*lr + 4.*exp(lr) - 4.)/pow(lr,3);
        ca52 +=  ((2.*exp(lr) + 2.)*lr + 4. - 4.*exp(lr))/pow(lr,3);
        ca53 +=  ((2.*exp(lr) + 2.)*lr + 4. - 4.*exp(lr))/pow(lr,3);
        ca54 +=  ((-exp(lr) - 3.)*lr - pow(lr,2) + 4.*exp(lr) - 4.)/pow(lr,3);
      }


      dfloat a21=real(ca21)/ (double) Nr;

      dfloat a31=real(ca31)/ (double) Nr;
      dfloat a32=real(ca32)/ (double) Nr;

      dfloat a41=real(ca41)/ (double) Nr;
      dfloat a43=real(ca43)/ (double) Nr;

      dfloat a51=real(ca51)/ (double) Nr;
      dfloat a52=real(ca52)/ (double) Nr;
      dfloat a53=real(ca53)/ (double) Nr;
      dfloat a54=real(ca54)/ (double) Nr;

      // first set non-semianalytic part of the integrator
      rkX[0] = 1.0; rkX[1] = std::exp(dfloat(0.5)*alpha); rkX[2] = std::exp(dfloat(0.5)*alpha); rkX[3] = std::exp(alpha); rkX[4] = std::exp(alpha);

      rkA[0+0*Nrk] = 0.0; rkA[1+0*Nrk] = 0.0; rkA[2+0*Nrk] = 0.0; rkA[3+0*Nrk] = 0.0; rkA[4+0*Nrk] = 0.0;
      rkA[0+1*Nrk] = a21; rkA[1+1*Nrk] = 0.0; rkA[2+1*Nrk] = 0.0; rkA[3+1*Nrk] = 0.0; rkA[4+1*Nrk] = 0.0;
      rkA[0+2*Nrk] = a31; rkA[1+2*Nrk] = a32; rkA[2+2*Nrk] = 0.0; rkA[3+2*Nrk] = 0.0; rkA[4+2*Nrk] = 0.0;
      rkA[0+3*Nrk] = a41; rkA[1+3*Nrk] = 0.0; rkA[2+3*Nrk] = a43; rkA[3+3*Nrk] = 0.0; rkA[4+3*Nrk] = 0.0;
      rkA[0+4*Nrk] = a51; rkA[1+4*Nrk] = a52; rkA[2+4*Nrk] = a53; rkA[3+4*Nrk] = a54; rkA[4+4*Nrk] = 0.0;

      rkE[0] = 0.0; rkE[1] = 0.0; rkE[2] = 0.0; rkE[3] = -a54; rkE[4] = a54;
    }
  }
}


/**************************************************/
/* PML version                                    */
/**************************************************/

sark4_pml::sark4_pml(dlong _Nelements, dlong NpmlElements, dlong _NhaloElements,
                     int _Np, int _Nfields, int _Npmlfields,
                     memory<dfloat> _lambda,
                     platform_t& _platform, comm_t _comm):
  pmlTimeStepperBase_t(_Nelements, NpmlElements, _NhaloElements,
                       _Np, _Nfields, _Npmlfields,
                       _platform, _comm),
  Np(_Np),
  Nfields(_Nfields),
  Nelements(_Nelements),
  NhaloElements(_NhaloElements) {

  lambda.malloc(Nfields);
  lambda.copyFrom(_lambda);

  //Nrk = 5;
  order = 4;
  embeddedOrder = 3;

  dlong Nlocal = Nelements*Np*Nfields;
  dlong Ntotal = (Nelements+NhaloElements)*Np*Nfields;

  o_rkq    = platform.malloc<dfloat>(Ntotal);
  o_rhsq   = platform.malloc<dfloat>(Nlocal);
  o_rkrhsq = platform.malloc<dfloat>(Nlocal*Nrk);
  o_rkpmlq    = platform.malloc<dfloat>(Npml);
  o_rhspmlq   = platform.malloc<dfloat>(Npml);
  o_rkrhspmlq = platform.malloc<dfloat>(Npml*Nrk);

  o_rkerr  = platform.malloc<dfloat>(Nlocal);

  o_saveq  = platform.malloc<dfloat>(Nlocal);
  o_savepmlq   = platform.malloc<dfloat>(Npml);

  const int blocksize=256;

  Nblock = (N+blocksize-1)/blocksize;
  h_errtmp = platform.hostMalloc<dfloat>(Nblock);
  o_errtmp = platform.malloc<dfloat>(Nblock);

  hlong gNtotal = Nlocal;
  comm.Allreduce(gNtotal);

  //copy base occa properties from platform
  properties_t kernelInfo = platform.props();

  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)blocksize;
  kernelInfo["defines/" "p_Nrk"]     = (int)Nrk;
  kernelInfo["defines/" "p_Np"]      = (int)Np;
  kernelInfo["defines/" "p_Nfields"] = (int)Nfields;

  rkUpdateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperSARK.okl",
                                    "sarkRkUpdate",
                                    kernelInfo);
  rkPmlUpdateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                      "timeStepperSARK.okl",
                                      "sarkRkPmlUpdate",
                                      kernelInfo);
  rkStageKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperSARK.okl",
                                    "sarkRkStage",
                                    kernelInfo);
  rkPmlStageKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                      "timeStepperSARK.okl",
                                      "sarkRkPmlStage",
                                      kernelInfo);
  rkErrorEstimateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperSARK.okl",
                                    "sarkErrorEstimate",
                                    kernelInfo);

  // Semi-Analytic Runge Kutta - order (3) 4 with PID timestep control
  dfloat _rkC[Nrk] = {0.0, 0.5, 0.5, 1.0, 1.0};
  rkC.malloc(Nrk);
  rkC.copyFrom(_rkC);

  h_rkX = platform.hostMalloc<dfloat>(Nfields*Nrk);
  h_rkA = platform.hostMalloc<dfloat>(Nfields*Nrk*Nrk);
  h_rkE = platform.hostMalloc<dfloat>(Nfields*Nrk);

  o_rkX = platform.malloc<dfloat>(Nfields*Nrk);
  o_rkA = platform.malloc<dfloat>(Nfields*Nrk*Nrk);
  o_rkE = platform.malloc<dfloat>(Nfields*Nrk);

  pmlrkA.malloc(Nrk*Nrk);

  dfloat _pmlrkA[Nrk*Nrk]
                  = {      0.0,      0.0,       0.0,      0.0,       0.0,
                           0.5,      0.0,       0.0,      0.0,       0.0,
                           0.0,      0.5,       0.0,      0.0,       0.0,
                           0.0,      0.0,       1.0,      0.0,       0.0,
                       1.0/6.0,  1.0/3.0,   1.0/3.0,   1.0/6.0,      0.0};
  pmlrkA.copyFrom(_pmlrkA);

  o_pmlrkA = platform.malloc<dfloat>(pmlrkA);

  dtMIN = 1E-9; //minumum allowed timestep
  ATOL = 1E-5;  //absolute error tolerance
  RTOL = 1E-4;  //relative error tolerance
  safe = 0.8;   //safety factor

  //error control parameters
  beta = 0.05;
  factor1 = 0.2;
  factor2 = 10.0;

  exp1 = 1.0/(embeddedOrder+1) - 0.75*beta;
  invfactor1 = 1.0/factor1;
  invfactor2 = 1.0/factor2;
  facold = 1E-4;
  sqrtinvNtotal = 1.0/sqrt(gNtotal);
}

void sark4_pml::Run(solver_t& solver,
                    deviceMemory<dfloat> &o_q,
                    deviceMemory<dfloat> &o_pmlq,
                    dfloat start, dfloat end) {

  dfloat time = start;

  solver.Report(time,0);

  dfloat outputInterval=0.0;
  solver.settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat outputTime = time + outputInterval;

  int tstep=0, allStep=0;

  //Compute Butcher Tableau
  UpdateCoefficients(Nfields, lambda, dt, h_rkX, h_rkA, h_rkE);

  // move data to platform
  h_rkX.copyTo(o_rkX, properties_t("async", true));
  h_rkA.copyTo(o_rkA, properties_t("async", true));
  h_rkE.copyTo(o_rkE, properties_t("async", true));

  while (time < end) {

    LIBP_ABORT("Time step became too small at time step = " << tstep,
               dt<dtMIN);
    LIBP_ABORT("Solution became unstable at time step = " << tstep,
               std::isnan(dt));

    //check for final timestep
    if (time+dt > end){
      dt = end-time;
    }

    Step(solver, o_q, o_pmlq, time, dt);

    // compute Dopri estimator
    dfloat err = Estimater(o_q);

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
        Backup(o_rkq, o_rkpmlq);

        // change dt to match output
        dt = outputTime-time;

        //Compute Butcher Tableau
        UpdateCoefficients(Nfields, lambda, dt, h_rkX, h_rkA, h_rkE);

        // move data to platform
        h_rkX.copyTo(o_rkX, properties_t("async", true));
        h_rkA.copyTo(o_rkA, properties_t("async", true));
        h_rkE.copyTo(o_rkE, properties_t("async", true));

        // time step to output
        Step(solver, o_q, o_pmlq, time, dt);

        // shift for output
        o_rkq.copyTo(o_q, properties_t("async", true));

        // output  (print from rkq)
        solver.Report(outputTime,tstep);

        // restore time step
        dt = savedt;

        // increment next output time
        outputTime += outputInterval;

        // accept saved rkq
        Restore(o_rkq, o_rkpmlq);
        AcceptStep(o_q, o_rkq, o_pmlq, o_rkpmlq);

      } else {
        // accept rkq
        AcceptStep(o_q, o_rkq, o_pmlq, o_rkpmlq);
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

    //Compute Butcher Tableau
    UpdateCoefficients(Nfields, lambda, dt, h_rkX, h_rkA, h_rkE);

    // move data to platform
    h_rkX.copyTo(o_rkX, properties_t("async", true));
    h_rkA.copyTo(o_rkA, properties_t("async", true));
    h_rkE.copyTo(o_rkE, properties_t("async", true));

    allStep++;
  }
}

void sark4_pml::Backup(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_pmlq) {
  o_saveq.copyFrom(o_q, N, properties_t("async", true));
  o_savepmlq.copyFrom(o_rkpmlq, Npml, properties_t("async", true));
}

void sark4_pml::Restore(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_pmlq) {
  o_saveq.copyTo(o_q, N, properties_t("async", true));
  o_savepmlq.copyTo(o_rkpmlq, Npml, properties_t("async", true));
}

void sark4_pml::AcceptStep(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rq,
                           deviceMemory<dfloat> &o_pmlq, deviceMemory<dfloat> &o_rpmlq) {
  o_q.copyFrom(o_rq, N, properties_t("async", true));
  o_pmlq.copyFrom(o_rpmlq, Npml, properties_t("async", true));
}

void sark4_pml::Step(solver_t& solver,
                     deviceMemory<dfloat> &o_q,
                     deviceMemory<dfloat> &o_pmlq,
                     dfloat time, dfloat _dt) {

  //RK step
  for(int rk=0;rk<Nrk;++rk){

    // t_rk = t + C_rk*_dt
    dfloat currentTime = time + rkC[rk]*_dt;

    //compute RK stage
    // rkq = x_{rk}*q + _dt sum_{i=0}^{rk-1} a_{rk,i}*rhsq_i
    rkStageKernel(Nelements,
                  rk,
                  _dt,
                  o_rkX,
                  o_rkA,
                  o_q,
                  o_rkrhsq,
                  o_rkq);
    if (Npml)
      rkPmlStageKernel(Npml,
                      rk,
                      _dt,
                      o_pmlrkA,
                      o_pmlq,
                      o_rkrhspmlq,
                      o_rkpmlq);

    //evaluate ODE rhs = f(q,t)
    solver.rhsf_pml(o_rkq, o_rkpmlq, o_rhsq, o_rhspmlq, currentTime);

    // update solution using Runge-Kutta
    // rkrhsq_rk = rhsq
    // if rk==6
    //   q = rkX_{rk}*q + _dt*sum_{i=0}^{rk} rkA_{rk,i}*rkrhs_i
    //   rkerr = _dt*sum_{i=0}^{rk} rkE_{rk,i}*rkrhs_i
    rkUpdateKernel(Nelements,
                   rk,
                   _dt,
                   o_rkX,
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
                         o_pmlrkA,
                         o_pmlq,
                         o_rhspmlq,
                         o_rkrhspmlq,
                         o_rkpmlq);
  }
}

dfloat sark4_pml::Estimater(deviceMemory<dfloat>& o_q){

  //Error estimation
  //E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
  //      DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
  rkErrorEstimateKernel(Nelements*Np*Nfields,
                        ATOL,
                        RTOL,
                        o_q,
                        o_rkq,
                        o_rkerr,
                        o_errtmp);

  h_errtmp.copyFrom(o_errtmp);
  dfloat err = 0;
  for(dlong n=0;n<Nblock;++n){
    err += h_errtmp[n];
  }
  comm.Allreduce(err);

  err = sqrt(err)*sqrtinvNtotal;

  return err;
}

} //namespace TimeStepper

} //namespace libp
