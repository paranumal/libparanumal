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

#include <math.h>
#include "core.hpp"
#include "timeStepper.hpp"
#include <complex>

namespace TimeStepper {

using std::complex;

sark4::sark4(dlong _Nelements, dlong _NhaloElements,
             int _Np, int _Nfields,
             dfloat *_lambda, solver_t& _solver):
  timeStepper_t(_Nelements, _NhaloElements, _Np, _Nfields, _solver),
  Np(_Np),
  Nfields(_Nfields),
  Nelements(_Nelements),
  NhaloElements(_NhaloElements) {

  lambda = (dfloat *) malloc(Nfields*sizeof(dfloat));
  memcpy(lambda, _lambda, Nfields*sizeof(dfloat));

  Nrk = 5;
  order = 4;
  embeddedOrder = 3;

  dlong Nlocal = Nelements*Np*Nfields;
  dlong Ntotal = (Nelements+NhaloElements)*Np*Nfields;

  o_rkq    = device.malloc(Ntotal*sizeof(dfloat));
  o_rhsq   = device.malloc(Nlocal*sizeof(dfloat));
  o_rkrhsq = device.malloc(Nlocal*Nrk*sizeof(dfloat));
  o_rkerr  = device.malloc(Nlocal*sizeof(dfloat));

  o_saveq  = device.malloc(Nlocal*sizeof(dfloat));

  Nblock = (N+BLOCKSIZE-1)/BLOCKSIZE;
  errtmp = (dfloat*) calloc(Nblock, sizeof(dfloat));
  o_errtmp = device.malloc(Nblock*sizeof(dfloat));

  hlong gNlocal = Nlocal;
  hlong gNtotal;
  MPI_Allreduce(&gNlocal, &gNtotal, 1, MPI_HLONG, MPI_SUM, comm);

  occa::properties kernelInfo = props; //copy base occa properties from solver

  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)BLOCKSIZE;
  kernelInfo["defines/" "p_Nrk"]     = (int)Nrk;
  kernelInfo["defines/" "p_Np"]      = (int)Np;
  kernelInfo["defines/" "p_Nfields"] = (int)Nfields;

  rkUpdateKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                    "timeStepperSARK.okl",
                                    "sarkRkUpdate",
                                    kernelInfo, comm);

  rkStageKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                    "timeStepperSARK.okl",
                                    "sarkRkStage",
                                    kernelInfo, comm);

  rkErrorEstimateKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                    "timeStepperSARK.okl",
                                    "sarkErrorEstimate",
                                    kernelInfo, comm);

  // Semi-Analytic Runge Kutta - order (3) 4 with PID timestep control
  dfloat _rkC[Nrk] = {0.0, 0.5, 0.5, 1.0, 1.0};
  rkC = (dfloat*) calloc(Nrk, sizeof(dfloat));
  memcpy(rkC, _rkC, Nrk*sizeof(dfloat));

  rkX = (dfloat*) occaHostMallocPinned(device, Nfields*Nrk*    sizeof(dfloat),
                                       NULL, o_rkX, h_rkX);
  rkA = (dfloat*) occaHostMallocPinned(device, Nfields*Nrk*Nrk*sizeof(dfloat),
                                       NULL, o_rkA, h_rkA);
  rkE = (dfloat*) occaHostMallocPinned(device, Nfields*Nrk*    sizeof(dfloat),
                                       NULL, o_rkE, h_rkE);

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

void sark4::Run(occa::memory &o_q, dfloat start, dfloat end) {

  dfloat time = start;

  int rank;
  MPI_Comm_rank(comm, &rank);

  solver.Report(time,0);

  dfloat outputInterval;
  settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat outputTime = time + outputInterval;

  int tstep=0, allStep=0;

  //Compute Butcher Tableau
  UpdateCoefficients();

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

        //Compute Butcher Tableau
        UpdateCoefficients();

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

    //Compute Butcher Tableau
    UpdateCoefficients();

    allStep++;
  }

  if (!rank)
    printf("%d accepted steps and %d total steps\n", tstep, allStep);
}

void sark4::Backup(occa::memory &o_Q) {
  o_saveq.copyFrom(o_Q, N*sizeof(dfloat));
}

void sark4::Restore(occa::memory &o_Q) {
  o_saveq.copyTo(o_Q, N*sizeof(dfloat));
}

void sark4::AcceptStep(occa::memory &o_q, occa::memory &o_rq) {
  o_q.copyFrom(o_rq, N*sizeof(dfloat));
}

void sark4::Step(occa::memory &o_q, dfloat time, dfloat _dt) {

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

dfloat sark4::Estimater(occa::memory& o_q){

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

void sark4::UpdateCoefficients() {

  const int Nr = 32;
  complex<dfloat> R[Nr];

  //contour integral integration points
  for(int ind=1; ind <= Nr; ++ind){
    const dfloat theta = (dfloat) (ind - 0.5) / (dfloat) Nr;
    complex<dfloat> z(0., M_PI* theta);
    R[ind-1] = exp(z);
  }

  for (int n=0;n<Nfields;n++) {

    if (lambda[n]==0) { //Zero exponential term, usual RK coefficients

      dfloat _rkX[Nrk]  = {1.0, 1.0, 1.0, 1.0, 1.0 };
      dfloat _rkA[Nrk*Nrk]
                    = {      0.0,      0.0,       0.0,      0.0,       0.0,
                             0.5,      0.0,       0.0,      0.0,       0.0,
                             0.0,      0.5,       0.0,      0.0,       0.0,
                             0.0,      0.0,       1.0,      0.0,       0.0,
                         1.0/6.0,  1.0/3.0,   1.0/3.0,   1.0/6.0,      0.0};
      dfloat _rkE[Nrk]= {    0.0,      0.0,       0.0,  -1.0/6.0,  1.0/6.0};

      memcpy(rkX+n*Nrk    ,_rkX,    Nrk*sizeof(dfloat));
      memcpy(rkA+n*Nrk*Nrk,_rkA,Nrk*Nrk*sizeof(dfloat));
      memcpy(rkE+n*Nrk    ,_rkE,    Nrk*sizeof(dfloat));

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
      dfloat _rkX[Nrk]  = {1.0, exp(0.5*alpha), exp(0.5*alpha), exp(alpha), exp(alpha) };
      dfloat _rkA[Nrk*Nrk]
                      ={   0.0,  0.0,  0.0,   0.0,  0.0,
                           a21,  0.0,  0.0,   0.0,  0.0,
                           a31,  a32,  0.0,   0.0,  0.0,
                           a41,  0.0,  a43,   0.0,  0.0,
                           a51,  a52,  a53,   a54,  0.0};
      dfloat _rkE[Nrk]= {  0.0,  0.0,  0.0,  -a54,  a54};

      memcpy(rkX+n*Nrk    ,_rkX,    Nrk*sizeof(dfloat));
      memcpy(rkA+n*Nrk*Nrk,_rkA,Nrk*Nrk*sizeof(dfloat));
      memcpy(rkE+n*Nrk    ,_rkE,    Nrk*sizeof(dfloat));
    }

    // move data to device
    // o_rkX.copyFrom(rkX, "async: true");
    // o_rkA.copyFrom(rkA, "async: true");
    // o_rkE.copyFrom(rkE, "async: true");
    o_rkX.copyFrom(rkX);
    o_rkA.copyFrom(rkA);
    o_rkE.copyFrom(rkE);
  }
}

sark4::~sark4() {
  if (o_rkq.size()) o_rkq.free();
  if (o_rkrhsq.size()) o_rkrhsq.free();
  if (o_rkerr.size()) o_rkerr.free();
  if (o_errtmp.size()) o_errtmp.free();
  if (o_rkX.size()) o_rkX.free();
  if (o_rkA.size()) o_rkA.free();
  if (o_rkE.size()) o_rkE.free();

  if (errtmp) free(errtmp);
  if (lambda) free(lambda);
  if (rkC) free(rkC);

  if (h_rkX.size()) h_rkX.free();
  if (h_rkA.size()) h_rkA.free();
  if (h_rkE.size()) h_rkE.free();

  rkUpdateKernel.free();
  rkStageKernel.free();
  rkErrorEstimateKernel.free();
}

/**************************************************/
/* PML version                                    */
/**************************************************/

sark4_pml::sark4_pml(dlong _Nelements, dlong _NpmlElements, dlong _NhaloElements,
            int _Np, int _Nfields, int _Npmlfields,
            dfloat *_lambda, solver_t& _solver):
  sark4(_Nelements, _NhaloElements, _Np, _Nfields, _lambda, _solver),
  Npml(_Npmlfields*_Np*_NpmlElements) {

  if (Npml) {
    dfloat *pmlq = (dfloat *) calloc(Npml,sizeof(dfloat));
    o_pmlq   = device.malloc(Npml*sizeof(dfloat), pmlq);
    free(pmlq);

    o_rkpmlq    = device.malloc(Npml*sizeof(dfloat));
    o_rhspmlq   = device.malloc(Npml*sizeof(dfloat));
    o_rkrhspmlq = device.malloc(Npml*Nrk*sizeof(dfloat));

    o_savepmlq   = device.malloc(Npml*sizeof(dfloat));

    occa::properties kernelInfo = props; //copy base occa properties from solver

    //add defines
    kernelInfo["defines/" "p_blockSize"] = (int)BLOCKSIZE;
    kernelInfo["defines/" "p_Nrk"]     = (int)Nrk;
    kernelInfo["defines/" "p_Np"]      = (int)Np;
    kernelInfo["defines/" "p_Nfields"] = (int)Nfields;

    rkPmlUpdateKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                      "timeStepperSARK.okl",
                                      "sarkRkPmlUpdate",
                                      kernelInfo, comm);

    rkPmlStageKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                      "timeStepperSARK.okl",
                                      "sarkRkPmlStage",
                                      kernelInfo, comm);

    // Semi-Analytic Runge Kutta - order (3) 4 with PID timestep control
    pmlrkA = (dfloat*) malloc(Nrk*Nrk*sizeof(dfloat));

    dfloat _pmlrkA[Nrk*Nrk]
                    = {      0.0,      0.0,       0.0,      0.0,       0.0,
                             0.5,      0.0,       0.0,      0.0,       0.0,
                             0.0,      0.5,       0.0,      0.0,       0.0,
                             0.0,      0.0,       1.0,      0.0,       0.0,
                         1.0/6.0,  1.0/3.0,   1.0/3.0,   1.0/6.0,      0.0};
    memcpy(pmlrkA, _pmlrkA, Nrk*Nrk*sizeof(dfloat));

    o_pmlrkA = device.malloc(Nrk*Nrk*sizeof(dfloat), pmlrkA);
  }
}

void sark4_pml::Backup(occa::memory &o_Q) {
  o_saveq.copyFrom(o_Q, N*sizeof(dfloat));
  if (Npml)
    o_savepmlq.copyFrom(o_rkpmlq, Npml*sizeof(dfloat));
}

void sark4_pml::Restore(occa::memory &o_Q) {
  o_saveq.copyTo(o_Q, N*sizeof(dfloat));
  if (Npml)
    o_savepmlq.copyTo(o_rkpmlq, Npml*sizeof(dfloat));
}

void sark4_pml::AcceptStep(occa::memory &o_q, occa::memory &o_rq) {
  o_q.copyFrom(o_rq, N*sizeof(dfloat));
  if (Npml)
    o_pmlq.copyFrom(o_rkpmlq, Npml*sizeof(dfloat));
}

void sark4_pml::Step(occa::memory &o_q, dfloat time, dfloat _dt) {

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

sark4_pml::~sark4_pml() {
  if (o_pmlq.size()) o_pmlq.free();
  if (o_rkpmlq.size()) o_rkpmlq.free();
  if (o_rhspmlq.size()) o_rhspmlq.free();
  if (o_rkrhspmlq.size()) o_rkrhspmlq.free();
  if (o_pmlrkA.size()) o_pmlrkA.free();

  if (o_savepmlq.size()) o_savepmlq.free();

  rkPmlUpdateKernel.free();
  rkPmlStageKernel.free();
}

} //namespace TimeStepper
