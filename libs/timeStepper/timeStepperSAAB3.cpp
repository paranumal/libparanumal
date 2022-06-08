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
#include <complex>

namespace libp {

namespace TimeStepper {

using std::complex;

saab3::saab3(dlong _Nelements, dlong _NhaloElements,
             int _Np, int _Nfields,
             memory<dfloat> _lambda,
             platform_t& _platform, comm_t _comm):
  timeStepperBase_t(_Nelements, _NhaloElements,
                    _Np, _Nfields, _platform, _comm),
  Np(_Np),
  Nfields(_Nfields),
  Nelements(_Nelements),
  NhaloElements(_NhaloElements) {

  lambda.malloc(Nfields);
  lambda.copyFrom(_lambda);

  Nstages = 3;
  shiftIndex = 0;

  o_rhsq = platform.malloc<dfloat>(Nstages*N);

  const int blocksize=256;

  properties_t kernelInfo = platform.props(); //copy base occa properties from solver

  kernelInfo["defines/" "p_blockSize"] = blocksize;
  kernelInfo["defines/" "p_Nstages"] = Nstages;
  kernelInfo["defines/" "p_Np"]      = (int)Np;
  kernelInfo["defines/" "p_Nfields"] = (int)Nfields;

  updateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperSAAB.okl",
                                    "saabUpdate",
                                    kernelInfo);

  h_saab_x = platform.hostMalloc<dfloat>(Nfields);
  o_saab_x = platform.malloc<dfloat>(Nfields);

  h_saab_a = platform.hostMalloc<dfloat>(Nfields*Nstages*Nstages);
  o_saab_a = platform.malloc<dfloat>(Nfields*Nstages*Nstages);
}

void saab3::Run(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat start, dfloat end) {

  dfloat time = start;

  solver.Report(time,0);

  dfloat outputInterval=0.0;
  solver.settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat outputTime = time + outputInterval;

  //Compute SAAB coefficients
  UpdateCoefficients();

  int tstep=0;
  int order=0;
  while (time < end) {
    Step(solver, o_q, time, dt, order);
    time += dt;
    tstep++;
    if (order<Nstages-1) order++;

    if (time>outputTime) {
      //report state
      solver.Report(time,tstep);
      outputTime += outputInterval;
    }
  }
}

void saab3::Step(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat time, dfloat _dt, int order) {

  //rhs at current index
  deviceMemory<dfloat> o_rhsq0 = o_rhsq + shiftIndex*N;

  //coefficients at current order
  deviceMemory<dfloat> o_X = o_saab_x;
  deviceMemory<dfloat> o_A = o_saab_a + order*Nstages;

  //evaluate ODE rhs = f(q,t)
  solver.rhsf(o_q, o_rhsq0, time);

  //update q
  updateKernel(Nelements,
               dt,
               shiftIndex,
               o_X,
               o_A,
               o_rhsq,
               o_q);

  //rotate index
  shiftIndex = (shiftIndex+Nstages-1)%Nstages;
}

void saab3::UpdateCoefficients() {

  const int Nr = 32;
  complex<dfloat> R[Nr];

  //contour integral integration points
  for(int ind=1; ind <= Nr; ++ind){
    const dfloat theta = (dfloat) (ind - 0.5) / (dfloat) Nr;
    complex<dfloat> z(0., M_PI* theta);
    R[ind-1] = exp(z);
  }

  for (int n=0;n<Nfields;n++) {

    if (lambda[n]==0) { //Zero exponential term, usual AB coefficients

      dfloat _saab_X[1]  = { 1.0 };
      dfloat _saab_A[Nstages*Nstages]
                    = {    1.0,      0.0,    0.0,
                         3./2.,   -1./2.,    0.0,
                       23./12., -16./12., 5./12.};

      h_saab_x.copyFrom(_saab_X,               1, n                );
      h_saab_a.copyFrom(_saab_A, Nstages*Nstages, n*Nstages*Nstages);

    } else {

      //Compute coefficients via contour integral
      dfloat alpha = lambda[n]*dt;

      // Initialize complex variables for contour integral
      complex<double> a11(0.,0.);
      complex<double> a21(0.,0.);
      complex<double> a22(0.,0.);
      complex<double> a31(0.,0.);
      complex<double> a32(0.,0.);
      complex<double> a33(0.,0.);

      for(int i = 0; i<Nr; ++i ){
        complex<double> lr = alpha  + R[i];

        a11 +=  (exp(lr) - 1.)/lr;

        a21 +=  (-2.*lr + (1.+lr)*exp(lr) - 1.)/pow(lr,2);
        a22 +=  (lr - exp(lr) + 1.)/pow(lr,2);

        a31 += (-2.5*lr - 3.*pow(lr,2) + (1.+pow(lr,2)+1.5*lr)*exp(lr) - 1.)/pow(lr,3);
        a32 += (4.*lr + 3.*pow(lr,2)- (2.*lr + 2.0)*exp(lr) + 2.)/pow(lr,3);
        a33 +=-(1.5*lr + pow(lr,2)- (0.5*lr + 1.)*exp(lr) + 1.)/pow(lr,3);
      }


      dfloat aa11=real(a11)/ (double) Nr;

      dfloat aa21=real(a21)/ (double) Nr;
      dfloat aa22=real(a22)/ (double) Nr;

      dfloat aa31=real(a31)/ (double) Nr;
      dfloat aa32=real(a32)/ (double) Nr;
      dfloat aa33=real(a33)/ (double) Nr;

      dfloat _saab_X[1]  = { std::exp(alpha) };
      dfloat _saab_A[Nstages*Nstages]
                      ={   aa11,   0.0,   0.0,
                           aa21,  aa22,   0.0,
                           aa31,  aa32,  aa33 };

      h_saab_x.copyFrom(_saab_X,               1, n                );
      h_saab_a.copyFrom(_saab_A, Nstages*Nstages, n*Nstages*Nstages);
    }

    // move data to platform
    h_saab_x.copyTo(o_saab_x);
    h_saab_a.copyTo(o_saab_a);
  }
}


/**************************************************/
/* PML version                                    */
/**************************************************/

saab3_pml::saab3_pml(dlong _Nelements, dlong _NpmlElements, dlong _NhaloElements,
                     int _Np, int _Nfields, int Npmlfields,
                     memory<dfloat> _lambda,
                     platform_t& _platform, comm_t _comm):
  saab3(_Nelements, _NhaloElements, _Np, _Nfields, _lambda, _platform, _comm),
  Npml(Npmlfields*_Np*_NpmlElements) {

  if (Npml) {
    memory<dfloat> pmlq(Npml,0.0);
    o_pmlq   = platform.malloc<dfloat>(pmlq);

    o_rhspmlq = platform.malloc<dfloat>(Nstages*Npml);

    properties_t kernelInfo = platform.props(); //copy base occa properties from solver

    const int blocksize=256;

    kernelInfo["defines/" "p_blockSize"] = blocksize;
    kernelInfo["defines/" "p_Nstages"] = Nstages;
    kernelInfo["defines/" "p_Np"]      = (int)Np;
    kernelInfo["defines/" "p_Nfields"] = (int)Nfields;

    pmlUpdateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                      "timeStepperSAAB.okl",
                                      "saabPmlUpdate",
                                      kernelInfo);

    // initialize AB time stepping coefficients
    dfloat _ab_a[Nstages*Nstages] = {
                             1.0,      0.0,    0.0,
                           3./2.,   -1./2.,    0.0,
                         23./12., -16./12., 5./12.};

    pmlsaab_a.malloc(Nstages*Nstages);
    pmlsaab_a.copyFrom(_ab_a);

    o_pmlsaab_a = platform.malloc<dfloat>(pmlsaab_a);
  }
}


void saab3_pml::Step(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat time, dfloat _dt, int order) {

  //rhs at current index
  deviceMemory<dfloat> o_rhsq0    = o_rhsq + shiftIndex*N;
  deviceMemory<dfloat> o_rhspmlq0;


  //coefficients at current order
  deviceMemory<dfloat> o_X = o_saab_x;
  deviceMemory<dfloat> o_A = o_saab_a + order*Nstages;
  deviceMemory<dfloat> o_pmlA;

  if (Npml) {
    o_rhspmlq0 = o_rhspmlq + shiftIndex*Npml;
    o_pmlA = o_pmlsaab_a + order*Nstages;
  }

  //evaluate ODE rhs = f(q,t)
  solver.rhsf_pml(o_q, o_pmlq, o_rhsq0, o_rhspmlq0, time);

  //update q
  updateKernel(Nelements,
               dt,
               shiftIndex,
               o_X,
               o_A,
               o_rhsq,
               o_q);
  //update pmlq
  if (Npml)
    pmlUpdateKernel(Npml,
                   dt,
                   shiftIndex,
                   o_pmlA,
                   o_rhspmlq,
                   o_pmlq);

  //rotate index
  shiftIndex = (shiftIndex+Nstages-1)%Nstages;
}

} //namespace TimeStepper

} //namespace libp
