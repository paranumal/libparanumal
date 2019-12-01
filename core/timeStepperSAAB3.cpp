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
#include <complex>

namespace TimeStepper {

using std::complex;

saab3::saab3(dlong _Nelements, dlong _NhaloElements,
             int _Np, int _Nfields,
             dfloat *_lambda, solver_t& _solver):
  timeStepper_t(_Nelements, _NhaloElements, _Np, _Nfields, _solver),
  Np(_Np),
  Nfields(_Nfields),
  Nelements(_Nelements),
  NhaloElements(_NhaloElements) {

  lambda = (dfloat *) malloc(Nfields*sizeof(dfloat));
  memcpy(lambda, _lambda, Nfields*sizeof(dfloat));

  Nstages = 3;
  shiftIndex = 0;

  o_rhsq = device.malloc(Nstages*N*sizeof(dfloat));

  occa::properties kernelInfo = props; //copy base occa properties from solver

  kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;
  kernelInfo["defines/" "p_Nstages"] = Nstages;
  kernelInfo["defines/" "p_Np"]      = (int)Np;
  kernelInfo["defines/" "p_Nfields"] = (int)Nfields;

  updateKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                    "timeStepperSAAB.okl",
                                    "saabUpdate",
                                    kernelInfo, comm);

  saab_x = (dfloat*) malloc(Nfields*sizeof(dfloat));
  o_saab_x = device.malloc(Nfields*sizeof(dfloat));

  saab_a = (dfloat*) malloc(Nfields*Nstages*Nstages*sizeof(dfloat));
  o_saab_a =  device.malloc(Nfields*Nstages*Nstages*sizeof(dfloat));
}

void saab3::Run(occa::memory &o_q, dfloat start, dfloat end) {

  dfloat time = start;

  solver.Report(time,0);

  dfloat outputInterval;
  settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat outputTime = time + outputInterval;

  //Compute SAAB coefficients
  UpdateCoefficients();

  int tstep=0;
  int order=0;
  while (time < end) {
    Step(o_q, time, dt, order);
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

void saab3::Step(occa::memory &o_q, dfloat time, dfloat _dt, int order) {

  //rhs at current index
  occa::memory o_rhsq0 = o_rhsq + shiftIndex*N*sizeof(dfloat);

  //coefficients at current order
  occa::memory o_X = o_saab_x;
  occa::memory o_A = o_saab_a + order*Nstages*sizeof(dfloat);

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

      memcpy(saab_x+n                ,_saab_X,    1*sizeof(dfloat));
      memcpy(saab_a+n*Nstages*Nstages,_saab_A,Nstages*Nstages*sizeof(dfloat));

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

      dfloat _saab_X[1]  = { exp(alpha) };
      dfloat _saab_A[Nstages*Nstages]
                      ={   aa11,   0.0,   0.0,
                           aa21,  aa22,   0.0,
                           aa31,  aa32,  aa33 };

      memcpy(saab_x+n                ,_saab_X,    1*sizeof(dfloat));
      memcpy(saab_a+n*Nstages*Nstages,_saab_A,Nstages*Nstages*sizeof(dfloat));
    }

    // move data to device
    o_saab_x.copyFrom(saab_x);
    o_saab_a.copyFrom(saab_a);
  }
}

saab3::~saab3() {
  if (o_rhsq.size()) o_rhsq.free();
  if (o_saab_x.size()) o_saab_x.free();
  if (o_saab_a.size()) o_saab_a.free();

  updateKernel.free();
}


/**************************************************/
/* PML version                                    */
/**************************************************/

saab3_pml::saab3_pml(dlong _Nelements, dlong _NpmlElements, dlong _NhaloElements,
                     int _Np, int _Nfields, int Npmlfields,
                     dfloat *_lambda, solver_t& _solver):
  saab3(_Nelements, _NhaloElements, _Np, _Nfields, _lambda, _solver),
  Npml(Npmlfields*_Np*_NpmlElements) {

  if (Npml) {
    dfloat *pmlq = (dfloat *) calloc(Npml,sizeof(dfloat));
    o_pmlq   = device.malloc(Npml*sizeof(dfloat), pmlq);
    free(pmlq);

    o_rhspmlq = device.malloc(Nstages*Npml*sizeof(dfloat));

    occa::properties kernelInfo = props; //copy base occa properties from solver

    kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;
    kernelInfo["defines/" "p_Nstages"] = Nstages;
    kernelInfo["defines/" "p_Np"]      = (int)Np;
    kernelInfo["defines/" "p_Nfields"] = (int)Nfields;

    pmlUpdateKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                      "timeStepperSAAB.okl",
                                      "saabPmlUpdate",
                                      kernelInfo, comm);

    // initialize AB time stepping coefficients
    dfloat _ab_a[Nstages*Nstages] = {
                             1.0,      0.0,    0.0,
                           3./2.,   -1./2.,    0.0,
                         23./12., -16./12., 5./12.};

    pmlsaab_a = (dfloat*) calloc(Nstages*Nstages, sizeof(dfloat));
    memcpy(pmlsaab_a, _ab_a, Nstages*Nstages*sizeof(dfloat));

    o_pmlsaab_a = device.malloc(Nstages*Nstages*sizeof(dfloat), pmlsaab_a);
  }
}


void saab3_pml::Step(occa::memory &o_q, dfloat time, dfloat _dt, int order) {

  //rhs at current index
  occa::memory o_rhsq0    = o_rhsq + shiftIndex*N*sizeof(dfloat);
  occa::memory o_rhspmlq0;


  //coefficients at current order
  occa::memory o_X = o_saab_x;
  occa::memory o_A = o_saab_a + order*Nstages*sizeof(dfloat);
  occa::memory o_pmlA;

  if (Npml) {
    o_rhspmlq0 = o_rhspmlq + shiftIndex*Npml*sizeof(dfloat);
    o_pmlA = o_pmlsaab_a + order*Nstages*sizeof(dfloat);
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

saab3_pml::~saab3_pml() {
  if (o_pmlq.size()) o_pmlq.free();
  if (o_rhspmlq.size()) o_rhspmlq.free();
  if (o_pmlsaab_a.size()) o_pmlsaab_a.free();

  pmlUpdateKernel.free();
}

} //namespace TimeStepper
