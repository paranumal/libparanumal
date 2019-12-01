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

/* Adams Bashforth, order 3 */
ab3::ab3(dlong Nelements, dlong NhaloElements,
         int Np, int Nfields, solver_t& _solver):
  timeStepper_t(Nelements, NhaloElements, Np, Nfields, _solver) {

  Nstages = 3;
  shiftIndex = 0;

  o_rhsq = device.malloc(Nstages*N*sizeof(dfloat));

  occa::properties kernelInfo = props; //copy base occa properties from solver

  kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;
  kernelInfo["defines/" "p_Nstages"] = Nstages;

  updateKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                    "timeStepperAB.okl",
                                    "abUpdate",
                                    kernelInfo, comm);

  // initialize AB time stepping coefficients
  dfloat _ab_a[Nstages*Nstages] = {
                           1.0,      0.0,    0.0,
                         3./2.,   -1./2.,    0.0,
                       23./12., -16./12., 5./12.};

  ab_a = (dfloat*) calloc(Nstages*Nstages, sizeof(dfloat));
  memcpy(ab_a, _ab_a, Nstages*Nstages*sizeof(dfloat));

  o_ab_a = device.malloc(Nstages*Nstages*sizeof(dfloat), ab_a);
}

void ab3::Run(occa::memory &o_q, dfloat start, dfloat end) {

  dfloat time = start;

  solver.Report(time,0);

  dfloat outputInterval;
  settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat outputTime = time + outputInterval;

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

void ab3::Step(occa::memory &o_q, dfloat time, dfloat _dt, int order) {

  //rhs at current index
  occa::memory o_rhsq0 = o_rhsq + shiftIndex*N*sizeof(dfloat);

  //A coefficients at current order
  occa::memory o_A = o_ab_a + order*Nstages*sizeof(dfloat);

  //evaluate ODE rhs = f(q,t)
  solver.rhsf(o_q, o_rhsq0, time);

  //update q
  updateKernel(N,
               dt,
               shiftIndex,
               o_A,
               o_rhsq,
               o_q);

  //rotate index
  shiftIndex = (shiftIndex+Nstages-1)%Nstages;
}

ab3::~ab3() {
  if (o_rhsq.size()) o_rhsq.free();
  if (o_ab_a.size()) o_ab_a.free();

  if (ab_a) free(ab_a);

  updateKernel.free();
}

/**************************************************/
/* PML version                                    */
/**************************************************/

/* Adams Bashforth, order 3 */
ab3_pml::ab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
                int Np, int Nfields, int Npmlfields, solver_t& _solver):
  ab3(Nelements, NhaloElements, Np, Nfields, _solver),
  Npml(NpmlElements*Np*Npmlfields) {

  if (Npml) {
    dfloat *pmlq = (dfloat *) calloc(Npml,sizeof(dfloat));
    o_pmlq   = device.malloc(Npml*sizeof(dfloat), pmlq);
    free(pmlq);

    o_rhspmlq = device.malloc(Nstages*Npml*sizeof(dfloat));
  }
}

void ab3_pml::Step(occa::memory &o_q, dfloat time, dfloat _dt, int order) {

  //rhs at current index
  occa::memory o_rhsq0    = o_rhsq    + shiftIndex*N*sizeof(dfloat);
  occa::memory o_rhspmlq0;
  if (Npml)    o_rhspmlq0 = o_rhspmlq + shiftIndex*Npml*sizeof(dfloat);

  //A coefficients at current order
  occa::memory o_A = o_ab_a + order*Nstages*sizeof(dfloat);

  //evaluate ODE rhs = f(q,t)
  solver.rhsf_pml(o_q, o_pmlq, o_rhsq0, o_rhspmlq0, time);

  //update q & pmlq
  updateKernel(N,
               dt,
               shiftIndex,
               o_A,
               o_rhsq,
               o_q);
  if (Npml)
    updateKernel(Npml,
                 dt,
                 shiftIndex,
                 o_A,
                 o_rhspmlq,
                 o_pmlq);

  //rotate index
  shiftIndex = (shiftIndex+Nstages-1)%Nstages;
}

ab3_pml::~ab3_pml() {
  if (o_pmlq.size()) o_pmlq.free();
  if (o_rhspmlq.size()) o_rhspmlq.free();
}

} //namespace TimeStepper
