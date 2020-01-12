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

/* Backward Difference Formula, order 3, with subcycling */
ssbdf3::ssbdf3(dlong Nelements, dlong NhaloElements,
                 int Np, int Nfields, solver_t& _solver):
  timeStepper_t(Nelements, NhaloElements, Np, Nfields, _solver) {

  Nstages = 3;
  shiftIndex = 0;

  o_qn   = device.malloc(Nstages*N*sizeof(dfloat)); //q history
  o_qhat = device.malloc(Nstages*N*sizeof(dfloat)); //F(q) history (explicit part)
  o_rhs  = device.malloc(N*sizeof(dfloat)); //rhs storage

  occa::properties kernelInfo = props; //copy base occa properties from solver

  kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;
  kernelInfo["defines/" "p_Nstages"] = Nstages;

  rhsKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                    "timeStepperSSBDF.okl",
                                    "ssbdfRHS",
                                    kernelInfo, comm);

  // initialize BDF time stepping coefficients
  dfloat _b[Nstages*(Nstages+1)] = {
                            1.,    1.,     0.,    0.,
                         3./2.,    2., -1./2.,    0.,
                        11./6.,    3., -3./2., 1./3.};

  ssbdf_b = (dfloat*) calloc(Nstages*(Nstages+1), sizeof(dfloat));
  memcpy(ssbdf_b, _b, Nstages*(Nstages+1)*sizeof(dfloat));

  o_ssbdf_b = device.malloc(Nstages*(Nstages+1)*sizeof(dfloat), ssbdf_b);
}

dfloat ssbdf3::getGamma() {
  return *(ssbdf_b + (Nstages-1)*(Nstages+1)); //first entry of last row of B
}

void ssbdf3::Run(occa::memory &o_q, dfloat start, dfloat end) {

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

void ssbdf3::Step(occa::memory &o_q, dfloat time, dfloat _dt, int order) {

  //BDF coefficients at current order
  occa::memory o_B = o_ssbdf_b + order*(Nstages+1)*sizeof(dfloat);
  dfloat *B = ssbdf_b + order*(Nstages+1);

  //put current q into history
  occa::memory o_qn0 = o_qn + shiftIndex*N*sizeof(dfloat);
  o_qn0.copyFrom(o_q, N*sizeof(dfloat));

  // Compute qhat = sum_i=1^s B_i qhat(t_n+1-i) by
  // where qhat(t) is the Lagrangian state of q
  // by subcycling each of the history states q(t_n-i)
  solver.rhs_subcycle_f(o_qn, o_qhat, time, _dt, B, order, shiftIndex, Nstages);

  //build rhs for implicit step and update history
  rhsKernel(N,
           _dt,
           shiftIndex,
           o_B,
           o_qhat,
           o_rhs);

  dfloat gamma = B[0]/_dt;

  //solve implicit part:
  // find q such that gamma*q - G(q) = rhs
  int iter = solver.rhs_imex_invg(o_rhs, o_q, gamma, time+_dt);

  //rotate index
  shiftIndex = (shiftIndex+Nstages-1)%Nstages;
}

ssbdf3::~ssbdf3() {
  if (o_rhs.size()) o_rhs.free();
  if (o_qn.size()) o_qn.free();
  if (o_qhat.size()) o_qhat.free();
  if (o_ssbdf_b.size()) o_ssbdf_b.free();

  if (ssbdf_b) free(ssbdf_b);

  rhsKernel.free();
}

} //namespace TimeStepper
