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

/* Backward Difference Formula, order 3, with extrapolation */
extbdf3::extbdf3(dlong Nelements, dlong NhaloElements,
                 int Np, int Nfields, solver_t& _solver):
  timeStepper_t(Nelements, NhaloElements, Np, Nfields, _solver) {

  Nstages = 3;
  shiftIndex = 0;

  o_qn = device.malloc(Nstages*N*sizeof(dfloat)); //q history
  o_rhs = device.malloc(N*sizeof(dfloat)); //rhs storage

  dfloat *F = (dfloat *) calloc(Nstages*N,sizeof(dfloat));
  o_F  = device.malloc(Nstages*N*sizeof(dfloat), F); //F(q) history (explicit part)
  free(F);

  occa::properties kernelInfo = props; //copy base occa properties from solver

  kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;
  kernelInfo["defines/" "p_Nstages"] = Nstages;

  rhsKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                    "timeStepperEXTBDF.okl",
                                    "extbdfRHS",
                                    kernelInfo, comm);

  // initialize EXT and BDF time stepping coefficients
  dfloat _a[Nstages*Nstages] = {
                            1.,    0.,     0.,
                            2.,   -1.,     0.,
                            3.,   -3.,     1.};
  dfloat _b[Nstages*(Nstages+1)] = {
                            1.,    1.,     0.,    0.,
                         3./2.,    2., -1./2.,    0.,
                        11./6.,    3., -3./2., 1./3.};

  extbdf_a = (dfloat*) calloc(Nstages*Nstages, sizeof(dfloat));
  extbdf_b = (dfloat*) calloc(Nstages*(Nstages+1), sizeof(dfloat));
  memcpy(extbdf_a, _a, Nstages*Nstages*sizeof(dfloat));
  memcpy(extbdf_b, _b, Nstages*(Nstages+1)*sizeof(dfloat));

  o_extbdf_a = device.malloc(Nstages*Nstages*sizeof(dfloat), extbdf_a);
  o_extbdf_b = device.malloc(Nstages*(Nstages+1)*sizeof(dfloat), extbdf_b);
}

dfloat extbdf3::getGamma() {
  return *(extbdf_b + (Nstages-1)*(Nstages+1)); //first entry of last row of B
}

void extbdf3::Run(occa::memory &o_q, dfloat start, dfloat end) {

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

void extbdf3::Step(occa::memory &o_q, dfloat time, dfloat _dt, int order) {

  //F(q) at current index
  occa::memory o_F0 = o_F + shiftIndex*N*sizeof(dfloat);

  //coefficients at current order
  occa::memory o_A = o_extbdf_a + order*Nstages*sizeof(dfloat);
  occa::memory o_B = o_extbdf_b + order*(Nstages+1)*sizeof(dfloat);
  dfloat *B = extbdf_b + order*(Nstages+1);

  //evaluate explicit part of rhs: F(q,t)
  solver.rhs_imex_f(o_q, o_F0, time);

  //build rhs for implicit step and update history
  rhsKernel(N,
           _dt,
           shiftIndex,
           o_A,
           o_B,
           o_q,
           o_F,
           o_qn,
           o_rhs);

  dfloat gamma = B[0]/_dt;

  //solve implicit part:
  // find q such that gamma*q - G(q) = rhs
  int iter = solver.rhs_imex_invg(o_rhs, o_q, gamma, time+_dt);

  //rotate index
  shiftIndex = (shiftIndex+Nstages-1)%Nstages;
}

extbdf3::~extbdf3() {
  if (o_rhs.size()) o_rhs.free();
  if (o_qn.size()) o_qn.free();
  if (o_F.size()) o_F.free();
  if (o_extbdf_a.size()) o_extbdf_a.free();
  if (o_extbdf_b.size()) o_extbdf_b.free();

  if (extbdf_a) free(extbdf_a);
  if (extbdf_b) free(extbdf_b);

  rhsKernel.free();
}

} //namespace TimeStepper
