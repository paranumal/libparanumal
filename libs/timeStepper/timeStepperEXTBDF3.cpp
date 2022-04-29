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

/* Backward Difference Formula, order 3, with extrapolation */
extbdf3::extbdf3(dlong Nelements, dlong NhaloElements,
                 int Np, int Nfields,
                 platform_t& _platform, comm_t _comm):
  timeStepperBase_t(Nelements, NhaloElements, Np, Nfields,
                    _platform, _comm) {

  Nstages = 3;
  shiftIndex = 0;

  memory<dfloat> qn(Nstages*N, 0.0);
  o_qn = platform.malloc<dfloat>(qn); //q history

  memory<dfloat> rhs(N,0.0);
  o_rhs = platform.malloc<dfloat>(rhs); //rhs storage

  o_F  = platform.malloc<dfloat>(qn); //F(q) history (explicit part)

  properties_t kernelInfo = platform.props(); //copy base occa properties from solver

  const int blocksize=256;

  kernelInfo["defines/" "p_blockSize"] = blocksize;
  kernelInfo["defines/" "p_Nstages"] = Nstages;

  rhsKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperEXTBDF.okl",
                                    "extbdfRHS",
                                    kernelInfo);

  // initialize EXT and BDF time stepping coefficients
  dfloat _a[Nstages*Nstages] = {
                            1.,    0.,     0.,
                            2.,   -1.,     0.,
                            3.,   -3.,     1.};
  dfloat _b[Nstages*(Nstages+1)] = {
                            1.,    1.,     0.,    0.,
                         3./2.,    2., -1./2.,    0.,
                        11./6.,    3., -3./2., 1./3.};

  extbdf_a.malloc(Nstages*Nstages);
  extbdf_b.malloc(Nstages*(Nstages+1));
  extbdf_a.copyFrom(_a);
  extbdf_b.copyFrom(_b);

  o_extbdf_a = platform.malloc<dfloat>(extbdf_a);
  o_extbdf_b = platform.malloc<dfloat>(extbdf_b);
}

dfloat extbdf3::GetGamma() {
  return extbdf_b[(Nstages-1)*(Nstages+1)]; //first entry of last row of B
}

void extbdf3::Run(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat start, dfloat end) {

  dfloat time = start;

  solver.Report(time,0);

  dfloat outputInterval=0.0;
  solver.settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat outputTime = time + outputInterval;

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

void extbdf3::Step(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat time, dfloat _dt, int order) {

  //F(q) at current index
  deviceMemory<dfloat> o_F0 = o_F + shiftIndex*N;

  //coefficients at current order
  deviceMemory<dfloat> o_A = o_extbdf_a + order*Nstages;
  deviceMemory<dfloat> o_B = o_extbdf_b + order*(Nstages+1);
  memory<dfloat> B = extbdf_b + order*(Nstages+1);

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
  solver.rhs_imex_invg(o_rhs, o_q, gamma, time+_dt);

  //rotate index
  shiftIndex = (shiftIndex+Nstages-1)%Nstages;
}

} //namespace TimeStepper

} //namespace libp
