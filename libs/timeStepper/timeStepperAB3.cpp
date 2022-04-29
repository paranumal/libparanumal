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

/* Adams Bashforth, order 3 */
ab3::ab3(dlong Nelements, dlong NhaloElements,
         int Np, int Nfields,
         platform_t& _platform, comm_t _comm):
  timeStepperBase_t(Nelements, NhaloElements, Np, Nfields,
                    _platform, _comm) {

  Nstages = 3;
  shiftIndex = 0;

  memory<dfloat> rhsq(Nstages*N,0.0);
  o_rhsq = platform.malloc<dfloat>(rhsq);

  properties_t kernelInfo = platform.props(); //copy base occa properties from solver

  const int blocksize=256;

  kernelInfo["defines/" "p_blockSize"] = blocksize;
  kernelInfo["defines/" "p_Nstages"] = Nstages;

  updateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperAB.okl",
                                    "abUpdate",
                                    kernelInfo);

  // initialize AB time stepping coefficients
  dfloat _ab_a[Nstages*Nstages] = {
                           1.0,      0.0,    0.0,
                         3./2.,   -1./2.,    0.0,
                       23./12., -16./12., 5./12.};

  ab_a.malloc(Nstages*Nstages);
  ab_a.copyFrom(_ab_a);

  o_ab_a = platform.malloc<dfloat>(ab_a);
}

void ab3::Run(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat start, dfloat end) {

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

void ab3::Step(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat time, dfloat _dt, int order) {

  //rhs at current index
  deviceMemory<dfloat> o_rhsq0 = o_rhsq + shiftIndex*N;

  //A coefficients at current order
  deviceMemory<dfloat> o_A = o_ab_a + order*Nstages;

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

/**************************************************/
/* PML version                                    */
/**************************************************/

/* Adams Bashforth, order 3 */
ab3_pml::ab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
                int Np, int Nfields, int Npmlfields,
                platform_t& _platform, comm_t _comm):
  ab3(Nelements, NhaloElements, Np, Nfields, _platform, _comm),
  Npml(NpmlElements*Np*Npmlfields) {

  if (Npml) {
    memory<dfloat> pmlq(Npml,0.0);
    o_pmlq   = platform.malloc<dfloat>(pmlq);
    o_rhspmlq = platform.malloc<dfloat>(Nstages*Npml);
  }
}

void ab3_pml::Step(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat time, dfloat _dt, int order) {

  //rhs at current index
  deviceMemory<dfloat> o_rhsq0 = o_rhsq + shiftIndex*N;
  deviceMemory<dfloat> o_rhspmlq0;
  if (Npml)    o_rhspmlq0 = o_rhspmlq + shiftIndex*Npml;

  //A coefficients at current order
  deviceMemory<dfloat> o_A = o_ab_a + order*Nstages;

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

} //namespace TimeStepper

} //namespace libp
