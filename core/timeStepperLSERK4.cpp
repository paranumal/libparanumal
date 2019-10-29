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

lserk4::lserk4(dlong _N, dlong _Nhalo, solver_t& _solver):
  timeStepper_t(_N, _Nhalo, _solver) {};

lserk4::~lserk4() {
  if (o_rhsq.size()) o_rhsq.free();
  if (o_resq.size()) o_resq.free();

  if (rka) free(rka);
  if (rkb) free(rkb);
  if (rkc) free(rkc);

  updateKernel.free();
}

void lserk4::Init() {
  o_resq = device.malloc(N*sizeof(dfloat));
  o_rhsq = device.malloc(N*sizeof(dfloat));

  occa::properties kernelInfo = props; //copy base occa properties from solver

  kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;

  updateKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                    "timeStepperLSERK4.okl",
                                    "lserk4Update",
                                    kernelInfo, comm);

  // initialize LSERK4 time stepping coefficients
  dfloat _rka[5] = {0.0,
       -567301805773.0/1357537059087.0 ,
       -2404267990393.0/2016746695238.0 ,
       -3550918686646.0/2091501179385.0  ,
       -1275806237668.0/842570457699.0};
  dfloat _rkb[5] = { 1432997174477.0/9575080441755.0 ,
        5161836677717.0/13612068292357.0 ,
        1720146321549.0/2090206949498.0  ,
        3134564353537.0/4481467310338.0  ,
        2277821191437.0/14882151754819.0};
  // added one more for advanced time step
  dfloat _rkc[6] = {0.0  ,
       1432997174477.0/9575080441755.0 ,
       2526269341429.0/6820363962896.0 ,
       2006345519317.0/3224310063776.0 ,
       2802321613138.0/2924317926251.0 ,
       1.0};

  Nrk = 5;
  rka = (dfloat*) calloc(Nrk, sizeof(dfloat));
  rkb = (dfloat*) calloc(Nrk, sizeof(dfloat));
  rkc = (dfloat*) calloc(Nrk+1, sizeof(dfloat));
  memcpy(rka, _rka, Nrk*sizeof(dfloat));
  memcpy(rkb, _rkb, Nrk*sizeof(dfloat));
  memcpy(rkc, _rkc, (Nrk+1)*sizeof(dfloat));
}

void lserk4::Run(occa::memory &o_q, dfloat start, dfloat end) {

  dfloat time = start;

  solver.Report(time,0);

  dfloat outputInterval;
  settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat outputTime = time + outputInterval;

  int tstep=0;

  while (time < end) {

    if (time<outputTime && time+dt>=outputTime) {

      //save current state
      occa::memory o_saveq = device.malloc(N*sizeof(dfloat));
      o_saveq.copyFrom(o_q, N*sizeof(dfloat));

      dfloat savedt = dt;

      dt = outputTime-time;

      //take small time step
      Step(o_q, time, dt);

      //report state
      solver.Report(outputTime,tstep);

      //restore previous state
      dt = savedt;
      o_q.copyFrom(o_saveq, N*sizeof(dfloat));
      o_saveq.free();

      outputTime += outputInterval;
    }

    //check for final timestep
    if (time+dt > end){
      dt = end-time;
    }

    Step(o_q, time, dt);
    time += dt;
    tstep++;
  }
}

void lserk4::Step(occa::memory &o_q, dfloat time, dfloat _dt) {

  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(int rk=0;rk<Nrk;++rk){

    dfloat currentTime = time + rkc[rk]*_dt;

    //evaluate ODE rhs = f(q,t)
    solver.rhsf(o_q, o_rhsq, currentTime);

    // update solution using Runge-Kutta
    updateKernel(N, _dt, rka[rk], rkb[rk],
                 o_rhsq, o_resq, o_q);
  }
}
