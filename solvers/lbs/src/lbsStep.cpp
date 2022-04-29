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

#include "lbs.hpp"

dfloat lbs_t::MaxWaveSpeed(){
  const dfloat vmax = sqrt(2.0f)*1.f/sqrt(3.f);
  return vmax;
}

//evaluate ODE rhs = f(q,t)
void lbs_t::rhsf(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

  // extract q trace halo and start exchange
  traceHalo.ExchangeStart(o_Q, 1);

  // compute volume contribution to lbs RHS
  rhsVolume(mesh.Nelements, o_Q, o_RHS, T);

  // complete trace halo exchange
  traceHalo.ExchangeFinish(o_Q, 1);

  // compute surface contribution to lbs RHS
  rhsSurface(mesh.Nelements, o_Q, o_RHS, T);
}

void lbs_t::rhsVolume(dlong N, deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

  // compute volume contribution to lbs RHS
  if (N){
    const dfloat dt    = timeStepper.GetTimeStep();
    const dfloat gamma = alpha/timeStepper.GetTimeStep();

    forcingKernel(N,
                  T,
                  dt,
                  gamma,
                  nu,
                  o_LBM,
                  mesh.o_x,
                  mesh.o_y,
                  mesh.o_z,
                  o_Q,
                  o_F,
                  o_U);

    collisionKernel(N,
                    T,
                    dt,
                    gamma,
                    nu,
                    o_LBM,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    o_F,
                    o_U,
                    o_Q);

    volumeKernel(N,
                 mesh.o_vgeo,
                 mesh.o_D,
                 mesh.o_x,
                 mesh.o_y,
                 mesh.o_z,
                 T,
                 nu,
                 gamma,
                 o_LBM,
                 o_Q,
                 o_U,
                 o_RHS);
  }
}


void lbs_t::rhsSurface(dlong N, deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
  
  const dfloat dt = timeStepper.GetTimeStep();
  // // compute volume contribution to lbs RHS
  if (N)
    surfaceKernel(N,
                  mesh.o_sgeo,
                  mesh.o_LIFT,
                  mesh.o_vmapM,
                  mesh.o_vmapP,
                  mesh.o_EToB,
                  mesh.o_x,
                  mesh.o_y,
                  mesh.o_z,
                  dt,
                  T,
                  nu,
                  o_LMAP, 
                  o_LBM,
                  o_F, 
                  o_U,
                  o_Q,
                  o_RHS);
}
