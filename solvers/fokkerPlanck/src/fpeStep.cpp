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

#include "fpe.hpp"

dfloat fpe_t::MaxWaveSpeed(deviceMemory<dfloat>& o_Q, const dfloat T){

  //Note: if this is on the critical path in the future, we should pre-allocate this
  deviceMemory<dfloat> o_maxSpeed = platform.malloc<dfloat>(mesh.Nelements);

  maxWaveSpeedKernel(mesh.Nelements,
                     mesh.o_vgeo,
                     mesh.o_sgeo,
                     mesh.o_vmapM,
                     mesh.o_EToB,
                     T,
                     mesh.o_x,
                     mesh.o_y,
                     mesh.o_z,
                     o_Q,
                     o_maxSpeed);

  const dfloat vmax = platform.linAlg().max(mesh.Nelements, o_maxSpeed, mesh.comm);

  return vmax;
}

//evaluate ODE rhs = f(q,t)
void fpe_t::rhsf(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
  Advection(o_Q, o_RHS, T);
  Diffusion(o_Q, o_RHS, T);
}

// Evaluation of rhs f function
void fpe_t::rhs_imex_f(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
  Advection(o_Q, o_RHS, T);
}

// Evaluation of rhs g function
void fpe_t::rhs_imex_g(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
  Diffusion(o_Q, o_RHS, T);
}

// Inversion of diffusion operator
//  Solves gamma*q - mu*Laplacian*q = rhs
void fpe_t::rhs_imex_invg(deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_Q, const dfloat gamma, const dfloat T){

  // rhs = MM*rhs/mu
  diffusionRhsKernel(mesh.Nelements,
                      mesh.o_vmapM,
                      mesh.o_vgeo,
                      mesh.o_sgeo,
                      mesh.o_EToB,
                      mesh.o_D,
                      mesh.o_LIFT,
                      mesh.o_MM,
                      tau,
                      mu,
                      mesh.o_x,
                      mesh.o_y,
                      mesh.o_z,
                      T,
                      o_RHS);

  int maxIter = 5000;
  // int verbose = settings.compareSetting("ELLIPTIC VERBOSE", "TRUE") ? 1 : 0;
  int verbose =0;

  //call the solver to solve -Laplacian*q + lambda*q = rhs
  dfloat tol = (sizeof(dfloat)==sizeof(double)) ? 1.0e-8 : 1.0e-5;
  elliptic.lambda = gamma/mu;
  int iter = elliptic.Solve(linearSolver, o_Q, o_RHS, tol, maxIter, verbose);

  if (mesh.rank==0){
    printf("\rSolver iterations: %3d.  ", iter); fflush(stdout);
  }
}

// Evolve rhs f function via a sub-timestepper
void fpe_t::rhs_subcycle_f(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_QHAT,
                           const dfloat T, const dfloat dt, const memory<dfloat> B,
                           const int order, const int shiftIndex, const int maxOrder) {

  //subcycle each Lagrangian state qhat by stepping dqhat/dt = F(qhat,t)

  //At each iteration of n, we step the partial sum
  // sum_i=n^order B[i]*Q(t-i*dt) from t-n*dt to t-(n-1)*dt
  //To keep BCs consistent, we step the scaled partial sum:
  // QHAT = sum_i=n^order B[i]*Q(t-i*dt)/(sum_i=n^order B[i])
  dfloat bSum = 0.0;

  dlong N = mesh.Nelements*mesh.Np;

  for (int n=order;n>=0;n--) { //for each history state, starting with oldest

    //q at t-n*dt
    deviceMemory<dfloat> o_Qn = o_Q + ((shiftIndex+n)%maxOrder)*N;

    //next scaled partial sum
    platform.linAlg().axpy(N, B[n+1]/(B[n+1]+bSum), o_Qn,
                              bSum/(B[n+1]+bSum), o_QHAT);
    bSum += B[n+1];

    subStepper.Run(subcycler, o_QHAT, T-n*dt, T-(n-1)*dt);
  }
}

void fpe_t::Advection(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T) {

  // extract q halo on DEVICE
  traceHalo.ExchangeStart(o_Q, 1);

  if (cubature)
    advectionVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubD,
                         mesh.o_cubPDT,
                         mesh.o_cubInterp,
                         mesh.o_cubProject,
                         T,
                         mesh.o_cubx,
                         mesh.o_cuby,
                         mesh.o_cubz,
                         o_Q,
                         o_RHS);
  else
    advectionVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_D,
                         T,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         o_Q,
                         o_RHS);

  traceHalo.ExchangeFinish(o_Q, 1);

  if (cubature)
    advectionSurfaceKernel(mesh.Nelements,
                          mesh.o_vgeo,
                          mesh.o_cubsgeo,
                          mesh.o_vmapM,
                          mesh.o_vmapP,
                          mesh.o_EToB,
                          mesh.o_intInterp,
                          mesh.o_intLIFT,
                          T,
                          mesh.o_intx,
                          mesh.o_inty,
                          mesh.o_intz,
                          o_Q,
                          o_RHS);
  else
    advectionSurfaceKernel(mesh.Nelements,
                          mesh.o_sgeo,
                          mesh.o_LIFT,
                          mesh.o_vmapM,
                          mesh.o_vmapP,
                          mesh.o_EToB,
                          T,
                          mesh.o_x,
                          mesh.o_y,
                          mesh.o_z,
                          o_Q,
                          o_RHS);
}

void fpe_t::Diffusion(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T) {

  //compute gradq and pack with q
  gradientKernel(mesh.Nelements,
                  mesh.o_vgeo,
                  mesh.o_D,
                  o_Q,
                  o_grad);

  traceHalo.ExchangeStart(o_grad, 4);

  if(mesh.NinternalElements)
    diffusionKernel(mesh.NinternalElements,
                    mesh.o_internalElementIds,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    tau,
                    mesh.o_vgeo,
                    mesh.o_sgeo,
                    mesh.o_EToB,
                    mesh.o_D,
                    mesh.o_LIFT,
                    T,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    mu,
                    o_grad,
                    o_RHS);

  traceHalo.ExchangeFinish(o_grad, 4);

  if(mesh.NhaloElements)
    diffusionKernel(mesh.NhaloElements,
                    mesh.o_haloElementIds,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    tau,
                    mesh.o_vgeo,
                    mesh.o_sgeo,
                    mesh.o_EToB,
                    mesh.o_D,
                    mesh.o_LIFT,
                    T,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    mu,
                    o_grad,
                    o_RHS);
}
