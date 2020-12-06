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

#include "ins.hpp"

// Inversion of diffusion operator
//  Solves gamma*U - mu*Laplacian*U = rhs
//  Afterwards, imposes incompressiblity via pressure problem
void ins_t::rhs_imex_invg(occa::memory& o_RHS, occa::memory& o_U, const dfloat gamma, const dfloat T){

  const dfloat dt = timeStepper->GetTimeStep();

  if (pressureIncrement) {
    //use current pressure in velocity RHS
    // RHS = RHS - grad P
    Gradient(-1.0, o_p, 1.0, o_RHS, T);

    //call velocty solver to solve
    // gamma*U - mu*Laplacian*U = RHS
    VelocitySolve(o_U, o_RHS, gamma, T);

    // rhsP = -Div U
    Divergence(-1.0, o_U, 0.0, o_rhsP, T);

    //remove old pressure gradient from U
    // U = U + dt*grad P
    Gradient(dt, o_p, 1.0, o_U, T);

    //call pressure increment solver to solve
    // -dt*Laplacian*PI = rhsP
    // P += PI
    PressureIncrementSolve(o_p, o_rhsP, dt, T, dt);

    //update velocity with new pressure correction
    // U = U - dt*grad P
    Gradient(-dt, o_p, 1.0, o_U, T);

  } else {
    //call velocty solver to solve
    // gamma*U - mu*Laplacian*U = RHS
    VelocitySolve(o_U, o_RHS, gamma, T);

    // rhsP = -Div U
    Divergence(-1.0, o_U, 0.0, o_rhsP, T);

    //call pressure solver to solve
    // -dt*Laplacian*P = rhsP
    PressureSolve(o_p, o_rhsP, dt, T);

    //update velocity with pressure correction
    // U = U - dt*grad P
    Gradient(-dt, o_p, 1.0, o_U, T);
  }


  if (mesh.rank==0 && mesh.dim==2) {
    printf("\nSolver iterations: U - %3d, V - %3d, P - %3d", NiterU, NiterV, NiterP); fflush(stdout);
  } else if (mesh.rank==0 && mesh.dim==3) {
    printf("\nSolver iterations: U - %3d, V - %3d, W - %3d, P - %3d", NiterU, NiterV, NiterW, NiterP); fflush(stdout);
  }
}

// Evaluation of rhs f function
void ins_t::rhs_imex_f(occa::memory& o_U, occa::memory& o_RHS, const dfloat T){
  // RHS = N(U)
  Advection(1.0, o_U, 0.0, o_RHS, T);
}

// Evolve rhs f function via a sub-timestepper
void ins_t::rhs_subcycle_f(occa::memory& o_U, occa::memory& o_UHAT,
                           const dfloat T, const dfloat dt, const dfloat* B,
                           const int order, const int shiftIndex, const int maxOrder) {

  //subcycle each Lagrangian state qhat by stepping dqhat/dt = F(qhat,t)

  if (order>=3)
    LIBP_ABORT("Subcycling supports only order 3 interpolation for now.")

  subcycler->order = order;
  subcycler->maxOrder = maxOrder;
  subcycler->shiftIndex = shiftIndex;
  subcycler->T0 = T;
  subcycler->dt = dt;

  dlong Ntot   = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*NVfields;

  if (cubature) {
    //if computing with cubature, interpolate U to all the cubature
    // nodes once to avoid repeatedly interpolating when time stepping

    dlong cNtot  = mesh.Nelements*mesh.cubNp*NVfields;
    dlong cNftot = (mesh.Nelements+mesh.totalHaloPairs)*mesh.intNfp*mesh.Nfaces*NVfields;

    if (subcycler->o_cUh.size() < maxOrder*cNtot*sizeof(dfloat)
      ||subcycler->o_sUh.size() < maxOrder*cNftot*sizeof(dfloat))
      LIBP_ABORT("Cubature history field allocated too small.")

    // interpolate current U to cubature nodes and insert in sperate history
    occa::memory o_Unow = o_U + shiftIndex*Ntot*sizeof(dfloat);
    occa::memory o_cUnow = subcycler->o_cUh + shiftIndex*cNtot*sizeof(dfloat);
    occa::memory o_sUnow = subcycler->o_sUh + shiftIndex*cNftot*sizeof(dfloat);

    // extract U halo
    vTraceHalo->ExchangeStart(o_Unow, 1, ogs_dfloat);

    //interpolate current U to cubature volume nodes and save
    advectionInterpolationVolumeKernel(mesh.Nelements,
                                       mesh.o_cubvgeo,
                                       mesh.o_cubInterp,
                                       o_Unow,
                                       o_cUnow);

    vTraceHalo->ExchangeFinish(o_Unow, 1, ogs_dfloat);

    //interpolate current U and halo to surface cubature nodes and save
    advectionInterpolationSurfaceKernel(mesh.Nelements+mesh.totalHaloPairs,
                                        mesh.o_cubsgeo,
                                        mesh.o_intInterp,
                                        o_Unow,
                                        o_sUnow);

  } else {
    subcycler->o_Uh = o_U; //history
  }


  //At each iteration of n, we step the partial sum
  // sum_i=n^order B[i]*U(t-i*dt) from t-n*dt to t-(n-1)*dt
  //To keep BCs consistent, we step the scaled partial sum:
  // UHAT = sum_i=n^order B[i]*U(t-i*dt)/(sum_i=n^order B[i])
  dfloat bSum = 0.0;

  dlong N = mesh.Nelements*mesh.Np*NVfields;

  for (int n=order;n>=0;n--) { //for each history state, starting with oldest

    //q at t-n*dt
    occa::memory o_Un = o_U + ((shiftIndex+n)%maxOrder)*Ntot*sizeof(dfloat);

    //next scaled partial sum
    linAlg.axpy(N, B[n+1]/(B[n+1]+bSum), o_Un,
                   bSum/(B[n+1]+bSum), o_UHAT);
    bSum += B[n+1];

    subStepper->Run(o_UHAT, T-n*dt, T-(n-1)*dt);
  }
}
