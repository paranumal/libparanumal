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

#include "elliptic.h"

void ellipticMultigridSmooth(void **args, occa::memory &o_r, occa::memory &o_x, bool xIsZero) {

  elliptic_t *elliptic = (elliptic_t *) args[0];
  agmgLevel *level = (agmgLevel *) args[1];

  occa::memory o_res = level->o_smootherResidual;

  if (xIsZero) {
    level->device_smoother(level->smootherArgs, o_r, o_x);
    return;
  }

  dfloat one = 1.; dfloat mone = -1.;

  //res = r-Ax
  level->device_Ax(level->AxArgs,o_x,o_res);
  elliptic->scaledAddKernel(level->Nrows,one, o_r, mone, o_res);

  //smooth the fine problem x = x + S(r-Ax)
  level->device_smoother(level->smootherArgs, o_res, o_res);
  elliptic->scaledAddKernel(level->Nrows,one, o_res, one, o_x);
}

void ellipticMultigridSmoothChebyshev(void **args, occa::memory &o_r, occa::memory &o_x, bool xIsZero) {

  elliptic_t *elliptic = (elliptic_t *) args[0];
  agmgLevel *level = (agmgLevel *) args[1];

  dfloat lambdaN = level->smoother_params[0];
  dfloat lambda1 = level->smoother_params[1];

  dfloat theta = 0.5*(lambdaN+lambda1);
  dfloat delta = 0.5*(lambdaN-lambda1);
  dfloat invTheta = 1.0/theta;
  dfloat sigma = theta/delta;
  dfloat rho_n = 1./sigma;
  dfloat rho_np1;

  dfloat one = 1., mone = -1., zero = 0.0;

  occa::memory o_res = level->o_smootherResidual;
  occa::memory o_Ad  = level->o_smootherResidual2;
  occa::memory o_d   = level->o_smootherUpdate;

  if(xIsZero){ //skip the Ax if x is zero
    //res = Sr
    level->device_smoother(level->smootherArgs, o_r, o_res);

    //d = invTheta*res
    elliptic->scaledAddKernel(level->Nrows, invTheta, o_res, zero, o_d);
  } else {
    //res = S(r-Ax)
    level->device_Ax(level->AxArgs,o_x,o_res);
    elliptic->scaledAddKernel(level->Nrows, one, o_r, mone, o_res);
    level->device_smoother(level->smootherArgs, o_res, o_res);

    //d = invTheta*res
    elliptic->scaledAddKernel(level->Nrows, invTheta, o_res, zero, o_d);
  }

  for (int k=0;k<level->ChebyshevIterations;k++) {
    //x_k+1 = x_k + d_k
    if (xIsZero&&(k==0))
      elliptic->scaledAddKernel(level->Nrows, one, o_d, zero, o_x);
    else
      elliptic->scaledAddKernel(level->Nrows, one, o_d, one, o_x);

    //r_k+1 = r_k - SAd_k
    level->device_Ax(level->AxArgs,o_d,o_Ad);
    level->device_smoother(level->smootherArgs, o_Ad, o_Ad);
    elliptic->scaledAddKernel(level->Nrows, mone, o_Ad, one, o_res);

    rho_np1 = 1.0/(2.*sigma-rho_n);
    dfloat rhoDivDelta = 2.0*rho_np1/delta;

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    elliptic->scaledAddKernel(level->Nrows, rhoDivDelta, o_res, rho_np1*rho_n, o_d);

    rho_n = rho_np1;
  }
  //x_k+1 = x_k + d_k
  elliptic->scaledAddKernel(level->Nrows, one, o_d, one, o_x);

}

void LocalPatch(void **args, occa::memory &o_r, occa::memory &o_Sr) {

  elliptic_t *elliptic = (elliptic_t*) args[0];
  mesh_t *mesh = elliptic->mesh;
  precon_t *precon = elliptic->precon;

  occaTimerTic(mesh->device,"approxBlockJacobiSolveKernel");
  precon->approxBlockJacobiSolverKernel(mesh->Nelements,
                            precon->o_patchesIndex,
                            precon->o_invAP,
                            precon->o_invDegreeAP,
                            o_r,
                            o_Sr);
  occaTimerToc(mesh->device,"approxBlockJacobiSolveKernel");
}

void dampedJacobi(void **args, occa::memory &o_r, occa::memory &o_Sr) {

  elliptic_t *elliptic = (elliptic_t *) args[0];
  mesh_t *mesh = elliptic->mesh;

  occa::memory o_invDiagA = elliptic->precon->o_invDiagA;

  elliptic->dotMultiplyKernel(mesh->Np*mesh->Nelements,o_invDiagA,o_r,o_Sr);
}