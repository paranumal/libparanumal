/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"
#include "parAlmond/parAlmondKernels.hpp"
#include "parAlmond/parAlmondCoarseSolver.hpp"

namespace libp {

namespace parAlmond {

void multigrid_t::kcycle(const int k, deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_X){

  //check for base level
  if(k==baseLevel) {
    coarseSolver->solve(o_RHS, o_X);
    return;
  }

  multigridLevel& level  = *levels[k];
  multigridLevel& levelC = *levels[k+1];
  deviceMemory<dfloat>& o_RHSC = o_rhs[k+1];
  deviceMemory<dfloat>& o_XC   = o_x[k+1];
  deviceMemory<dfloat>& o_RES  = o_scratch;

  const dlong mCoarse = levelC.Nrows;

  //apply smoother to x and then compute res = rhs-Ax
  level.smooth(o_RHS, o_X, true);

  level.residual(o_RHS, o_X, o_RES);

  // rhsC = P^T res
  level.coarsen(o_RES, o_RHSC);

  if(k+1>NUMKCYCLES) {
    vcycle(k+1, o_RHSC, o_XC);
  } else{
    // first inner krylov iteration
    kcycle(k+1, o_RHSC, o_XC);

    deviceMemory<dfloat>& o_CK   = o_ck[k+1];
    deviceMemory<dfloat>& o_VK   = o_vk[k+1];
    deviceMemory<dfloat>& o_WK   = o_wk[k+1];

    // ck = xC, vk = A*ck
    // alpha1=ck*rhsC, rho1=ck*Ack, norm_rhs=sqrt(rhsC*rhsC)
    // rhsC = rhsC - (alpha1/rho1)*vkp1
    // norm_rtilde = sqrt(rhsC*rhsC)
    dfloat rho1, alpha1, norm_rhs, norm_rhstilde;
    kcycleOp1(levelC, o_XC, o_RHSC, o_CK, o_VK,
              alpha1, rho1, norm_rhs, norm_rhstilde);

    if(norm_rhstilde < KCYCLETOL*norm_rhs){
      // xC = (alpha1/rho1)*xC
      platform.linAlg().scale(mCoarse, alpha1/rho1, o_XC);
    } else{

      // second inner krylov iteration
      kcycle(k+1, o_RHSC, o_XC);

      // wk = A*xC
      // gamma=xC*Ack, beta=xC*AxC, alpha2=xC*rhsC
      // rho2=beta - gamma*gamma/rho1
      // xC = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ck + (alpha2/rho2)*xC
      kcycleOp2(levelC, o_XC, o_RHSC, o_CK, o_VK, o_WK, alpha1, rho1);
    }
  }

  // x = x + P xC
  level.prolongate(o_XC, o_X);

  level.smooth(o_RHS, o_X, false);
}


void multigrid_t::kcycleOp1(multigridLevel& level,
                           deviceMemory<dfloat>& o_X,  deviceMemory<dfloat>& o_RHS,
                           deviceMemory<dfloat>& o_CK, deviceMemory<dfloat>& o_VK,
                           dfloat& alpha1, dfloat& rho1,
                           dfloat& norm_rhs, dfloat& norm_rhstilde) {

  //ck = x
  platform.linAlg().axpy(level.Nrows, 1.0, o_X, 0.0, o_CK);

  // vk = A*ck
  level.Operator(o_CK,o_VK);

  // alpha1=ck*rhsC, rho1=ck*Ack, norm_rhs=sqrt(rhsC*rhsC)
  if(ktype == PCG)
    kcycleCombinedOp1(level, o_CK, o_RHS, o_VK, alpha1, rho1, norm_rhs);

  if(ktype == GMRES)
    kcycleCombinedOp1(level, o_VK, o_RHS, o_VK, alpha1, rho1, norm_rhs);

  norm_rhs = sqrt(norm_rhs);

  // rhs = rhs - (alpha1/rho1)*vk
  const dfloat a = -(alpha1)/(rho1);
  norm_rhstilde = sqrt(vectorAddInnerProd(level, a, o_VK, 1.0, o_RHS));
}

void multigrid_t::kcycleOp2(multigridLevel& level,
                            deviceMemory<dfloat>& o_X,  deviceMemory<dfloat>& o_RHS,
                            deviceMemory<dfloat>& o_CK, deviceMemory<dfloat>& o_VK, deviceMemory<dfloat>& o_WK,
                            const dfloat alpha1, const dfloat rho1) {

  if(std::abs(rho1) > (dfloat) 1e-20){
    // wk = A*x
    level.Operator(o_X,o_WK);

    // gamma=xC*Ack, beta=xC*AxC, alpha2=xC*rhsC
    dfloat gamma, beta, alpha2;

    if(ktype == PCG)
      kcycleCombinedOp2(level, o_X, o_VK, o_WK, o_RHS, gamma, beta, alpha2);

    if(ktype == GMRES)
      kcycleCombinedOp2(level, o_WK, o_VK, o_WK, o_RHS, gamma, beta, alpha2);

    const dfloat rho2 = beta - gamma*gamma/rho1;

    if(std::abs(rho2) > (dfloat) 1e-20){
      // x = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ck + (alpha2/rho2)*x
      const dfloat a = alpha1/rho1 - gamma*alpha2/(rho1*rho2);
      const dfloat b = alpha2/rho2;

      platform.linAlg().axpy(level.Nrows, a, o_CK, b, o_X);
    }
  }
}

// returns aDotbc[0] = a\dot b, aDotbc[1] = a\dot c, aDotbc[2] = b\dot b,
void multigrid_t::kcycleCombinedOp1(multigridLevel& level,
                                    deviceMemory<dfloat>& o_a,
                                    deviceMemory<dfloat>& o_b,
                                    deviceMemory<dfloat>& o_c,
                                    dfloat& aDotb,
                                    dfloat& aDotc,
                                    dfloat& bDotb) {

  const dlong N = level.Nrows;
  dlong numBlocks = std::min(N, PARALMOND_NBLOCKS);

  kcycleCombinedOp1Kernel(numBlocks,N,o_a,o_b,o_c,o_reductionScratch);

  if (numBlocks>0) {
    reductionScratch.copyFrom(o_reductionScratch,3*numBlocks);
  } else {
    reductionScratch[0] = 0.0;
    reductionScratch[1] = 0.0;
    reductionScratch[2] = 0.0;
  }

  for(dlong i=1; i<numBlocks; i++) {
    reductionScratch[0] += reductionScratch[3*i+0];
    reductionScratch[1] += reductionScratch[3*i+1];
    reductionScratch[2] += reductionScratch[3*i+2];
  }
  comm.Allreduce(reductionScratch, Comm::Sum, 3);
  aDotb = reductionScratch[0];
  aDotc = reductionScratch[1];
  bDotb = reductionScratch[2];
}

// returns aDotbcd[0] = a\dot b, aDotbcd[1] = a\dot c, aDotbcd[2] = a\dot d,
void multigrid_t::kcycleCombinedOp2(multigridLevel& level,
                                    deviceMemory<dfloat>& o_a,
                                    deviceMemory<dfloat>& o_b,
                                    deviceMemory<dfloat>& o_c,
                                    deviceMemory<dfloat>& o_d,
                                    dfloat& aDotb,
                                    dfloat& aDotc,
                                    dfloat& aDotd) {

  const dlong N = level.Nrows;
  dlong numBlocks = std::min(N, PARALMOND_NBLOCKS);

  kcycleCombinedOp2Kernel(numBlocks,N,o_a,o_b,o_c,o_d,o_reductionScratch);

  if (numBlocks>0) {
    reductionScratch.copyFrom(o_reductionScratch,3*numBlocks);
  } else {
    reductionScratch[0] = 0.0;
    reductionScratch[1] = 0.0;
    reductionScratch[2] = 0.0;
  }

  for(dlong i=1; i<numBlocks; i++) {
    reductionScratch[0] += reductionScratch[3*i+0];
    reductionScratch[1] += reductionScratch[3*i+1];
    reductionScratch[2] += reductionScratch[3*i+2];
  }
  comm.Allreduce(reductionScratch, Comm::Sum, 3);
  aDotb = reductionScratch[0];
  aDotc = reductionScratch[1];
  aDotd = reductionScratch[2];
}

// y = beta*y + alpha*x, and return y\dot y
dfloat multigrid_t::vectorAddInnerProd(multigridLevel& level,
                                      const dfloat alpha, deviceMemory<dfloat>& o_X,
                                      const dfloat beta,  deviceMemory<dfloat>& o_Y){

  const dlong N = level.Nrows;
  dlong numBlocks = std::min(N, PARALMOND_NBLOCKS);

  vectorAddInnerProdKernel(numBlocks,N,alpha,beta,o_X,o_Y,o_reductionScratch);

  if (numBlocks>0) {
    reductionScratch.copyFrom(o_reductionScratch,numBlocks);
  } else {
    reductionScratch[0] = 0.0;
  }

  // #pragma omp parallel for reduction(+:result)
  for (dlong i=1; i<numBlocks; i++) {
    reductionScratch[0] += reductionScratch[i];
  }
  comm.Allreduce(reductionScratch, Comm::Sum, 1);
  return reductionScratch[0];
}

} //namespace parAlmond

} //namespace libp
