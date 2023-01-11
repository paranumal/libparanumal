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

void multigrid_t::kcycle(const int k, deviceMemory<pfloat>& o_rhs, deviceMemory<pfloat>& o_x){

  //check for base level
  if(k==baseLevel) {
    coarseSolver->solve(o_rhs, o_x);
    return;
  }

  multigridLevel& level  = *levels[k];
  multigridLevel& levelC = *levels[k+1];

  const dlong mCoarse = levelC.Nrows;

  //apply smoother to x and then compute res = rhs-Ax
  level.smooth(o_rhs, o_x, true);

  deviceMemory<pfloat> o_rhsC = platform.reserve<pfloat>(levelC.Ncols);
  deviceMemory<pfloat> o_xC   = platform.reserve<pfloat>(levelC.Ncols);
  deviceMemory<pfloat> o_res  = platform.reserve<pfloat>(level.Ncols);

  level.residual(o_rhs, o_x, o_res);

  // rhsC = P^T res
  level.coarsen(o_res, o_rhsC);
  o_res.free();

  if(k+1>NUMKCYCLES) {
    vcycle(k+1, o_rhsC, o_xC);
  } else{
    deviceMemory<pfloat> o_ck = platform.reserve<pfloat>(levelC.Ncols);
    deviceMemory<pfloat> o_vk = platform.reserve<pfloat>(levelC.Nrows);
    deviceMemory<pfloat> o_wk = platform.reserve<pfloat>(levelC.Nrows);

    // first inner krylov iteration
    kcycle(k+1, o_rhsC, o_xC);

    // ck = xC, vk = A*ck
    // alpha1=ck*rhsC, rho1=ck*Ack, norm_rhs=sqrt(rhsC*rhsC)
    // rhsC = rhsC - (alpha1/rho1)*vkp1
    // norm_rtilde = sqrt(rhsC*rhsC)
    pfloat rho1, alpha1, norm_rhs, norm_rhstilde;
    kcycleOp1(levelC, o_xC, o_rhsC, o_ck, o_vk,
              alpha1, rho1, norm_rhs, norm_rhstilde);

    if(norm_rhstilde < KCYCLETOL*norm_rhs){
      // xC = (alpha1/rho1)*xC
      platform.linAlg().scale(mCoarse, alpha1/rho1, o_xC);
    } else{

      // second inner krylov iteration
      kcycle(k+1, o_rhsC, o_xC);

      // wk = A*xC
      // gamma=xC*Ack, beta=xC*AxC, alpha2=xC*rhsC
      // rho2=beta - gamma*gamma/rho1
      // xC = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ck + (alpha2/rho2)*xC
      kcycleOp2(levelC, o_xC, o_rhsC, o_ck, o_vk, o_wk, alpha1, rho1);
    }
  }

  // x = x + P xC
  level.prolongate(o_xC, o_x);
  o_xC.free(); o_rhsC.free();

  level.smooth(o_rhs, o_x, false);
}


void multigrid_t::kcycleOp1(multigridLevel& level,
                           deviceMemory<pfloat>& o_x,  deviceMemory<pfloat>& o_rhs,
                           deviceMemory<pfloat>& o_ck, deviceMemory<pfloat>& o_vk,
                           pfloat& alpha1, pfloat& rho1,
                           pfloat& norm_rhs, pfloat& norm_rhstilde) {

  //ck = x
  platform.linAlg().axpy(level.Nrows, 1.0, o_x, 0.0, o_ck);

  // vk = A*ck
  level.Operator(o_ck,o_vk);

  // alpha1=ck*rhsC, rho1=ck*Ack, norm_rhs=sqrt(rhsC*rhsC)
  if(ktype == PCG)
    kcycleCombinedOp1(level, o_ck, o_rhs, o_vk, alpha1, rho1, norm_rhs);

  if(ktype == GMRES)
    kcycleCombinedOp1(level, o_vk, o_rhs, o_vk, alpha1, rho1, norm_rhs);

  norm_rhs = sqrt(norm_rhs);

  // rhs = rhs - (alpha1/rho1)*vk
  const pfloat a = -(alpha1)/(rho1);
  norm_rhstilde = sqrt(vectorAddInnerProd(level, a, o_vk, 1.0, o_rhs));
}

void multigrid_t::kcycleOp2(multigridLevel& level,
                            deviceMemory<pfloat>& o_x,  deviceMemory<pfloat>& o_rhs,
                            deviceMemory<pfloat>& o_ck, deviceMemory<pfloat>& o_vk,
                            deviceMemory<pfloat>& o_wk,
                            const pfloat alpha1, const pfloat rho1) {

  if(std::abs(rho1) > (pfloat) 1e-20){
    // wk = A*x
    level.Operator(o_x,o_wk);

    // gamma=xC*Ack, beta=xC*AxC, alpha2=xC*rhsC
    pfloat gamma, beta, alpha2;

    if(ktype == PCG)
      kcycleCombinedOp2(level, o_x, o_vk, o_wk, o_rhs, gamma, beta, alpha2);

    if(ktype == GMRES)
      kcycleCombinedOp2(level, o_wk, o_vk, o_wk, o_rhs, gamma, beta, alpha2);

    const pfloat rho2 = beta - gamma*gamma/rho1;

    if(std::abs(rho2) > (pfloat) 1e-20){
      // x = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ck + (alpha2/rho2)*x
      const pfloat a = alpha1/rho1 - gamma*alpha2/(rho1*rho2);
      const pfloat b = alpha2/rho2;

      platform.linAlg().axpy(level.Nrows, a, o_ck, b, o_x);
    }
  }
}

// returns aDotbc[0] = a\dot b, aDotbc[1] = a\dot c, aDotbc[2] = b\dot b,
void multigrid_t::kcycleCombinedOp1(multigridLevel& level,
                                    deviceMemory<pfloat>& o_a,
                                    deviceMemory<pfloat>& o_b,
                                    deviceMemory<pfloat>& o_c,
                                    pfloat& aDotb,
                                    pfloat& aDotc,
                                    pfloat& bDotb) {

  const dlong N = level.Nrows;
  dlong numBlocks = (N+blockSize-1)/blockSize;
  numBlocks = std::min(numBlocks, PARALMOND_NBLOCKS);

  //pinned scratch buffer
  pinnedMemory<pfloat> h_scratch = platform.hostReserve<pfloat>(3);
  deviceMemory<pfloat> o_scratch = platform.reserve<pfloat>(3*PARALMOND_NBLOCKS);

  kcycleCombinedOp1Kernel(numBlocks,N,o_a,o_b,o_c,o_scratch);

  if (numBlocks>0) {
    h_scratch.copyFrom(o_scratch, 3, 0, properties_t("async", true));
    platform.finish();
  } else {
    h_scratch[0] = 0.0;
    h_scratch[1] = 0.0;
    h_scratch[2] = 0.0;
  }

  comm.Allreduce(h_scratch, Comm::Sum, 3);
  aDotb = h_scratch[0];
  aDotc = h_scratch[1];
  bDotb = h_scratch[2];
}

// returns aDotbcd[0] = a\dot b, aDotbcd[1] = a\dot c, aDotbcd[2] = a\dot d,
void multigrid_t::kcycleCombinedOp2(multigridLevel& level,
                                    deviceMemory<pfloat>& o_a,
                                    deviceMemory<pfloat>& o_b,
                                    deviceMemory<pfloat>& o_c,
                                    deviceMemory<pfloat>& o_d,
                                    pfloat& aDotb,
                                    pfloat& aDotc,
                                    pfloat& aDotd) {

  const dlong N = level.Nrows;
  dlong numBlocks = (N+blockSize-1)/blockSize;
  numBlocks = std::min(numBlocks, PARALMOND_NBLOCKS);

  //pinned scratch buffer
  pinnedMemory<pfloat> h_scratch = platform.hostReserve<pfloat>(3);
  deviceMemory<pfloat> o_scratch = platform.reserve<pfloat>(3*PARALMOND_NBLOCKS);

  kcycleCombinedOp2Kernel(numBlocks,N,o_a,o_b,o_c,o_d,o_scratch);

  if (numBlocks>0) {
    h_scratch.copyFrom(o_scratch, 3, 0, properties_t("async", true));
    platform.finish();
  } else {
    h_scratch[0] = 0.0;
    h_scratch[1] = 0.0;
    h_scratch[2] = 0.0;
  }

  comm.Allreduce(h_scratch, Comm::Sum, 3);
  aDotb = h_scratch[0];
  aDotc = h_scratch[1];
  aDotd = h_scratch[2];
}

// y = beta*y + alpha*x, and return y\dot y
pfloat multigrid_t::vectorAddInnerProd(multigridLevel& level,
                                      const pfloat alpha, deviceMemory<pfloat>& o_x,
                                      const pfloat beta,  deviceMemory<pfloat>& o_Y){

  const dlong N = level.Nrows;
  dlong numBlocks = (N+blockSize-1)/blockSize;
  numBlocks = std::min(numBlocks, PARALMOND_NBLOCKS);

  //pinned scratch buffer
  pinnedMemory<pfloat> h_scratch = platform.hostReserve<pfloat>(1);
  deviceMemory<pfloat> o_scratch = platform.reserve<pfloat>(PARALMOND_NBLOCKS);

  vectorAddInnerProdKernel(numBlocks,N,alpha,beta,o_x,o_Y,o_scratch);

  if (numBlocks>0) {
    h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
    platform.finish();
  } else {
    h_scratch[0] = 0.0;
  }

  comm.Allreduce(h_scratch, Comm::Sum, 1);
  return h_scratch[0];
}

} //namespace parAlmond

} //namespace libp
