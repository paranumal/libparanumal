/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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
#include "parAlmond/parAlmondMultigrid.hpp"

namespace parAlmond {

void multigrid_t::kcycle(const int k, occa::memory& o_RHS, occa::memory& o_X){

  //check for base level
  if(k==baseLevel) {
    coarseSolver->solve(o_RHS, o_X);
    return;
  }

  multigridLevel *level  = levels[k];
  multigridLevel *levelC = levels[k+1];
  occa::memory& o_RHSC = o_rhs[k+1];
  occa::memory& o_XC   = o_x[k+1];
  occa::memory& o_RES   = o_scratch;

  const dlong mCoarse = levelC->Nrows;

  //apply smoother to x and then compute res = rhs-Ax
  level->smooth(o_RHS, o_X, true);

  level->residual(o_RHS, o_X, o_RES);

  // rhsC = P^T res
  level->coarsen(o_RES, o_RHSC);

  if(k+1>NUMKCYCLES) {
    vcycle(k+1, o_RHSC, o_XC);
  } else{
    // first inner krylov iteration
    kcycle(k+1, o_RHSC, o_XC);

    occa::memory& o_CK   = o_ck[k+1];
    occa::memory& o_VK   = o_vk[k+1];
    occa::memory& o_WK   = o_wk[k+1];

    // ck = xC, vk = A*ck
    // alpha1=ck*rhsC, rho1=ck*Ack, norm_rhs=sqrt(rhsC*rhsC)
    // rhsC = rhsC - (alpha1/rho1)*vkp1
    // norm_rtilde = sqrt(rhsC*rhsC)
    dfloat rho1, alpha1, norm_rhs, norm_rhstilde;
    kcycleOp1(levelC, o_XC, o_RHSC, o_CK, o_VK,
              &alpha1, &rho1, &norm_rhs, &norm_rhstilde);

    if(norm_rhstilde < KCYCLETOL*norm_rhs){
      // xC = (alpha1/rho1)*xC
      platform.linAlg.scale(mCoarse, alpha1/rho1, o_XC);
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
  level->prolongate(o_XC, o_X);

  level->smooth(o_RHS, o_X, false);
}


void multigrid_t::kcycleOp1(multigridLevel* level,
                           occa::memory& o_X,  occa::memory& o_RHS,
                           occa::memory& o_CK, occa::memory& o_VK,
                           dfloat *alpha1, dfloat *rho1,
                           dfloat *norm_rhs, dfloat *norm_rhstilde) {

  //ck = x
  platform.linAlg.axpy(level->Nrows, 1.0, o_X, 0.0, o_CK);

  // vk = A*ck
  level->Operator(o_CK,o_VK);

  // alpha1=ck*rhsC, rho1=ck*Ack, norm_rhs=sqrt(rhsC*rhsC)
  dfloat rho[3];
  if(ktype == PCG)
    kcycleCombinedOp1(level, rho, o_CK, o_RHS, o_VK);

  if(ktype == GMRES)
    kcycleCombinedOp1(level, rho, o_VK, o_RHS, o_VK);

  *alpha1 = rho[0];
  *rho1   = rho[1];
  *norm_rhs = sqrt(rho[2]);

  // rhs = rhs - (alpha1/rho1)*vk
  const dfloat a = -(*alpha1)/(*rho1);
  *norm_rhstilde = sqrt(vectorAddInnerProd(level, a, o_VK, 1.0, o_RHS));
}

void multigrid_t::kcycleOp2(multigridLevel* level,
                            occa::memory& o_X,  occa::memory& o_RHS,
                            occa::memory& o_CK, occa::memory& o_VK, occa::memory& o_WK,
                            const dfloat alpha1, const dfloat rho1) {

  if(fabs(rho1) > (dfloat) 1e-20){
    // wk = A*x
    level->Operator(o_X,o_WK);

    // gamma=xC*Ack, beta=xC*AxC, alpha2=xC*rhsC
    dfloat rho[3];

    if(ktype == PCG)
      kcycleCombinedOp2(level, rho, o_X, o_VK, o_WK, o_RHS);

    if(ktype == GMRES)
      kcycleCombinedOp2(level, rho, o_WK, o_VK, o_WK, o_RHS);

    const dfloat gamma  = rho[0];
    const dfloat beta   = rho[1];
    const dfloat alpha2 = rho[2];


    const dfloat rho2 = beta - gamma*gamma/rho1;

    if(fabs(rho2) > (dfloat) 1e-20){
      // x = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ck + (alpha2/rho2)*x
      const dfloat a = alpha1/rho1 - gamma*alpha2/(rho1*rho2);
      const dfloat b = alpha2/rho2;

      platform.linAlg.axpy(level->Nrows, a, o_CK, b, o_X);
    }
  }
}

// returns aDotbc[0] = a\dot b, aDotbc[1] = a\dot c, aDotbc[2] = b\dot b,
void multigrid_t::kcycleCombinedOp1(multigridLevel* level,
                                    dfloat *aDotbc, occa::memory& o_a,
                                    occa::memory& o_b, occa::memory& o_c) {

  const dlong N = level->Nrows;
  dfloat result[3] = {0.,0.,0.};
  dlong numBlocks = (N < PARALMOND_NBLOCKS) ? N : PARALMOND_NBLOCKS;

  if (level->weighted) {
    kcycleWeightedCombinedOp1Kernel(numBlocks,N,o_a,o_b,o_c,level->o_weight,
                                    o_reductionScratch);
  } else {
    kcycleCombinedOp1Kernel(numBlocks,N,o_a,o_b,o_c,o_reductionScratch);
  }
  o_reductionScratch.copyTo(reductionScratch,3*numBlocks*sizeof(dfloat),0);

  for(dlong i=0; i<numBlocks; i++) {
    result[0] += ((dfloat*)reductionScratch)[3*i+0];
    result[1] += ((dfloat*)reductionScratch)[3*i+1];
    result[2] += ((dfloat*)reductionScratch)[3*i+2];
  }
  MPI_Allreduce(result,aDotbc,3,MPI_DFLOAT,MPI_SUM,platform.comm);
}

// returns aDotbcd[0] = a\dot b, aDotbcd[1] = a\dot c, aDotbcd[2] = a\dot d,
void multigrid_t::kcycleCombinedOp2(multigridLevel* level, dfloat *aDotbcd,
                                    occa::memory& o_a, occa::memory& o_b,
                                    occa::memory& o_c, occa::memory& o_d) {

  const dlong N = level->Nrows;
  dfloat result[3] = {0.,0.,0.};
  dlong numBlocks = (N < PARALMOND_NBLOCKS) ? N : PARALMOND_NBLOCKS;

  if (level->weighted) {
    kcycleWeightedCombinedOp2Kernel(numBlocks,N,o_a,o_b,o_c,o_d,
                                    level->o_weight,o_reductionScratch);
  } else {
    kcycleCombinedOp2Kernel(numBlocks,N,o_a,o_b,o_c,o_d,o_reductionScratch);
  }
  o_reductionScratch.copyTo(reductionScratch,3*numBlocks*sizeof(dfloat),0);

  for(dlong i=0; i<numBlocks; i++) {
    result[0] += ((dfloat*)reductionScratch)[3*i+0];
    result[1] += ((dfloat*)reductionScratch)[3*i+1];
    result[2] += ((dfloat*)reductionScratch)[3*i+2];
  }
  MPI_Allreduce(result,aDotbcd,3,MPI_DFLOAT,MPI_SUM,platform.comm);
}

// y = beta*y + alpha*x, and return y\dot y
dfloat multigrid_t::vectorAddInnerProd(multigridLevel* level,
                                      const dfloat alpha, occa::memory& o_X,
                                      const dfloat beta,  occa::memory& o_Y){

  const dlong N = level->Nrows;
  dfloat result = 0.;
  dfloat gresult = 0.;
  dlong numBlocks = (N < PARALMOND_NBLOCKS) ? N : PARALMOND_NBLOCKS;

  if (level->weighted) {
    vectorAddWeightedInnerProdKernel(numBlocks,N,alpha,beta,o_X,o_Y,
                                      level->o_weight,o_reductionScratch);
  } else {
    vectorAddInnerProdKernel(numBlocks,N,alpha,beta,o_X,o_Y,o_reductionScratch);
  }
  o_reductionScratch.copyTo(reductionScratch,numBlocks*sizeof(dfloat),0);

  // #pragma omp parallel for reduction(+:result)
  for (dlong i=0; i<numBlocks; i++) {
    result += ((dfloat*)reductionScratch)[i];
  }
  MPI_Allreduce(&result,&gresult,1,MPI_DFLOAT,MPI_SUM,platform.comm);
  return gresult;
}

} //namespace parAlmond
