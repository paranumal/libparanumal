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

#include "elliptic.hpp"

template class ellipticStoppingCriteria<float>;
template class ellipticStoppingCriteria<double>;

template <typename T>
ellipticStoppingCriteria<T>::ellipticStoppingCriteria(elliptic_t *_elliptic, occa::kernel *_addBCKernel){
  
  elliptic = _elliptic;

  addBCKernel = _addBCKernel;
  
  mesh_t &mesh = elliptic->mesh;
  
  occa::properties kernelInfo = mesh.props; //copy base occa properties

  int blockMax = 256;
  if (elliptic->platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1,blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  blockSize=256;
  kernelInfo["defines/" "p_blockSize"] = blockSize;
  
  std::string fileName;
  std::string kernelName;
  std::string suffix = mesh.GetSuffix();
  std::string oklFilePrefix = DELLIPTIC "/okl/";
  std::string oklFileSuffix = ".okl";
  
  //  kernelInfo["defines/pfloat"] = pfloatString;
  //  std::cout << "ellipticStoppingCriteria build " << kernelInfo << std::endl;
  
  occa::properties kernelInfoNp = kernelInfo;
  occa::properties kernelInfoNfp = kernelInfo;

  fileName   = oklFilePrefix + "/ellipticStrongResidual" + suffix + oklFileSuffix;
  kernelName = "ellipticStrongResidual" + suffix;
  strongVolumeResidualKernel = elliptic->platform.buildKernel(fileName, kernelName, kernelInfoNp);

  fileName   = oklFilePrefix + "/ellipticStoppingCriteria" + oklFileSuffix;
  kernelName = "ellipticStoppingCriteria";
  stoppingCriteriaKernel = elliptic->platform.buildKernel(fileName, kernelName, kernelInfoNp);

  o_qL = elliptic->platform.malloc<T>(mesh.Nelements*mesh.Np);
  o_bL = elliptic->platform.malloc<T>(mesh.Nelements*mesh.Np);
  o_RnL = elliptic->platform.malloc<T>(mesh.Nelements*mesh.Np);

  // TW: this needs to be checked
  if(elliptic->disc_c0){
    Ngall = elliptic->ogsMasked.Ngather +  elliptic->gHalo.Nhalo;
  }else{
    Ngall = mesh.Nelements*mesh.Np;
  }
  Nblocks = (Ngall+blockSize-1)/blockSize;
  
  o_Rn     = elliptic->platform.malloc<T>(Ngall);
  o_normRn = elliptic->platform.malloc<T>(Nblocks);
  o_normFn = elliptic->platform.malloc<T>(Nblocks);

  normRn.malloc(Nblocks);
  normFn.malloc(Nblocks);
  
  //  int maxIterations = 5000; // this needs to be made robust
  //  etaHistory = (T*) calloc(maxIterations, sizeof(T));

  //setup cubature for error estimate (make this optional later)
  mesh.CubatureSetup();
  mesh.CubaturePhysicalNodes();

  occa::properties kernelInfoCubNp = kernelInfo;

  dlong p_NblockC = (mesh.dim==3) ? 1:(512/(mesh.cubNp));
    
  kernelInfoCubNp["defines/" "p_NblockC"] = p_NblockC;
  kernelInfoCubNp["defines/" "p_cubNp"]   = mesh.cubNp;
  kernelInfoCubNp["defines/" "p_cubNq"]   = mesh.cubNq;

  kernelInfoCubNp["defines/" "p_NblockC"] = p_NblockC;

  // TW: H1 error estimate stuff

  NblocksC = (mesh.Nelements+p_NblockC-1)/p_NblockC;
  o_errH1 = elliptic->platform.malloc<T>(NblocksC);
  
  std::string dataFileName;
  elliptic->settings.getSetting("DATA FILE", dataFileName);
  kernelInfoCubNp["includes"] += dataFileName;

  //  std::cout  << kernelInfoCubNp;

#if 0
  fileName   = oklFilePrefix + "/ellipticCubatureH1Error" + suffix + oklFileSuffix;
  kernelName = "ellipticCubatureH1Error" + suffix;
  cubatureH1ErrorKernel = elliptic->platform.buildKernel(fileName, kernelName, kernelInfoCubNp);
#endif
}

template <typename T>
void ellipticStoppingCriteria<T>::setLocalRHS(deviceMemory<T> &o_bLin){
  o_bL.copyFrom(o_bLin);
}

template <typename T>
void ellipticStoppingCriteria<T>::reset(){

  //  eta_old = 1;  
  eta = 0;

  if(elliptic->disc_c0 && addBCKernel){
    //fill masked nodes with BC data
    addBCKernel[0](elliptic->mesh.Nelements,
		   elliptic->mesh.o_x,
		   elliptic->mesh.o_y,
		   elliptic->mesh.o_z,
		   elliptic->o_mapB,
		   o_qL);
  }
  
  currentIteration = 0;
}

template <typename T>
T ellipticStoppingCriteria<T>::errorEstimate(deviceMemory<T> &o_q, deviceMemory<T> &o_r){

  mesh_t &mesh = elliptic->mesh;

  if(elliptic->disc_c0){  
    elliptic->ogsMasked.Scatter(o_qL, o_q, 1, ogs::NoTrans); // what is the 1 ?
  }else{
    o_q.copyTo(o_qL);
  }
  
  // compute local components of etaV2 = sum_e  h^2_e*|| lambda*q - laplacian*q - f ||^2_e
  strongVolumeResidualKernel(mesh.Nelements,
			     mesh.o_wJ,
			     mesh.o_ggeo,
			     mesh.o_vgeo, // this for tensor-product elements
			     mesh.o_D,
			     mesh.o_strongS,
			     mesh.o_MM,
			     elliptic->lambda,
			     o_qL,
			     o_bL,
			     o_RnL);

  if(elliptic->disc_c0){  
    // TW: double check this ??
    elliptic->ogsMasked.Gather(o_Rn, o_RnL, 1, ogs::Add, ogs::Trans);

    // max(|| o_Rn ||_{l2}, || o_Rn - o_r ||_{l2} )
    stoppingCriteriaKernel(Nblocks, Ngall, o_Rn, o_r, o_normRn, o_normFn);
  }else{
    // max(|| o_Rn ||_{l2}, || o_Rn - o_r ||_{l2} )
    stoppingCriteriaKernel(Nblocks, Ngall, o_RnL, o_r, o_normRn, o_normFn);
  }

  // probably should stack normRn and normFn
  T normRn2 = elliptic->platform.linAlg().sum(Nblocks, o_normRn, mesh.comm);
  T normFn2 = elliptic->platform.linAlg().sum(Nblocks, o_normFn, mesh.comm);

//  std::cout << "normRn2, normFn2 = " << sqrt(normRn2) << ", " << sqrt(normFn2) << std::endl;
  
  eta = std::max(normRn2, normFn2); // as iter=>infty, normRn2=>normF2 since r=>0
  eta = sqrt(eta);
  //  eta = sqrt(normRn2+normFn2);

  return eta;
}

// note bL excludes Neumann data
template <typename T>
int ellipticStoppingCriteria<T>::stopTest(int iteration,
                                       deviceMemory<T> &o_q,
                                       deviceMemory<T> &o_r,
                                       T rdotr, T TOL){

  mesh_t &mesh = elliptic->mesh;

  {
#if 0
    memory<T> errH1tmp(NblocksC);
    o_errH1.copyTo(errH1tmp);
    for(int n=0;n<NblocksC;++n){
      std::cout << "errH1(" << n <<  ")=" << errH1tmp[n] << std::endl;
    }

    memory<T> cubw(mesh.cubNp);
    mesh.o_cubw.copyTo(cubw);
    for(int n=0;n<mesh.cubNp;++n){
      std::cout << "from H1 err cubw(" << n <<  ")=" << cubw[n] << std::endl;
    }
#endif
  }

  
  if(rdotr<=TOL){
    printf("stopping due to small rdotr\n");
    return 1;
  }

  int stopStep = 1;
  elliptic->settings.getSetting("STOPPING CRITERIA STEP", stopStep);
  if((iteration%stopStep)){
    return 0;
  }

  // only recompute if eta does not look converged
  T etaConvergenceFactor = 1e-1;
  //  printf("errorNorm: eta=%g, eta_old=%g\n", eta, eta_old);
  //  if(fabs(eta-eta_old)>etaConvergenceFactor*eta){
  {
    //    eta_old = eta;
    eta = errorEstimate(o_q, o_r);

    //    etaHistory[currentIteration++] = eta;
    
  }

  // https://onlinelibrary.wiley.com/doi/pdf/10.1002/cnm.1120
  //    T fudgeFactor = 7.5e-3;
  //    T fudgeFactor = 0.01;
  //  T fudgeFactor = 0.1; // 0.03 2d
  T fudgeFactor = 1.e-1; // 0.03 2d

  //  printf("estimate eta: %g\n", eta);
  elliptic->platform.linAlg().set(NblocksC, (T)0.0, o_errH1);

#if 0
  cubatureH1ErrorKernel(mesh.Nelements,
			mesh.o_cubw, mesh.o_cubx, mesh.o_cuby, mesh.o_cubz, mesh.o_cubwJ,
			mesh.o_vgeo, mesh.o_D, mesh.o_cubInterp, elliptic->lambda,
			o_qL, o_errH1);
  
  T errH1 = elliptic->platform.linAlg().sum(NblocksC, o_errH1, mesh.comm);
#else
  T errH1 = 0;
#endif
  int recommendStop = sqrt(rdotr)<(fudgeFactor*eta);
  printf("CG (Stopping Criteria): iter=%d, errH1=%g, norm(r)=%g, eta=%g, recommendStop=%d\n",
	 iteration, sqrt(errH1), sqrt(rdotr), eta, recommendStop);
  
  //  if(sqrt(rdotr)<(fudgeFactor*eta) && fabs(eta-eta_old)<etaConvergenceFactor*eta_old){
  if(recommendStop){
    // TW WHAT TO DO HERE ?
    printf("stopping due to sqrt(rdotr)  %g < %g*eta %g \n", sqrt(rdotr), fudgeFactor, eta);
    //    printf("STOPPING CONDITIONI WOULD STOP HERE\n");
    //    return 1;
    return 1;
  }

  //  eta_old = eta;

  return 0;
}

