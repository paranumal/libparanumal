#include "agmg.h"

void buildAlmondKernels(parAlmond_t *parAlmond){

  int rank, size;
  rank = agmg::rank;
  size = agmg::size;
  
  occa::properties kernelInfo;
 kernelInfo["defines"].asObject();
 kernelInfo["includes"].asArray();
 kernelInfo["header"].asArray();
 kernelInfo["flags"].asObject();


  kernelInfo["defines/" "bdim"]= AGMGBDIM;
  kernelInfo["defines/" "simd"]= SIMDWIDTH;

  if(sizeof(dlong)==4){
    kernelInfo["defines/" "dlong"]="int";
  }
  if(sizeof(dlong)==8){
    kernelInfo["defines/" "dlong"]="long long int";
  }

  if(sizeof(dfloat) == sizeof(double)){
    kernelInfo["defines/" "dfloat"]= "double";
    kernelInfo["defines/" "dfloat4"]= "double4";
  }
  else if(sizeof(dfloat) == sizeof(float)){
    kernelInfo["defines/" "dfloat"]= "float";
    kernelInfo["defines/" "dfloat4"]= "float4";
  }

  kernelInfo["defines/" "p_RDIMX"]= RDIMX;
  kernelInfo["defines/" "p_RDIMY"]= RDIMY;

  kernelInfo["includes"] += DPWD "/okl/twoPhaseReduction.h";

  if(parAlmond->device.mode()=="OpenCL"){
    //    parAlmond->device.setCompilerFlags("-cl-opt-disable");
    kernelInfo["compiler_flags"] += "-cl-opt-disable";
  }

  if(parAlmond->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo["compiler_flags"] += "--ftz=true";
    kernelInfo["compiler_flags"] += "--prec-div=false";
    kernelInfo["compiler_flags"] += "--prec-sqrt=false";
    kernelInfo["compiler_flags"] += "--use_fast_math";
    kernelInfo["compiler_flags"] += "--fmad=true"; // compiler option for cuda
  }

  if (rank==0) printf("Compiling parALMOND Kernels \n");

  for (int r=0;r<size;r++) {
    if (r==rank) {
      parAlmond->ellAXPYKernel = parAlmond->device.buildKernel(DPWD "/okl/ellAXPY.okl",
           "ellAXPY", kernelInfo);

      parAlmond->ellZeqAXPYKernel = parAlmond->device.buildKernel(DPWD "/okl/ellAXPY.okl",
              "ellZeqAXPY", kernelInfo);

      parAlmond->ellJacobiKernel = parAlmond->device.buildKernel(DPWD "/okl/ellAXPY.okl",
              "ellJacobi", kernelInfo);

      parAlmond->cooAXKernel = parAlmond->device.buildKernel(DPWD "/okl/cooAX.okl",
             "cooAXKernel", kernelInfo);

      parAlmond->scaleVectorKernel = parAlmond->device.buildKernel(DPWD "/okl/scaleVector.okl",
             "scaleVectorKernel", kernelInfo);

      parAlmond->sumVectorKernel = parAlmond->device.buildKernel(DPWD "/okl/sumVector.okl",
             "sumVectorKernel", kernelInfo);

      parAlmond->addScalarKernel = parAlmond->device.buildKernel(DPWD "/okl/addScalar.okl",
             "addScalarKernel", kernelInfo);

      parAlmond->vectorAddKernel = parAlmond->device.buildKernel(DPWD "/okl/vectorAdd.okl",
             "vectorAddKernel", kernelInfo);

      parAlmond->vectorAddKernel2 = parAlmond->device.buildKernel(DPWD "/okl/vectorAdd.okl",
              "vectorAddKernel2", kernelInfo);

      parAlmond->setVectorKernel = parAlmond->device.buildKernel(DPWD "/okl/setVector.okl",
              "setVectorKernel", kernelInfo);

      parAlmond->dotStarKernel = parAlmond->device.buildKernel(DPWD "/okl/dotStar.okl",
               "dotStarKernel", kernelInfo);

      parAlmond->simpleDotStarKernel = parAlmond->device.buildKernel(DPWD "/okl/dotStar.okl",
               "simpleDotStarKernel", kernelInfo);

      parAlmond->haloExtract = parAlmond->device.buildKernel(DPWD "/okl/haloExtract.okl",
                "haloExtract", kernelInfo);

      parAlmond->agg_interpolateKernel = parAlmond->device.buildKernel(DPWD "/okl/agg_interpolate.okl",
                 "agg_interpolate", kernelInfo);

      parAlmond->innerProdKernel = parAlmond->device.buildKernel(DPWD "/okl/innerProduct.okl",
                 "innerProductKernel", kernelInfo);

      parAlmond->vectorAddInnerProdKernel = parAlmond->device.buildKernel(DPWD "/okl/vectorAddInnerProduct.okl",
                 "vectorAddInnerProductKernel", kernelInfo);

      parAlmond->kcycleCombinedOp1Kernel = parAlmond->device.buildKernel(DPWD "/okl/kcycleCombinedOp.okl",
                 "kcycleCombinedOp1Kernel", kernelInfo);

      parAlmond->kcycleCombinedOp2Kernel = parAlmond->device.buildKernel(DPWD "/okl/kcycleCombinedOp.okl",
                 "kcycleCombinedOp2Kernel", kernelInfo);

      parAlmond->vectorAddWeightedInnerProdKernel = parAlmond->device.buildKernel(DPWD "/okl/vectorAddInnerProduct.okl",
                 "vectorAddWeightedInnerProductKernel", kernelInfo);

      parAlmond->kcycleWeightedCombinedOp1Kernel = parAlmond->device.buildKernel(DPWD "/okl/kcycleCombinedOp.okl",
                 "kcycleWeightedCombinedOp1Kernel", kernelInfo);

      parAlmond->kcycleWeightedCombinedOp2Kernel = parAlmond->device.buildKernel(DPWD "/okl/kcycleCombinedOp.okl",
                 "kcycleWeightedCombinedOp2Kernel", kernelInfo);
    }
    MPI_Barrier(agmg::comm);
  }
}
