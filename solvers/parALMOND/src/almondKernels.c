#include "agmg.h"

void buildAlmondKernels(parAlmond_t *parAlmond){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  occa::kernelInfo kernelInfo;

  kernelInfo.addDefine("bdim", AGMGBDIM);
  kernelInfo.addDefine("simd", SIMDWIDTH);

  if(sizeof(dlong)==4){
    kernelInfo.addDefine("dlong","int");
  }
  if(sizeof(dlong)==8){
    kernelInfo.addDefine("dlong","long long int");
  }

  if(sizeof(dfloat) == sizeof(double)){
    kernelInfo.addDefine("dfloat", "double");
    kernelInfo.addDefine("dfloat4", "double4");
  }
  else if(sizeof(dfloat) == sizeof(float)){
    kernelInfo.addDefine("dfloat", "float");
    kernelInfo.addDefine("dfloat4", "float4");
  }

  kernelInfo.addDefine("p_RDIMX", RDIMX);
  kernelInfo.addDefine("p_RDIMY", RDIMY);

  kernelInfo.addInclude(DPWD "/okl/twoPhaseReduction.h");

  if(parAlmond->device.mode()=="OpenCL"){
    //    parAlmond->device.setCompilerFlags("-cl-opt-disable");
    kernelInfo.addCompilerFlag("-cl-opt-disable");
  }

  if(parAlmond->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo.addCompilerFlag("--ftz=true");
    kernelInfo.addCompilerFlag("--prec-div=false");
    kernelInfo.addCompilerFlag("--prec-sqrt=false");
    kernelInfo.addCompilerFlag("--use_fast_math");
    kernelInfo.addCompilerFlag("--fmad=true"); // compiler option for cuda
  }

  if (rank==0) printf("Compiling parALMOND Kernels \n");

  for (int r=0;r<size;r++) {
    if (r==rank) {
      parAlmond->ellAXPYKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/ellAXPY.okl",
           "ellAXPY", kernelInfo);

      parAlmond->ellZeqAXPYKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/ellAXPY.okl",
              "ellZeqAXPY", kernelInfo);

      parAlmond->ellJacobiKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/ellAXPY.okl",
              "ellJacobi", kernelInfo);

      parAlmond->cooAXKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/cooAX.okl",
             "cooAXKernel", kernelInfo);

      parAlmond->scaleVectorKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/scaleVector.okl",
             "scaleVectorKernel", kernelInfo);

      parAlmond->sumVectorKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/sumVector.okl",
             "sumVectorKernel", kernelInfo);

      parAlmond->addScalarKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/addScalar.okl",
             "addScalarKernel", kernelInfo);

      parAlmond->vectorAddKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/vectorAdd.okl",
             "vectorAddKernel", kernelInfo);

      parAlmond->vectorAddKernel2 = parAlmond->device.buildKernelFromSource(DPWD "/okl/vectorAdd.okl",
              "vectorAddKernel2", kernelInfo);

      parAlmond->setVectorKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/setVector.okl",
              "setVectorKernel", kernelInfo);

      parAlmond->dotStarKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/dotStar.okl",
               "dotStarKernel", kernelInfo);

      parAlmond->simpleDotStarKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/dotStar.okl",
               "simpleDotStarKernel", kernelInfo);

      parAlmond->haloExtract = parAlmond->device.buildKernelFromSource(DPWD "/okl/haloExtract.okl",
                "haloExtract", kernelInfo);

      parAlmond->agg_interpolateKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/agg_interpolate.okl",
                 "agg_interpolate", kernelInfo);

      parAlmond->innerProdKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/innerProduct.okl",
                 "innerProductKernel", kernelInfo);

      parAlmond->vectorAddInnerProdKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/vectorAddInnerProduct.okl",
                 "vectorAddInnerProductKernel", kernelInfo);

      parAlmond->kcycleCombinedOp1Kernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/kcycleCombinedOp.okl",
                 "kcycleCombinedOp1Kernel", kernelInfo);

      parAlmond->kcycleCombinedOp2Kernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/kcycleCombinedOp.okl",
                 "kcycleCombinedOp2Kernel", kernelInfo);

      parAlmond->vectorAddWeightedInnerProdKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/vectorAddInnerProduct.okl",
                 "vectorAddWeightedInnerProductKernel", kernelInfo);

      parAlmond->kcycleWeightedCombinedOp1Kernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/kcycleCombinedOp.okl",
                 "kcycleWeightedCombinedOp1Kernel", kernelInfo);

      parAlmond->kcycleWeightedCombinedOp2Kernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/kcycleCombinedOp.okl",
                 "kcycleWeightedCombinedOp2Kernel", kernelInfo);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
