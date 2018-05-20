#include "agmg.h"

void buildAlmondKernels(parAlmond_t *parAlmond){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  occa::kernelInfo defs;

  defs.addDefine("bdim", AGMGBDIM);
  defs.addDefine("simd", SIMDWIDTH);

  if(sizeof(dlong)==4){
    defs.addDefine("dlong","int");
  }
  if(sizeof(dlong)==8){
    defs.addDefine("dlong","long long int");
  }

  if(sizeof(dfloat) == sizeof(double)){
    defs.addDefine("dfloat", "double");
    defs.addDefine("dfloat4", "double4");
  }
  else if(sizeof(dfloat) == sizeof(float)){
    defs.addDefine("dfloat", "float");
    defs.addDefine("dfloat4", "float4");
  }

  defs.addInclude(DPWD "/okl/twoPhaseReduction.h");

  if(parAlmond->device.mode()=="OpenCL")
    parAlmond->device.setCompilerFlags("-cl-opt-disable");

  if(parAlmond->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    defs.addCompilerFlag("--ftz=true");
    defs.addCompilerFlag("--prec-div=false");
    defs.addCompilerFlag("--prec-sqrt=false");
    defs.addCompilerFlag("--use_fast_math");
    defs.addCompilerFlag("--fmad=true"); // compiler option for cuda
  }

  if (rank==0) printf("Compiling parALMOND Kernels \n");

  for (int r=0;r<size;r++) {
    if (r==rank) {
      parAlmond->ellAXPYKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/ellAXPY.okl",
           "ellAXPY", defs);

      parAlmond->ellZeqAXPYKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/ellAXPY.okl",
              "ellZeqAXPY", defs);

      parAlmond->ellJacobiKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/ellAXPY.okl",
              "ellJacobi", defs);

      parAlmond->cooAXKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/cooAX.okl",
             "cooAXKernel", defs);

      parAlmond->scaleVectorKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/scaleVector.okl",
             "scaleVectorKernel", defs);

      parAlmond->sumVectorKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/sumVector.okl",
             "sumVectorKernel", defs);

      parAlmond->addScalarKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/addScalar.okl",
             "addScalarKernel", defs);

      parAlmond->vectorAddKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/vectorAdd.okl",
             "vectorAddKernel", defs);

      parAlmond->vectorAddKernel2 = parAlmond->device.buildKernelFromSource(DPWD "/okl/vectorAdd.okl",
              "vectorAddKernel2", defs);

      parAlmond->setVectorKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/setVector.okl",
              "setVectorKernel", defs);

      parAlmond->dotStarKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/dotStar.okl",
               "dotStarKernel", defs);

      parAlmond->simpleDotStarKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/dotStar.okl",
               "simpleDotStarKernel", defs);

      parAlmond->haloExtract = parAlmond->device.buildKernelFromSource(DPWD "/okl/haloExtract.okl",
                "haloExtract", defs);

      parAlmond->agg_interpolateKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/agg_interpolate.okl",
                 "agg_interpolate", defs);

      parAlmond->innerProdKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/innerProduct.okl",
                 "innerProductKernel", defs);

      parAlmond->vectorAddInnerProdKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/vectorAddInnerProduct.okl",
                 "vectorAddInnerProductKernel", defs);

      parAlmond->kcycleCombinedOp1Kernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/kcycleCombinedOp.okl",
                 "kcycleCombinedOp1Kernel", defs);

      parAlmond->kcycleCombinedOp2Kernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/kcycleCombinedOp.okl",
                 "kcycleCombinedOp2Kernel", defs);

      parAlmond->vectorAddWeightedInnerProdKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/vectorAddInnerProduct.okl",
                 "vectorAddWeightedInnerProductKernel", defs);

      parAlmond->kcycleWeightedCombinedOp1Kernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/kcycleCombinedOp.okl",
                 "kcycleWeightedCombinedOp1Kernel", defs);

      parAlmond->kcycleWeightedCombinedOp2Kernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/kcycleCombinedOp.okl",
                 "kcycleWeightedCombinedOp2Kernel", defs);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
