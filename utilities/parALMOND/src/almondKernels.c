#include "parAlmond.h"

void buildAlmondKernels(almond_t *almond){

  occa::kernelInfo defs;

  defs.addDefine("bdim", AGMGBDIM);
  defs.addDefine("simd", SIMDWIDTH);

  if(sizeof(dfloat) == sizeof(double)){
    defs.addDefine("datafloat", "double");
    defs.addDefine("datafloat4", "double4");
  }
  else if(sizeof(dfloat) == sizeof(float)){
    defs.addDefine("datafloat", "float");
    defs.addDefine("datafloat4", "float4");
  }

  defs.addInclude(DPWD "/okl/workgroup_operations.h");
  defs.addInclude(DPWD "/okl/workgroup_segreduce.h");
  defs.addInclude(DPWD "/okl/workgroup_reduce.h");
  defs.addInclude(DPWD "/okl/workgroup_reduce_max.h");

  if(almond->device.mode()=="OpenCL")
    almond->device.setCompilerFlags("-cl-opt-disable");

  printf("Compiling parALMOND Kernels \n");

  almond->ellAXPYKernel = almond->device.buildKernelFromSource(DPWD "/okl/ellAXPY.okl",
		   "ellAXPY", defs);

  almond->ellZeqAXPYKernel = almond->device.buildKernelFromSource(DPWD "/okl/ellAXPY.okl",
		      "ellZeqAXPY", defs);

  almond->ellJacobi1Kernel = almond->device.buildKernelFromSource(DPWD "/okl/ellAXPY.okl",
		      "ellJacobi1", defs);

  almond->cooAXKernel1 = almond->device.buildKernelFromSource(DPWD "/okl/cooAX.okl",
			 "spmv_coo_kernel1", defs);

  almond->cooAXKernel2 = almond->device.buildKernelFromSource(DPWD "/okl/cooAX.okl",
			 "spmv_coo_kernel2", defs);

  almond->copyKernel = almond->device.buildKernelFromSource(DPWD "/okl/copy.okl",
		      "copyKernel", defs);

  almond->scaleVectorKernel = almond->device.buildKernelFromSource(DPWD "/okl/scaleVector.okl",
			   "scaleVectorKernel", defs);

  almond->partialInnerProdKernel = almond->device.buildKernelFromSource(DPWD "/okl/partialInnerProd.okl",
			    "partialInnerProd", defs);

  almond->vectorAddKernel = almond->device.buildKernelFromSource(DPWD "/okl/vectorAdd.okl",
			   "vectorAddKernel", defs);

  almond->vectorAddKernel2 = almond->device.buildKernelFromSource(DPWD "/okl/vectorAdd.okl",
			    "vectorAddKernel2", defs);

  almond->dotStarKernel = almond->device.buildKernelFromSource(DPWD "/okl/dotStar.okl",
			 "dotStarKernel", defs);

  almond->simpleDotStarKernel = almond->device.buildKernelFromSource(DPWD "/okl/dotStar.okl",
			       "simpleDotStarKernel", defs);

  almond->haloExtract = almond->device.buildKernelFromSource(DPWD "/okl/haloExtract.okl",
      "haloExtract", defs);

  almond->agg_interpolateKernel = almond->device.buildKernelFromSource(DPWD "/okl/agg_interpolate.okl",
             "agg_interpolate", defs);
}
