#include "parAlmond.h"

void buildAlmondKernels(parAlmond_t *parAlmond){

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

  if(parAlmond->device.mode()=="OpenCL")
    parAlmond->device.setCompilerFlags("-cl-opt-disable");

  printf("Compiling parALMOND Kernels \n");

  parAlmond->ellAXPYKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/ellAXPY.okl",
		   "ellAXPY", defs);

  parAlmond->ellZeqAXPYKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/ellAXPY.okl",
		      "ellZeqAXPY", defs);

  parAlmond->ellJacobi1Kernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/ellAXPY.okl",
		      "ellJacobi1", defs);

  parAlmond->cooAXKernel1 = parAlmond->device.buildKernelFromSource(DPWD "/okl/cooAX.okl",
			 "spmv_coo_kernel1", defs);

  parAlmond->cooAXKernel2 = parAlmond->device.buildKernelFromSource(DPWD "/okl/cooAX.okl",
			 "spmv_coo_kernel2", defs);

  parAlmond->dcsrAXPYKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/dcsrAXPY.okl",
       "dcsrAXPY", defs);
  
  parAlmond->dcsrZeqAXPYKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/dcsrAXPY.okl",
       "dcsrZeqAXPY", defs);
  
  parAlmond->dcsrJacobiKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/dcsrAXPY.okl",
       "dcsrJacobi", defs);

  parAlmond->copyKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/copy.okl",
		      "copyKernel", defs);

  parAlmond->scaleVectorKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/scaleVector.okl",
			   "scaleVectorKernel", defs);

  parAlmond->partialInnerProdKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/partialInnerProd.okl",
			    "partialInnerProd", defs);

  parAlmond->vectorAddKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/vectorAdd.okl",
			   "vectorAddKernel", defs);

  parAlmond->vectorAddKernel2 = parAlmond->device.buildKernelFromSource(DPWD "/okl/vectorAdd.okl",
			    "vectorAddKernel2", defs);

  parAlmond->dotStarKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/dotStar.okl",
			 "dotStarKernel", defs);

  parAlmond->simpleDotStarKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/dotStar.okl",
			       "simpleDotStarKernel", defs);

  parAlmond->haloExtract = parAlmond->device.buildKernelFromSource(DPWD "/okl/haloExtract.okl",
      "haloExtract", defs);

  parAlmond->agg_interpolateKernel = parAlmond->device.buildKernelFromSource(DPWD "/okl/agg_interpolate.okl",
             "agg_interpolate", defs);
}
