#include <unistd.h>
#include <stdio.h>
#include "occa.hpp"

#ifdef OCCA_VERSION_1_0

#include <occa/modes/opencl.hpp>
#include <occa/modes/cuda.hpp>
#include <occa/modes/hip.hpp>

void *occaHostMallocPinned(occa::device &device, size_t size, void *source, occa::memory &mem){

  occa::properties props;

#if OCCA_CUDA_ENABLED
  if(device.mode()=="CUDA"){
    props["cuda/mapped"] = true;
    mem = device.malloc(size, source, props);

    void *ptr = occa::cuda::getMappedPtr(mem);

    return ptr;
  }
#endif

#if OCCA_OPENCL_ENABLED
  if(device.mode()=="OPENCL"){
    props["opencl/mapped"] = true;

    mem = device.malloc(size, source, props);

    void *ptr = occa::opencl::getMappedPtr(mem);
    return ptr;
  }
#endif

#if OCCA_HIP_ENABLED
  if(device.mode()=="HIP"){
    props["hip/mapped"] = true;

    mem = device.malloc(size, source, props);

    void *ptr = occa::hip::getMappedPtr(mem);
    return ptr;
  }
#endif
  
  return NULL;
  
}
#else

void *occaHostMallocPinned(occa::device &device, size_t size, void *source, occa::memory &mem){

  mem = device.mappedAlloc(size, source);

  void *ptr = mem.getMappedPointer();
  
  return ptr;
}
  
#endif
