#include <unistd.h>
#include "occa.hpp"

#ifdef OCCA_VERSION_1_0

#include <occa/modes/opencl/utils.hpp>
#include <occa/modes/cuda/utils.hpp>
#include <occa/modes/hip/utils.hpp>

void *occaHostMallocPinned(occa::device &device, size_t size, void *source, occa::memory &mem){

  occa::properties props;

  if(device.mode()=="cuda"){
    props["cuda/mapped"] = true;
    mem = device.malloc(size, source, props);

    void *ptr = occa::cuda::getMappedPtr(mem);
    return ptr;
  }

  if(device.mode()=="opencl"){
    props["opencl/mapped"] = true;

    mem = device.malloc(size, source, props);

    void *ptr = occa::opencl::getMappedPtr(mem);
    return ptr;
  }

  if(device.mode()=="hip"){
    props["hip/mapped"] = true;

    mem = device.malloc(size, source, props);

    void *ptr = occa::hip::getMappedPtr(mem);
    return ptr;
  }

  
  return NULL;
  
}
#else

void *occaHostMallocPinned(occa::device &device, size_t size, void *source, occa::memory &mem){

  mem = device.mappedAlloc(size, source);

  void *ptr = device.getMappedPointer();

  return ptr;
}
  
#endif
