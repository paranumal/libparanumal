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

#include <unistd.h>
#include <stdio.h>
#include "occa.hpp"

#if 0

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
