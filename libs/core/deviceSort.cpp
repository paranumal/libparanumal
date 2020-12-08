/* The MIT License (MIT)
 *
 * Copyright (c) 2014-2018 David Medina and Tim Warburton
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 */

#include "mesh.hpp"

// was 512
#define BLOCK_SIZE 1024

void deviceSort_t::sort(const int    entries,
			occa::memory o_list) {

  //Nothing to sort
  if (entries < 2) return;

  bitonicSortSharedKernel(entries, o_list);

  for (int size = 4*BLOCK_SIZE; size < 2*entries; size <<= 1){
    //swap step
    bitonicSwapGlobalKernel(entries, o_list, size);
    
    for (int stride = size / 4; stride > 0; stride >>= 1){
      if (stride >= 2*BLOCK_SIZE) {
        bitonicMergeGlobalKernel(entries, o_list, stride);
      } else {
        bitonicMergeSharedKernel(entries, o_list, stride);
        break;
      }
    }
  }
}

deviceSort_t::deviceSort_t(platform_t &platform, const char *entryType, const char *entryHeader, occa::properties props){

  // Compile the kernel at run-time
  occa::settings()["kernel/verbose"] = true;

  occa::properties kernelInfo = props;
#if 0
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();
#endif
  
  kernelInfo["includes"] += entryType; // "entry.h";
  kernelInfo["includes"] += entryHeader; // "compareEntry.h";
  kernelInfo["defines/BLOCK_SIZE"] = (int)BLOCK_SIZE;
  
  bitonicSortSharedKernel = platform.buildKernel(LIBCORE_DIR "/okl/bitonicSortStructs.okl",
                                               "bitonicSortShared", kernelInfo);
  bitonicSwapGlobalKernel = platform.buildKernel(LIBCORE_DIR "/okl/bitonicSortStructs.okl",
                                               "bitonicSwapGlobal", kernelInfo);
  bitonicMergeGlobalKernel = platform.buildKernel(LIBCORE_DIR "/okl/bitonicSortStructs.okl",
                                               "bitonicMergeGlobal", kernelInfo);
  bitonicMergeSharedKernel = platform.buildKernel(LIBCORE_DIR "/okl/bitonicSortStructs.okl",
                                               "bitonicMergeShared", kernelInfo);
}
