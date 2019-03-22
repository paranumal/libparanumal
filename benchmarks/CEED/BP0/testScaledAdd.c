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


#include <stdio.h>
#include <stdlib.h>

#include "occa.hpp"

#if 1
#define dfloat double
#define dfloatString "double"
#else
#define dfloat float
#define dfloatString "float"
#endif

int main(int argc, char **argv){

  /* HOST INFO */
  int N;
  int maxN = 1024*1024*16;
  int minN = 1024;
  int stepN = 1024*8;
  int blockSize = 1024;
  int maxNblocks = (maxN+blockSize-1)/blockSize;
  
  dfloat *h_a = (dfloat*) calloc(maxN, sizeof(dfloat));
  dfloat *h_b = (dfloat*) calloc(maxN, sizeof(dfloat));
  dfloat *h_adotb = (dfloat*) calloc(maxNblocks, sizeof(dfloat));
  dfloat *h_result = (dfloat*) calloc(maxNblocks, sizeof(dfloat));
  for(int n=0;n<maxN;++n){
    h_a[n] = 1.;
    h_b[n] = 1.;
  }

  dfloat alpha = 1., beta = 1.;
  
  /* DEVICE INFO */
  char deviceConfig[BUFSIZ];
  sprintf(deviceConfig, "mode = CUDA, deviceID = 0");
  
  occa::device device;
  device.setup(deviceConfig);
  
  occa::memory o_a = device.malloc(maxN*sizeof(dfloat), h_a);
  occa::memory o_b = device.malloc(maxN*sizeof(dfloat), h_b);

  /* KERNEL INFO */
  occa::kernelInfo kernelInfo;
  kernelInfo.addDefine("dlong", "int");
  kernelInfo.addDefine("dfloat", dfloatString);
  
  /* KERNEL BUILD */
  occa::kernel scaledAddKernel
    = device.buildKernelFromSource(LIBP_DIR "/okl/scaledAdd.okl",
				   "scaledAdd",
				   kernelInfo);

  printf("N Elapsed Bandwidth\n");
  for(int N=minN;N<=maxN;N+=stepN){

    int Ntests = 10;

    occa::streamTag startTag = device.tagStream();
    
    for(int test=0;test<Ntests;++test)
      scaledAddKernel(N, alpha, o_a, beta, o_b);
    
    occa::streamTag stopTag = device.tagStream();
    
    double elapsed = device.timeBetween(startTag, stopTag)/(double)Ntests;
    dfloat BW = sizeof(dfloat)*(N*3)/(double)(1024*1024*1024*elapsed); // neglect write out
    printf("%d %g %g\n", N, elapsed, BW);
    
  }

  return 0;
}
