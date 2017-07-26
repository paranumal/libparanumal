
#include <stdio.h>
#include <stdlib.h>

#include "occa.hpp"

#define iint int
#define iintString "int"

#if 0
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

  /* DEVICE INFO */
  char deviceConfig[BUFSIZ];
  sprintf(deviceConfig, "mode = CUDA, deviceID = 0");
  
  occa::device device;
  device.setup(deviceConfig);
  
  occa::memory o_a = device.malloc(maxN*sizeof(dfloat), h_a);
  occa::memory o_b = device.malloc(maxN*sizeof(dfloat), h_b);
  occa::memory o_adotb = device.malloc(maxNblocks*sizeof(dfloat), h_adotb);

  /* KERNEL INFO */
  occa::kernelInfo kernelInfo;
  kernelInfo.addDefine("p_blockSize", blockSize);
  kernelInfo.addDefine("iint", iintString);
  kernelInfo.addDefine("dfloat", dfloatString);
  kernelInfo.addParserFlag("automate-add-barriers", "disabled");
  
  /* KERNEL BUILD */
  occa::kernel innerProductKernel
    = device.buildKernelFromSource(DHOLMES "/okl/innerProduct.okl",
				   "innerProductAtomic",
				   kernelInfo);

  printf("N Elapsed Bandwidth\n");
  for(int N=minN;N<=maxN;N+=stepN){
    // make each block load more than one value per thread
    int Nblocks = ((N+3)/4+blockSize-1)/blockSize;
    
    o_adotb.copyFrom(h_adotb);
    int Ntests = 10;

    occa::streamTag startTag = device.tagStream();
    
    for(int test=0;test<Ntests;++test)
      innerProductKernel(N, Nblocks, o_a, o_b, o_adotb);
    
    occa::streamTag stopTag = device.tagStream();

    o_adotb.copyTo(h_result);
    
    double elapsed = device.timeBetween(startTag, stopTag)/(double)Ntests;
    dfloat BW = sizeof(dfloat)*(N*2)/(double)(1024*1024*1024*elapsed); // neglect write out
    printf("%d %g %g %d\n", N, elapsed, BW, (iint)(h_result[0]/Ntests));
    
  }

  return 0;
}
