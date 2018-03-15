
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
    = device.buildKernelFromSource(DHOLMES "/okl/scaledAdd.okl",
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
