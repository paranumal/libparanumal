
#include <stdio.h>
#include <stdlib.h>

#include "occa.hpp"

#define iint int
#define iintString "int"

#if 1
#define dfloat double
#define dfloatString "double"
#else
#define dfloat float
#define dfloatString "float"
#endif

int main(int argc, char **argv){

 
 
  
  

 
  

  int Ntests = 10;
  
    int sz = atoi(argv[1]);
    
//for (int sz = 60000; sz<61000; sz+=1000)
//{
	//first, define memory
	printf("===================== size is %d  ================================\n", sz);
	
	 int maxNf = 1024*sz*3;
  int maxNd = 1024*sz;

//  int maxNblocks = (maxN+blockSize-1)/blockSize;
  
  double *h_a1 = (double*) calloc(maxNd, sizeof(double));
  double*h_a2 = (double*) calloc(maxNd, sizeof(double));
  
  float *h_b1 = (float*) calloc(maxNf, sizeof(float));
  float *h_b2 = (float*) calloc(maxNf, sizeof(float));
  
  
  double *h_oa1 = (double*) calloc(maxNd, sizeof(double));
  
  float *h_ob1 = (float*) calloc(maxNf, sizeof(float));
   printf("maxNd = %d maxfD = %d\n",maxNd, maxNf);
 
  for(int n=0;n<maxNd;++n){
    h_a1[n] = (double)rand() / RAND_MAX;
    h_a2[n] = (double)rand() / RAND_MAX;
    
   
  }
  
    for(int n=0;n<maxNf;++n){
    h_b1[n] = (float)rand() / RAND_MAX;
    h_b2[n] = (float)rand() / RAND_MAX;
  }
   /* DEVICE INFO */
  char deviceConfig[BUFSIZ];
  sprintf(deviceConfig, "mode = CUDA, deviceID = 0");
  
    occa::device device;
  device.setup(deviceConfig);
	 occa::memory o_a1 = device.malloc(maxNd*sizeof(double), h_a1);
  occa::memory o_a2 = device.malloc(maxNd*sizeof(double), h_a2);
  
  occa::memory o_b1 = device.malloc(maxNf*sizeof(float), h_b1);
  occa::memory o_b2 = device.malloc(maxNf*sizeof(float), h_b2);
  
  occa::memory o_oa1 = device.malloc(maxNd*sizeof(double), h_oa1);
  occa::memory o_ob1 = device.malloc(maxNf*sizeof(float), h_ob1);

  /* KERNEL INFO */
  occa::kernelInfo kernelInfo;
  kernelInfo.addDefine("p_exDoubles", sz);
  kernelInfo.addDefine("p_exFloats", sz*3);
  kernelInfo.addDefine("p_exFloatsLess", sz*2);
  kernelInfo.addDefine("p_exBlocksize", 1024);
  
  kernelInfo.addParserFlag("automate-add-barriers", "disabled");
  
  /* KERNEL BUILD */
  occa::kernel floatAndDoubleKernel
    = device.buildKernelFromSource(DHOLMES "/okl/testFloats.okl",
				   "testFloats_v0",
				   kernelInfo);
	
	
    occa::streamTag startTag = device.tagStream();
    
    for(int test=0;test<Ntests;++test){
      floatAndDoubleKernel(o_a1,o_a2, o_b1, o_b2, o_oa1, o_ob1);
      o_oa1.copyTo(h_oa1);
        o_ob1.copyTo(h_ob1);
      
    }
    
    occa::streamTag stopTag = device.tagStream();
     double elapsed = device.timeBetween(startTag, stopTag)/(double)Ntests;
     
     //  dfloat BW = sizeof(dfloat)*(N*2)/(double)(1024*1024*1024*elapsed); // neglect write out
    //   printf("%d %g %g %d\n", N, elapsed, BW, (iint)(h_result[0]/Ntests));
    
    
    
//}
  //  o_adotb.copyTo(h_result);
    
   // double elapsed = device.timeBetween(startTag, stopTag)/(double)Ntests;
  //  dfloat BW = sizeof(dfloat)*(N*2)/(double)(1024*1024*1024*elapsed); // neglect write out
  //  printf("%d %g %g %d\n", N, elapsed, BW, (iint)(h_result[0]/Ntests));
    
  //}

  return 0;
}
