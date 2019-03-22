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

#define int int
#define intString "int"

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
    = device.buildKernelFromSource(LIBP_DIR "/okl/testFloats.okl",
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
    //   printf("%d %g %g %d\n", N, elapsed, BW, (int)(h_result[0]/Ntests));
    
    
    
//}
  //  o_adotb.copyTo(h_result);
    
   // double elapsed = device.timeBetween(startTag, stopTag)/(double)Ntests;
  //  dfloat BW = sizeof(dfloat)*(N*2)/(double)(1024*1024*1024*elapsed); // neglect write out
  //  printf("%d %g %g %d\n", N, elapsed, BW, (int)(h_result[0]/Ntests));
    
  //}

  return 0;
}
