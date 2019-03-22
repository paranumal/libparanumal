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
	  // N regulates the size of an array
	  int  N = atoi(argv[1]);
    //size of register array (in doubles) per thread
    int reg = atoi(argv[2]);
    
    int gjNq = N+2;
    int Nq = N+1;
    int Np = (N+1)*(N+1)*(N+1);
    int gjNp = (gjNq)*(gjNq)*(gjNq);
   
    
  
     int Nelements = 512;
        int Nbytes =((Np*2 +7*gjNp))*Nelements/2;
         printf("N =%d Nq=%d Np=%d gjNp = %d Nbytes %d\n", N, Nq, Np, gjNp, Nbytes);
  //Nbytes =((sizeof(dfloat)*(mesh->Np*2 +7*gjNp)/2));
  double *h_in = (double*) calloc(Nbytes, sizeof(double));
  double *h_out = (double*) calloc(Nbytes, sizeof(double));
  
  
  for(int n=0;n<Nbytes;++n){
    h_in[n] = (double)rand() / RAND_MAX;
    
   
  }
 
  char deviceConfig[BUFSIZ];
  sprintf(deviceConfig, "mode = CUDA, deviceID = 0");
  
    occa::device device;
  device.setup(deviceConfig);
	occa::memory o_in = device.malloc(Nbytes*sizeof(double), h_in);
  occa::memory o_out = device.malloc(Nbytes*sizeof(double), h_out);

 
  occa::kernelInfo kernelInfo;
  kernelInfo.addDefine("p_reg", reg);
    kernelInfo.addDefine("p_Nq", N+1);
      kernelInfo.addDefine("p_Np", (N+1)*(N+1)*(N+1));
kernelInfo.addDefine("p_gjNp", (Nq+1)*(Nq+1)*(Nq+1));

kernelInfo.addDefine("elPerBlock", Np*2+7*gjNp);
  
  kernelInfo.addParserFlag("automate-add-barriers", "disabled");
  kernelInfo.addCompilerFlag("  --compiler-options -O3");
  
  occa::kernel testSharedKernel
    = device.buildKernelFromSource(LIBP_DIR "/okl/testRegistersAndShared.okl",
				   "testSharedRegisters_v0",
				   kernelInfo);
	
	
    occa::streamTag startTag = device.tagStream();
    
    for(int test=0;test<Ntests;++test){
      testSharedKernel(Nbytes, o_in,o_out);     
    }
    
    occa::streamTag stopTag = device.tagStream();
 /*   o_out.copyTo(h_out);

for (int n=0; n<Nbytes; ++n)
{
	if (h_in[n]!=h_out[n])
	printf("n=%d in %f out %f \n", n, h_in[n], h_out[n]);
}
  */  
     double copyElapsed = device.timeBetween(startTag, stopTag);
     
    Nbytes =(sizeof(double)/2)*(Np*2 +7*gjNp); 
     double copyBandwidth = Nelements*((Nbytes*Ntests*2.)/(1024.*1024.*1024.*copyElapsed));
     
     printf("time %8.8f bw %17.15E \n", copyElapsed, copyBandwidth);
     
}
