
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
	  // N regulates the size of an array
	  int  N = atoi(argv[1]);
    //size of register array (in doubles) per thread
    int reg = atoi(argv[2]);
    
    int gjNq = N+2;
    int Nq = N+1;
    int Np = (N+1)*(N+1)*(N+1);
    int gjNp = (gjNq)*(gjNq)*(gjNq);
   
    
  
     int Nelements = 512;
        iint Nbytes =((Np*2 +7*gjNp))*Nelements/2;
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
    = device.buildKernelFromSource(DHOLMES "/okl/testRegistersAndShared.okl",
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
