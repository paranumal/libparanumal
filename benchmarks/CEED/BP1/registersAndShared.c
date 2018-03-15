
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
  //size of shared (in doubles)
    int sh = atoi(argv[2]);
    //size of register array (in doubles) per thread
    int reg = atoi(argv[3]);
    
    int gjNq = N+2;
    int Nq = N;
    int Np = (Nq+1)*(Nq+1)*(Nq+1);
    
  
     int Nelements = 512;
        int Nbytes =(Np*2 +7*gjNp)/2))*Nelements;
  
  double *h_in = (double*) calloc(Nbytes, sizeof(double));
  double *h_out = (double*) calloc(Nbytes, sizeof(double));
  
  
  for(int n=0;n<Nbytes;++n){
    h_a1[n] = (double)rand() / RAND_MAX;
    h_a2[n] = (double)rand() / RAND_MAX;
    
   
  }
 
  char deviceConfig[BUFSIZ];
  sprintf(deviceConfig, "mode = CUDA, deviceID = 0");
  
    occa::device device;
  device.setup(deviceConfig);
	occa::memory o_in = device.malloc(Nbytes*sizeof(double), h_in);
  occa::memory o_out = device.malloc(Nbytes*sizeof(double), h_out);

 
  occa::kernelInfo kernelInfo;
  kernelInfo.addDefine("p_shared", sh);
  kernelInfo.addDefine("p_reg", reg);
  kernelInfo.addDefine("p_Nblocks", (int) Nbytes/1024+1);
  kernelInfo.addDefine("p_Nthreads", 1024);
  
  kernelInfo.addParserFlag("automate-add-barriers", "disabled");
  kernelInfo.addCompilerFlag("-O0");
  
  occa::kernel testSharedKernel
    = device.buildKernelFromSource(DHOLMES "/okl/testSharedRegisters.okl",
				   "testSharedRegisters_v0",
				   kernelInfo);
	
	
    occa::streamTag startTag = device.tagStream();
    
    for(int test=0;test<Ntests;++test){
      testSharedKernel(Nbytes, o_in,o_out);
      o_out.copyTo(h_out);
      
      
    }
    
    occa::streamTag stopTag = device.tagStream();
     double elapsed = device.timeBetween(startTag, stopTag)/(double)Ntests;
     
     printf("time %8.8f \n", elapsed);
  return 0;
}
