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

/*
  example usage: in advection in combination with autoTester:

  rm bandwidthTesterResults.m bandwidthTesterScript.m

  ../../utilities/autoTester/autoTester "../../utilities/bandwidthTester/bandwidthTester ./advectionMain" setups/setupTemplateHex3D.rc

*/

#include <unistd.h>
#include "bandwidthTester.h"

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);

  if(argc!=2){

    printf("usage: ./bandwidthTester setupFile \n");
    exit(-1);
  }
  
  char *setupFile = strdup(argv[1]);

  setupAide options(setupFile);
  
  // now benchmark memory on device
  occa::device device;

  char deviceConfig[BUFSIZ];

  int device_id = 0;

  options.getArgs("DEVICE NUMBER" ,device_id);

  // read thread model/device/platform from options
  if(options.compareArgs("THREAD MODEL", "CUDA")){
    sprintf(deviceConfig, "mode: 'CUDA', device_id: %d",device_id);
  }
  else if(options.compareArgs("THREAD MODEL", "HIP")){
    sprintf(deviceConfig, "mode: 'HIP', device_id: %d",device_id);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenCL")){
    int plat;
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "mode: 'OpenCL', device_id: %d, platform_id: %d", device_id, plat);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenMP")){
    sprintf(deviceConfig, "mode: 'OpenMP' ");
  }
  else{
    sprintf(deviceConfig, "mode: 'Serial' ");
  }

  device.setup(deviceConfig);

  // profile big copy

  long long int minDataBytes = 1024, maxDataBytes = 1024000000, stepDataBytes = 1024;
  
  options.getArgs("MINIMUM DATA BYTES", minDataBytes);
  options.getArgs("MAXIMUM DATA BYTES", maxDataBytes);
  options.getArgs("STEP DATA BYTES",   stepDataBytes);

  long long int bytes = maxDataBytes;
  occa::memory o_a = device.malloc(bytes);
  occa::memory o_b = device.malloc(bytes);
  occa::memory o_c = device.malloc(bytes);

  printf("\%\% BYTES COPIED, ELAPSED TIME (s), BANDWIDTH (GB/s) \n");
  printf("memcpyBW = [\n");
  int Ntests = 10;
  
  for(long long int streamedBytes=minDataBytes;streamedBytes<maxDataBytes;streamedBytes+=stepDataBytes){

    device.finish();

    occa::streamTag start = device.tagStream();
    
    for(int n=0;n<Ntests;++n){
      o_a.copyFrom(o_b, streamedBytes/2, 0);
      o_b.copyFrom(o_a, streamedBytes/2, 0);
    }
    
    occa::streamTag end = device.tagStream();
    
    device.finish();
    
    double elapsed = device.timeBetween(start, end);
    elapsed /= (2*Ntests);
    
    double bw = (streamedBytes/elapsed)/1.e9;
    
    printf("%lld %lg %lg\n", streamedBytes, elapsed, bw);
  }  

  printf("];\n");
  
#define dfloat double

  occa::properties props;
  props["defines/dfloat"] = (sizeof(dfloat)==4) ? "float": "double";

  char fileName[BUFSIZ];
  sprintf(fileName, "bandwidthTester.okl");
  
  // test throughput from kernels
  int Nkernels = 3;
  for(int knl=0;knl<Nkernels;++knl){

    char kernelName[BUFSIZ];
    sprintf(kernelName, "bandwidthTesterK%02d", knl);
    
    occa::kernel kernel = device.buildKernel(fileName, kernelName, props);

    printf("memcpyKNL%02d = [\n", knl);
    
    for(long long int streamedBytes=minDataBytes;streamedBytes<maxDataBytes;streamedBytes+=stepDataBytes){
      int Ndfloats = 1;
      int bytes = 1;
      
      switch(knl){
      case 0: // write
      case 2: // read
	Ndfloats = streamedBytes/sizeof(dfloat);
	bytes = Ndfloats*sizeof(dfloat);
	break;
      case 1: // read + write
	Ndfloats = streamedBytes/(2*sizeof(dfloat));
	bytes = 2*Ndfloats*sizeof(dfloat);
	break;
      }
      
      device.finish();
      
      occa::streamTag start = device.tagStream();
      
      for(int n=0;n<Ntests;++n){
	kernel(Ndfloats, o_a, o_b, o_c);
      }

      occa::streamTag end = device.tagStream();
      
      device.finish();
      
      double elapsed = device.timeBetween(start, end);
      elapsed /= Ntests;
      
      double bw = (bytes/elapsed)/1.e9;
      
      printf("%lld %lg %lg\n", bytes, elapsed, bw);
      
    }

    printf("];\n");
    
  }
  
  
  o_a.free();
  o_b.free();
  o_c.free();
    
  MPI_Finalize();
  
  return 0;
}
