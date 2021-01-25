
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

#include "platform.hpp"

occa::memory platform_t::mallocConvert(dlong N,
				       const dfloat *dataIn,
				       const char *typeOutString){

  // if dfloat and typeOut are same type then just malloc and copy
  if(!strcmp(typeOutString, dfloatString)){
    return device.malloc(N*sizeof(dfloat), dataIn);
  }

  printf("Converting array from HOST %s to DEVICE %s\n",
	 dfloatString, typeOutString);
  
  // size of out type
  int sizeOut = sizeof(dfloat);
#if 0
  if(!strcmp(typeOutString,"half"))
    sizeOut = 2; // hack in lieu of sizeof(half)
#endif
  if(!strcmp(typeOutString,"float"))
    sizeOut = sizeof(float);
  if(!strcmp(typeOutString,"double"))
    sizeOut = sizeof(double);
  

  // copy data to device for conversion on device
  occa::memory o_dataIn  = device.malloc(N*sizeof(dfloat), dataIn);
  occa::memory o_dataOut = device.malloc(N*sizeOut);

  char kernelCstr[BUFSIZ];
  
  sprintf(kernelCstr,
	  "@kernel void convert(%s N, %s *dataIn, %s *dataOut){  for(%s n=0;n<N;++n;@tile(256,@outer,@inner)){ dataOut[n] = (%s) dataIn[n]; } }",
	  dlongString, dfloatString, typeOutString, dlongString, typeOutString);

  occa::properties kernelInfo;
  kernelInfo["includes"].asArray();
  
  occa::kernel convertKernel = buildKernelFromString(kernelCstr, "convert", kernelInfo);  

  convertKernel(N, o_dataIn, o_dataOut);

  return o_dataOut;
  
}
