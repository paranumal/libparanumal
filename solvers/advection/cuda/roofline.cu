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

  rm rooflineResults.m rooflineScript.m

  ../../utilities/autoTester/autoTester "../../utilities/roofline/roofline ./advectionMain" 

*/

#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cuda.h"
#define mymax(a,b) ((a>b)?(a):(b))

__global__ void copyKernel(const int N, const float * __restrict__ a, float * __restrict__ b){

  int n= threadIdx.x + blockIdx.x*blockDim.x;

  if(n<N)
    b[n] = a[n];

}

int main(int argc, char **argv){

  if(argc!=2){

    printf("usage: ./roofline executable \n");
    exit(-1);
  }
  
  char *executable = strdup(argv[1]);

  char cmd1[BUFSIZ], cmd2[BUFSIZ];
  
  sprintf(cmd1, "nvprof -u ms --csv --metrics dram_read_throughput,dram_write_throughput,dram_write_transactions,dram_read_transactions,flop_dp_efficiency,flop_count_dp %s 2> out1", executable);
  printf("cmd1 = `%s`\n", cmd1);
  system(cmd1);

  sprintf(cmd2, "nvprof -u ms --csv  %s  2> out2 ", executable);
  printf("cmd2 = `%s`\n", cmd2);
  system(cmd2);

  char buf[BUFSIZ];
  char *token;
  char *rest = buf;
  
  FILE *fp2 = fopen("out2", "r");
  
  int Nkernels = 0;
  int maxNkernels = 1000;
  char **kernelNames = (char**) calloc(maxNkernels, sizeof(char*));
  double *kernelTime = (double*) calloc(maxNkernels, sizeof(double));
  long long int *kernelFlopCount = (long long int*) calloc(maxNkernels, sizeof(long long int));
  double *kernelFlopEfficiency = (double*) calloc(maxNkernels, sizeof(double));
  double *kernelReadThroughput = (double*) calloc(maxNkernels, sizeof(double));
  double *kernelWriteThroughput = (double*) calloc(maxNkernels, sizeof(double));
  long long int *kernelBytesRead = (long long int*) calloc(maxNkernels, sizeof(long long int));
  long long int *kernelBytesWritten = (long long int*) calloc(maxNkernels, sizeof(long long int));
  double *kernelMaxEmpiricalBandwidth = (double*) calloc(maxNkernels, sizeof(double));
  char *line;

  do{
    line = fgets(buf, BUFSIZ, fp2);
    rest = buf;
    if(line && strstr(rest, "GPU activities") && strstr(rest, "advection") ){ // AHEM
      strtok(rest, ",");
      for(int n=0;n<6;++n){
	token = strtok(NULL, ",");
      }
      kernelTime[Nkernels] = atof(token)/1.e3;// convert from ns to s


      token = strtok(NULL, "\"");
      kernelNames[Nkernels] = strdup(strtok(token, "("));

      printf("kernel Name = %s, took %lg seconds\n", kernelNames[Nkernels], kernelTime[Nkernels]);
      
      ++Nkernels;
      
    }
  }while(line!=NULL);

  typedef struct{
    double dram_read_throughput;
    double dram_write_throughput;
    double flop_dp_efficiency;
  } performance;
  
  FILE *fp1 = fopen("out1", "r");

  do{
    fgets(buf, BUFSIZ, fp1);
  }while(!strstr(buf, "Invocations")); // slight dangerous

  do{
    
    // assume "Device","Kernel","Invocations","Metric Name","Metric Description","Min","Max","Avg"
    line = fgets(buf, BUFSIZ, fp1);

    if(line && strstr(buf, "advection")){//  AHEM

      
      int knl;
      for(knl = 0;knl<Nkernels;++knl){
	if(strstr(buf, kernelNames[knl])){
	  //	  printf("buf=%s\n", buf);
	  
	  rest = strdup(buf);
	  
	  token = strtok(rest, "\"");
	  for(int n=0;n<6;++n){
	    token = strtok(NULL, "\"");
	    //	    printf("token=%s\n", token);
	  }
	  token = strtok(NULL, ",");
	  token = strtok(NULL, ",");
	  token = strtok(NULL, ",");
	  //	  printf("token#=%s\n", token);
	  
	  double val;
	  long long int cnt;
	  
	  if(strstr(buf, "flop_dp_efficiency")){
	    sscanf(token, "%lf", &val);

	    kernelFlopEfficiency[knl] = val;
	  }

	  if(strstr(buf, "dram_read_throughput")){
	    sscanf(token, "%lf", &val);
	    //	    printf("dram_read %lf\n", val);
	    kernelReadThroughput[knl] = val;
	  }

	  if(strstr(buf, "dram_write_throughput")){
	    sscanf(token, "%lf", &val);
	    //	    printf("dram_write %lf\n", val);
	    kernelWriteThroughput[knl] = val;
	  }

	  if(strstr(buf, "flop_count_dp")){
	    sscanf(token, "%lld", &cnt);
	    kernelFlopCount[knl] = cnt;
	  }

	  if(strstr(buf, "dram_read_transactions")){
	    sscanf(token, "%lld", &cnt);
	    kernelBytesRead[knl] = cnt;
	  }

	  if(strstr(buf, "dram_write_transactions")){
	    sscanf(token, "%lld", &cnt);
	    kernelBytesWritten[knl] = cnt;
	  }
	  
	  break;
	}
      }
    }
    
  }while(line);

  fclose(fp1);
  fclose(fp2);


  // now benchmark memory on device
  // profile big copy
  double maxBWest = 0;
  {

    long long int bytes = 2*1024*1024*(long long int)1024;

    void *o_a, *o_b;

    cudaMalloc(&o_a, bytes/2);
    cudaMalloc(&o_b, bytes/2);

    cudaDeviceSynchronize();

    cudaEvent_t start, end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);

    cudaEventRecord(start);    

    dim3 B(256,1,1);
    dim3 G( ((bytes/8) + 255)/256, 1, 1);
    
    int Ntests = 1;
    for(int n=0;n<Ntests;++n){
#if 0
      cudaMemcpy(o_a, o_b, (bytes/2), cudaMemcpyDeviceToDevice);
      cudaMemcpy(o_b, o_a, (bytes/2), cudaMemcpyDeviceToDevice);
#endif

      copyKernel <<< G, B >>> (bytes/8, (float*) o_a, (float*) o_b);
      copyKernel <<< G, B >>> (bytes/8, (float*) o_b, (float*) o_a);
      
    }

    cudaEventRecord(end);
    cudaEventSynchronize(end);

    cudaDeviceSynchronize();

    float elapsed;
    cudaEventElapsedTime(&elapsed, start, end);
    elapsed /= (Ntests*1000.);

    maxBWest = 2*bytes/(elapsed*1.e9);
    cudaFree(&o_a);
    cudaFree(&o_b);
  }
    

  int knl;
  for(knl = 0;knl<Nkernels;++knl){

    char resultsName[BUFSIZ];
    sprintf(resultsName, "%s.dat", kernelNames[knl]);
    
    FILE *fpResults = fopen(resultsName, "a");
    fprintf(fpResults, "%%%%  arithmeticIntensity, perf, kernelMaxEmpiricalBandwidth, maxEstimateGflops, maxEstimatedBandwidth\n");
    
    //    long long int bytes = (kernelReadThroughput[knl]+kernelWriteThroughput[knl])*kernelTime[knl]*1.e9;
    long long int bytes = (kernelBytesRead[knl]+kernelBytesWritten[knl])*32;
    long long int flops = kernelFlopCount[knl];
    double arithmeticIntensity = (double)flops/bytes;
    double perf = (flops/kernelTime[knl])/1.e9; // convert to GFLOPS/s

    printf("perf = %lf, eff = %lf\n",perf, kernelFlopEfficiency[knl]);
    double  maxGFLOPSest = 100*perf/kernelFlopEfficiency[knl]; // since efficiency is given in percent

    void *o_a, *o_b;

    cudaMalloc(&o_a, bytes/2);
    cudaMalloc(&o_b, bytes/2);

    cudaDeviceSynchronize();

    cudaEvent_t start, end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);

    cudaEventRecord(start);    


    dim3 B(256,1,1);
    dim3 G( ((bytes/8) + 255)/256, 1, 1);
    
    int Ntests = 1;
    for(int n=0;n<Ntests;++n){
#if 0
      cudaMemcpy(o_a, o_b, (bytes/2), cudaMemcpyDeviceToDevice);
      cudaMemcpy(o_b, o_a, (bytes/2), cudaMemcpyDeviceToDevice);
#endif
      
      copyKernel <<< G, B >>> (bytes/8, (float*) o_a, (float*) o_b);
      copyKernel <<< G, B >>> (bytes/8, (float*) o_b, (float*) o_a);
    }

    cudaEventRecord(end);
    cudaEventSynchronize(end);

    cudaDeviceSynchronize();

    float elapsed;
    cudaEventElapsedTime(&elapsed, start, end);
    elapsed /= (Ntests*1000.);

//    maxBWest = 2*bytes/(elapsed*1.e9);
    cudaFree(&o_a);
    cudaFree(&o_b);
    
    kernelMaxEmpiricalBandwidth[knl] = (2.*bytes/elapsed)/1.e9; // convert max empirical bw for this vector size to GB/s
    
    fprintf(fpResults, "%lg, %lg, %lg, %lg, %lg\n", arithmeticIntensity, perf, kernelMaxEmpiricalBandwidth[knl],
	    maxGFLOPSest, maxBWest);

    fflush(fpResults);
    fclose(fpResults);
  }

  
  return 0;
}
