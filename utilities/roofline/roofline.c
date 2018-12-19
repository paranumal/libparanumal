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

  ../../utilities/autoTester/autoTester "../../utilities/roofline/roofline ./advectionMain" setups/setupTemplateHex3D.rc

*/

#include <unistd.h>
#include "roofline.h"

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);

  if(argc!=3){

    printf("usage: ./roofline executable setupFile \n");
    exit(-1);
  }
  
  char *executable = strdup(argv[1]);
  char *setupFile = strdup(argv[2]);

  setupAide options(setupFile);
  
  char cmd1[BUFSIZ], cmd2[BUFSIZ];
  
  sprintf(cmd1, "nvprof -u ms --csv --metrics dram_read_throughput,dram_write_throughput,dram_write_transactions,dram_read_transactions,flop_dp_efficiency,flop_count_dp %s %s 2> out1", executable, setupFile);
  printf("cmd1 = `%s`\n", cmd1);
  system(cmd1);

  sprintf(cmd2, "nvprof -u ms --csv  %s %s 2> out2 ", executable, setupFile);
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
    if(line && strstr(rest, "GPU activities") && strstr(rest, "_occa") ){
      strtok(rest, ",");
      for(int n=0;n<4;++n){
	token = strtok(NULL, ",");
      }
      kernelTime[Nkernels] = atof(token)/1.e3;// convert from ns to s

      for(int n=0;n<3;++n){
	token = strtok(NULL, ",");
      }
      kernelNames[Nkernels] = strdup(strtok(token, "\""));

      //      printf("kernel Name = %s, took %lg seconds\n", kernelNames[Nkernels], kernelTime[Nkernels]);
      
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

    if(line && strstr(buf, "_occa")){

      
      int knl;
      for(knl = 0;knl<Nkernels;++knl){
	if(strstr(buf, kernelNames[knl])){
	  
	  rest = strdup(buf);
	  
	  token = strtok(rest, ",");
	  for(int n=0;n<7;++n){
	    token = strtok(NULL, ",");
	  }

	  double val;
	  long long int cnt;
	  
	  if(strstr(buf, "flop_dp_efficiency")){
	    sscanf(token, "%lf", &val);

	    kernelFlopEfficiency[knl] = val;
	  }

	  if(strstr(buf, "dram_read_throughput")){
	    sscanf(token, "%lf", &val);
	    kernelReadThroughput[knl] = val;
	  }

	  if(strstr(buf, "dram_write_throughput")){
	    sscanf(token, "%lf", &val);		  
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
  double maxBWest = 0;
  {

    long long int bytes = 2*1024*1024*(long long int)1024;
    occa::memory o_a = device.malloc(bytes/2);
    occa::memory o_b = device.malloc(bytes/2);

    device.finish();
    
    occa::streamTag start = device.tagStream();

    int Ntests = 10;
    for(int n=0;n<Ntests;++n){
      o_a.copyFrom(o_b);
      o_b.copyFrom(o_a);
    }

    occa::streamTag end = device.tagStream();
    
    device.finish();
    
    double elapsed = device.timeBetween(start, end);
    elapsed /= Ntests;

    maxBWest = 2*bytes/(elapsed*1.e9);
    o_a.free();
    o_b.free();
  }
    

  int knl;
  for(knl = 0;knl<Nkernels;++knl){

    char resultsName[BUFSIZ];
    sprintf(resultsName, "%s.dat", kernelNames[knl]);
    
    FILE *fpResults = fopen(resultsName, "a");
    fprintf(fpResults, "\%\%  arithmeticIntensity, perf, kernelMaxEmpiricalBandwidth, maxEstimatedBandwidth\n");
    
    //    long long int bytes = (kernelReadThroughput[knl]+kernelWriteThroughput[knl])*kernelTime[knl]*1.e9;
    long long int bytes = (kernelBytesRead[knl]+kernelBytesWritten[knl])*32;
    long long int flops = kernelFlopCount[knl];
    double arithmeticIntensity = (double)flops/bytes;
    double perf = (flops/kernelTime[knl])/1.e9; // convert to GFLOPS/s

    printf("perf = %lf, eff = %lf\n",perf, kernelFlopEfficiency[knl]);
    double  maxGFLOPSest = 100*perf/kernelFlopEfficiency[knl]; // since efficiency is given in percent
    
    occa::memory o_a = device.malloc(bytes/2);
    occa::memory o_b = device.malloc(bytes/2);

    occa::streamTag start = device.tagStream();

    int Ntests = 10;
    for(int n=0;n<Ntests;++n)
      o_a.copyFrom(o_b);

    occa::streamTag end = device.tagStream();
    
    device.finish();

    double elapsed = device.timeBetween(start, end);
    elapsed /= Ntests;
    
    kernelMaxEmpiricalBandwidth[knl] = (bytes/elapsed)/1.e9; // convert max empirical bw for this vector size to GB/s
    
    fprintf(fpResults, "%lg, %lg, %lg, %lg, %lg\n", arithmeticIntensity, perf, kernelMaxEmpiricalBandwidth[knl],
	    maxGFLOPSest, maxBWest);

    o_a.free();
    o_b.free();
    fflush(fpResults);
    fclose(fpResults);
  }

#if 0
  FILE *fpMatlab = fopen("rooflineScript.m", "w");

  for(knl = 0;knl<Nkernels;++knl){
    fprintf(fpMatlab, "%s = []\n", kernelNames[knl]);
  }
  fprintf(fpMatlab, "\n rooflineResults;\n");

  printf("maxGFLOPSest=%g\n", maxGFLOPSest);
  
  for(knl = 0;knl<Nkernels;++knl){
    fprintf(fpMatlab, "figure(%d);\n", knl+1);
    fprintf(fpMatlab, "scatter(%s(:,1),%s(:,2), 'filled', 'ko');\n", kernelNames[knl], kernelNames[knl]);

    fprintf(fpMatlab, "hold on;\n");
    fprintf(fpMatlab, "plot([0,20], [0,%lg]);", maxBWest*20);
    fprintf(fpMatlab, "plot([0,40], [%lg,%lg]);", maxGFLOPSest, maxGFLOPSest);
    
    fprintf(fpMatlab, "set(gca, 'FontSize', 14);\n");
    fprintf(fpMatlab, "xlabel('Arithmetic Intensity (FP64 flops)/deviceMemoryBytes', 'FontSize', 14);\n");
    fprintf(fpMatlab, "ylabel('FP64 GFLOPS/s', 'FontSize', 14);\n");
    fprintf(fpMatlab, "title(\"%s\", 'FontSize', 16, 'Interpreter', 'None');\n", kernelNames[knl]);

    // superimpose roofline
    fprintf(fpMatlab, "hold on; \n");
    fprintf(fpMatlab, "plot(%s(:,1),%s(:,1).*%s(:,3),'r*', 'LineWidth', 2); \n",
	    kernelNames[knl], kernelNames[knl], kernelNames[knl]);
    fprintf(fpMatlab, "axis([0,ceil(max(%s(:,1))),0,max(max(%s(:,1).*%s(:,3),%s(:,2)))]);\n",
	    kernelNames[knl],kernelNames[knl], kernelNames[knl], kernelNames[knl]);
  }


  
  fclose(fpMatlab);
#endif
  
  MPI_Finalize();
  
  return 0;
}
