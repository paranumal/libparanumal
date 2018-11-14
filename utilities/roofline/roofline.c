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


#include "roofline.h"

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);

  if(argc!=3){

    printf("usage: ./roofline executable setupFile \n");
    exit(-1);
  }
  
  char *executable = strdup(argv[1]);
  char *setupFile = strdup(argv[2]);

  char cmd1[BUFSIZ], cmd2[BUFSIZ];
  
  sprintf(cmd1, "nvprof --csv --metrics dram_read_throughput,dram_write_throughput,flop_dp_efficiency %s %s 2> out1", executable, setupFile);
  printf("cmd1 = `%s`\n", cmd1);
  system(cmd1);

  sprintf(cmd2, "nvprof --csv  %s %s 2> out2 ", executable, setupFile);
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
  char *line;

  do{
    line = fgets(buf, BUFSIZ, fp2);
    rest = buf;
    if(line && strstr(rest, "GPU activities") && strstr(rest, "_occa") ){
      strtok(rest, ",");
      for(int n=0;n<4;++n){
	token = strtok(NULL, ",");
      }
      kernelTime[Nkernels] = atof(token)/1000.f;// convert from ms to s

      for(int n=0;n<3;++n){
	token = strtok(NULL, ",");
      }
      kernelNames[Nkernels] = strdup(strtok(token, "\""));

      //      printf("kernel Name = %s, took %lg seconds\n", kernelNames[Nkernels], kernelTime[Nkernels]);
      
      ++Nkernels;
      
    }
  }while(line!=NULL);
  
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
	  sscanf(token, "%lf", &val);

	  if(strstr(buf, "flop_dp_efficiency")) 
	    printf("match on kernel %d, %s consumed %lf %% of peak GFLOPS/s in %g seconds\n", knl, kernelNames[knl], val, kernelTime[knl]);

	  if(strstr(buf, "dram_read_throughput")) 
	    printf("match on kernel %d, %s consumed %lf GB/s of read throughput in %g seconds\n", knl, kernelNames[knl], val, kernelTime[knl]);

	  if(strstr(buf, "dram_write_throughput")) 
	    printf("match on kernel %d, %s consumed %lf GB/s of write throughput in %g seconds\n", knl, kernelNames[knl], val, kernelTime[knl]);

	  break;
	}
      }
    }
    
  }while(line);

  fclose(fp1);
  fclose(fp2);
  
  MPI_Finalize();
  
  return 1;
}
