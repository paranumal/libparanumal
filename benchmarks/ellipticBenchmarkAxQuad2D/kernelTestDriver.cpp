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
#include <math.h>

#define dlong int
#define dlongString "int"

#if 1
#define dfloat double
#define dfloatString "double"
#else
#define dfloat float
#define dfloatString "float"
#endif

void randCalloc(occa::device &device, int sz, dfloat **pt, occa::memory &o_pt){

  *pt = (dfloat*) calloc(sz, sizeof(dfloat));

  dfloat sum = 0;
  for(int n=0;n<sz;++n){
    (*pt)[n] = drand48()-0.5;
    sum += pow((*pt)[n],2);
  }

  o_pt = device.malloc(sz*sizeof(dfloat), *pt);

}

int main(int argc, char **argv){

  int NKernels = 1;

  // default to 512 elements if no arg is given
  dlong E = (argc>=2) ? atoi(argv[1]):512;
  int p_N = (argc>=3) ? atoi(argv[2]):5;
  int option = (argc>=4) ? atoi(argv[3]):1;  
  int p_Ne = (argc>=5) ? atoi(argv[4]):1;
  int p_Nb = (argc>=6) ? atoi(argv[5]):1;

  int kMin = (argc>=7) ? atoi(argv[6]):0;
  int kMax = (argc>=8) ? atoi(argv[7]):8;

  int p_Np = (p_N+1)*(p_N+1);
  int p_Nfp = (p_N+1);
  int p_Nfaces = 4;
  int p_NfacesNfp = p_Nfaces*p_Nfp;
  int p_Nq = p_N+1;
  int p_Nggeo = 7;
  
  printf("==============================================================\n");
  printf("===================== BASIC INFO =============================\n");
  printf("==============================================================\n");
  printf("Number of elements : %d\n", E);
  printf("Polynomial degree  : %d\n", p_N);
  printf("Nodes per element  : %d\n", p_Np);
  printf("Elements per block : %d\n", p_Ne);
  printf("Outputs per thread : %d\n", p_Nb);
  printf("Running kernels    : %d to %d \n", kMin, kMax);
  printf("==============================================================\n");
  printf("==============================================================\n");
  printf("==============================================================\n");
  printf("\n\n");
  int BSIZE  = p_Np;

  // number of geometric factors
  int Nvgeo =  p_Nggeo;
  int Niter = 10, it;
  double gflops;
  
  gflops = p_Np*(4*p_Nq + 8 + 4*p_Nq); // FIX LATER
  gflops *= Niter;

  // build some dummy storage & parameters
  double results[15];
  double roofline[15];
  double timeData[15];
  occa::device device;
  occa::kernel ellipticAxKernel[15];
  occa::kernel correctRes;
  occa::kernelInfo kernelInfo;

  // specify device configuration

  // device.setup("mode = Serial");
  // device.setup("mode = OpenMP  , schedule = compact, chunk = 10");
  // device.setup("mode = OpenCL  , platformID = 0, deviceID = 0");
  device.setup("mode = CUDA    , deviceID = 0");
  // device.setup("mode = Pthreads, threadCount = 4, schedule = compact, pinnedCores = [0, 0, 1, 1]");
  // device.setup("mode = COI     , deviceID = 0");

  // compiler variables to be passed to backend compiler by OCCA
  kernelInfo.addDefine("dfloat", dfloatString);
  kernelInfo.addDefine("dlong", dlongString);
  kernelInfo.addDefine("p_N", p_N);
  kernelInfo.addDefine("p_Nq", p_Nq);
  kernelInfo.addDefine("p_Np", p_Np);
  kernelInfo.addDefine("p_Nfp", p_Nfp);
  kernelInfo.addDefine("p_Nfaces", p_Nfaces);
  kernelInfo.addDefine("p_NfacesNfp", p_Nfaces*p_Nfp);
  kernelInfo.addDefine("BSIZE", (int)BSIZE);
  kernelInfo.addDefine("p_Nggeo", p_Nggeo);
  kernelInfo.addDefine("p_Ne", p_Ne);
  kernelInfo.addDefine("p_Nb", p_Nb);

  dlong elementOffset = 0;
  dfloat    *ggeo, *DT,   *q, *Aq;
  occa::memory o_ggeo, o_DT, o_q, o_Aq; 
  
  srand(12345);
  
  kernelInfo.addDefine("p_G00ID", 0);
  kernelInfo.addDefine("p_G01ID", 1);
  kernelInfo.addDefine("p_G02ID", 2);
  kernelInfo.addDefine("p_G11ID", 3);
  kernelInfo.addDefine("p_G12ID", 4);
  kernelInfo.addDefine("p_G22ID", 5);
  kernelInfo.addDefine("p_GWJID", 6);

  randCalloc(device, E*p_Nggeo*BSIZE*BSIZE, &ggeo, o_ggeo);
  randCalloc(device, BSIZE*BSIZE, &DT, o_DT);
  randCalloc(device, BSIZE*E, &q, o_q);
  randCalloc(device, BSIZE*E, &Aq, o_Aq);

  char buf[200];
  for (int k =kMin; k<kMax+1; k++){
    printf("compiling Quad2D kernel %d ...\n", k);
    sprintf(buf, "ellipticPartialAxQuad2D_Ref%d", k); 
    ellipticAxKernel[k] = device.buildKernelFromSource("ellipticAxQuad2D.okl", buf, kernelInfo);
  }
  
  occa::initTimer(device);

  // queue Ax kernels
  for (int k =kMin;k<=kMax; k++){
    dfloat lambda = 1.;
    occa::streamTag startTag = device.tagStream();

    for(it=0;it<Niter;++it){
      ellipticAxKernel[k](E,
			  elementOffset,
			  o_ggeo,
			  o_DT,
			  lambda,
			  o_q,       
			  o_Aq);
    }
    
    //adjust GFLOPS for Ref2 and Ref3
    occa::streamTag stopTag = device.tagStream();
    double elapsed = device.timeBetween(startTag, stopTag);

    gflops = p_Np*(4*p_Nq + 8 + 4*p_Nq); // FIX LATER
    gflops *=Niter;      

    printf("\n\nKERNEL %d  ================================================== \n\n", k);
    printf("OCCA elapsed time = %g\n", elapsed);
    timeData[k] = elapsed/Niter;
    printf("number of flops = %f time = %f \n", gflops, elapsed);
    results[k] =E*gflops/(elapsed*1000*1000*1000); 
    //elapsed/Niter;
    //
    printf("OCCA: estimated time = %17.15f gflops = %17.17f\n", results[k], E*gflops/(elapsed*1000*1000*1000));
    printf("GFL %17.17f \n",E*gflops/(elapsed*1000*1000*1000) );      
    // compute l2 of data
    o_Aq.copyTo(Aq);
    dfloat normAq = 0;

    for(int n=0;n<E*BSIZE;++n)
      normAq += Aq[n]*Aq[n];
    normAq = sqrt(normAq);

    printf("OCCA: normAq = %17.15lf\n", normAq);
    
    for(int n=0;n<E*BSIZE;++n){
      Aq[n] = 0.0f;
    }
    
    o_Aq.copyFrom(Aq);
    
  }//for
  
  
  //printf("\n\nBWfromCopy%d = [", E);
  for (int k=kMin; k<=kMax; k++){
    printf("==== this is kernel %d \n", k);
    int p_Nq = k+1;
    int p_gjNq = k+2;
    int Niter = 10;
    double gflops;
    long long int Nbytes;

    Nbytes = p_Nq*p_Nq*sizeof(dfloat) + E*p_Np*p_Nggeo*sizeof(dfloat)+E*sizeof(int)+E*p_Np*2*sizeof(dfloat);

    
    Nbytes /= 2;
    gflops = p_Np*(4*p_Nq + 8 + 4*p_Nq); // FIX LATER

    occa::memory o_foo = device.malloc(Nbytes);
    occa::memory o_bah = device.malloc(Nbytes);

    occa::streamTag startCopy = device.tagStream();
    for(int it=0;it<Niter;++it){
      o_bah.copyTo(o_foo);
    }
    occa::streamTag endCopy = device.tagStream();
    double copyElapsed = device.timeBetween(startCopy, endCopy);
    double copyBandwidth = ((Nbytes*Niter*2.)/(1000.*1000.*1000.*copyElapsed));
    printf("copytime = %16.17f \n", copyElapsed);
    printf("pNp = %d p_Nfp = %d Nbytes = %d \n", p_Np, p_Nfp, Nbytes);    
    printf("copy BW = %f gflops = %f bytes = %d \n", copyBandwidth, gflops, Nbytes);
    roofline[k] = (copyBandwidth*gflops*E)/(2*Nbytes);

    o_foo.free();
    o_bah.free();
  }

  //printf("];\n\n");
  
  printf("\n\nROOFLINE = [");
  for (int k=kMin; k<=kMax; k++){
    
    printf(" %16.17f ", roofline[k]);
  }
  
  printf("]\n\n");


  printf("\n\nResults(:,%d)  = [", p_N);
  for (int k=kMin; k<=kMax; k++){

    printf(" %16.17f ", results[k]);
  }

  printf("];\n\n");

  printf("\n\ntimeData(:, %d) = [",p_N);
  for (int k=kMin; k<=kMax; k++){

    printf(" %16.17f ", timeData[k]);
  }

  printf("];\n\n");

  printf("\n\ngigaNodes(:, %d) = [",p_N);
  for (int k=kMin; k<=kMax; k++){

    printf(" %16.17f ", (p_Np*E/10e9)/timeData[k]);
  }

  printf("];\n\n");


  exit(0);
  return 0;

}
