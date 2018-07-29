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

#if 1
#define datafloat double
#define datafloatString "double"
#else
#define datafloat float
#define datafloatString "float"
#endif

void randCalloc(occa::device &device, int sz, datafloat **pt, occa::memory &o_pt){

  *pt = (datafloat*) calloc(sz, sizeof(datafloat));

  datafloat sum = 0;
  for(int n=0;n<sz;++n){
    (*pt)[n] = drand48()-0.5;
    sum += pow((*pt)[n],2);
  }

  o_pt = device.malloc(sz*sizeof(datafloat), *pt);

}

int main(int argc, char **argv){

  int NKernels = 2;

  // default to 512 elements if no arg is given
  int E = (argc>=2) ? atoi(argv[1]):512;
  int p_N = (argc>=3) ? atoi(argv[2]):5;
  int p_Ne = (argc>=4) ? atoi(argv[3]):1;
  int p_Nb = (argc>=5) ? atoi(argv[4]):1;

  int kMin = (argc>=6) ? atoi(argv[5]):0;
  int kMax = (argc>=7) ? atoi(argv[6]):1;


  int p_Np = ((p_N+1)*(p_N+2)*(p_N+3))/6;
  int p_Nfp = ((p_N+1)*(p_N+2))/2;
  int p_Nfaces = 4;
  int p_NfacesNfp = p_Nfaces*p_Nfp;

  printf("==============================================================\n");
  printf("===================== BASIC INFO =============================\n");
  printf("==============================================================\n");
  printf("Number of elements : %d\n", E);
  printf("Polynomial degree  : %d\n", p_N);
  printf("Nodes per element  : %d\n", p_Np);
  printf("Elements per block : %d\n", p_Nb);
  printf("Outputs per thread : %d\n", p_Ne);
  printf("Running kernels    : %d to %d \n", kMin, kMax);
  printf("==============================================================\n");
  printf("==============================================================\n");
  printf("==============================================================\n");
  printf("\n\n");
  int BSIZE  = p_Np;


  // number of geometric factors
  int Nvgeo =  10;
  int Nsgeo = 6;

  int pad;

  int Niter = 10, it;
  double gflops;
  gflops = p_Np*p_Np*2;

  gflops *= Niter;

  // build some dummy storage & parameters
  double results3D[15];
  double roofline[15];
  double timeData[15];
  occa::device device;
  occa::kernel Tet3Dkernel[15];
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
  kernelInfo.addDefine("datafloat", datafloatString);
  kernelInfo.addDefine("p_N", p_N);
  kernelInfo.addDefine("p_Np", p_Np);
  kernelInfo.addDefine("p_Nfp", p_Nfp);
  kernelInfo.addDefine("p_Nfaces", p_Nfaces);
  kernelInfo.addDefine("p_NfacesNfp", p_Nfaces*p_Nfp);
  kernelInfo.addDefine("BSIZE", (int)BSIZE);
  kernelInfo.addDefine("p_Nvgeo", 11);
  kernelInfo.addDefine("p_Nsgeo", 7);
  kernelInfo.addDefine("p_Nggeo", 7);
  kernelInfo.addDefine("p_Ne", p_Ne);
  kernelInfo.addDefine("p_Nb", p_Nb);


  /*
   *void ellipticPartialPreassembledAxTet3D_Ref0(const int Nelements,
   const int elementOffset,
   const datafloat * restrict ggeo,
   const datafloat * restrict SS,
   const datafloat  * restrict q,
   datafloat * restrict Sq)
   *
   * */


  datafloat  *ggeo, *SS,  *q, *Aq;
  occa::memory o_ggeo, o_SS, o_q, o_Aq; 

  srand(12345);

  kernelInfo.addDefine("p_G00ID", 0);
  kernelInfo.addDefine("p_G01ID", 1);
  kernelInfo.addDefine("p_G02ID", 2);
  kernelInfo.addDefine("p_G11ID", 3);
  kernelInfo.addDefine("p_G12ID", 4);
  kernelInfo.addDefine("p_G22ID", 5);
  kernelInfo.addDefine("p_GWJID", 6);
  randCalloc(device, E*7, &ggeo, o_ggeo);
  randCalloc(device, BSIZE*BSIZE*E, &SS, o_SS);
  randCalloc(device, BSIZE*E, &q, o_q);
  randCalloc(device, BSIZE*E, &Aq, o_Aq);

  char buf[200];
  for (int i =kMin; i<kMax+1; i++){
    printf("compiling preassembled matrix Tet kernel %d ...\n", i);
    sprintf(buf, "ellipticPartialPreassembledAxTet3D_Ref%d", i); 
    Tet3Dkernel[i] = device.buildKernelFromSource("ellipticPreassembledAxTet3D.okl", buf, kernelInfo);
  }
  int elementOffset = 0;
  occa::initTimer(device);
  occa::streamTag startTag, stopTag;

  // queue Ax kernels
  for (int i =kMin;i<=kMax; i++){

    device.finish();
    startTag = device.tagStream();
    
    // launch kernel
    for(it=0;it<Niter;++it){
      Tet3Dkernel[i](E,
		     elementOffset,
		     o_SS,
		     o_q,
		     o_Aq);
    }//for
  

    stopTag = device.tagStream();
    device.finish();
    double elapsed = device.timeBetween(startTag, stopTag);
    printf("\n\nKERNEL %d  ================================================== \n\n", i);

    printf("OCCA elapsed time = %g\n", elapsed);
    timeData[i] = elapsed/Niter;
    results3D[i] =E*gflops/(elapsed*1.e9); 

    printf("OCCA: estimated time = %17.15f gflops = %17.17f\n", results3D[i], E*gflops/(elapsed*1.e9));
    printf("GFL %17.17f \n",E*gflops/(elapsed*1.e9) );      
    // compute l2 of data
    o_Aq.copyTo(Aq);
    datafloat normAq = 0;

    for(int n=0;n<E*BSIZE;++n)
      normAq += Aq[n]*Aq[n];
    normAq = sqrt(normAq);

    printf("OCCA: normAq = %17.15lf\n", normAq);

    for(int n=0;n<E*BSIZE;++n){
      Aq[n] = 0.0f;
    }

    o_Aq.copyFrom(Aq);

  }//fore


  //printf("\n\nBWfromCopy%d = [", E);
  for (int k=kMin; k<=kMax; k++){
    printf("==== this is kernel %d \n", k);
    int p_Nq = k+1;
    int p_gjNq = k+2;
    int Niter = 10;
    double gflops;
    long long int Nbytes;

    Nbytes =  7*p_Np*p_Np*sizeof(datafloat) + E*p_Np*p_Np*sizeof(datafloat);

    Nbytes /= 2;

    occa::memory o_foo = device.malloc(Nbytes);
    occa::memory o_bah = device.malloc(Nbytes);

    occa::streamTag startCopy = device.tagStream();
    for(int it=0;it<Niter;++it){
      o_bah.copyTo(o_foo);
    }
    occa::streamTag endCopy = device.tagStream();
    double copyElapsed = device.timeBetween(startCopy, endCopy);
    double copyBandwidth = ((Nbytes*Niter*2.)/(1.e9*copyElapsed));
    printf("copytime = %16.17f \n", copyElapsed);
    printf("pNp = %d p_Nfp = %d Nbytes = %d \n", p_Np, p_Nfp, Nbytes);    
    printf("copy BW = %f gflops = %f bytes = %d \n", copyBandwidth, gflops, Nbytes);

    roofline[k] = (copyBandwidth*gflops*E)/(2.*Nbytes);

    o_foo.free();
    o_bah.free();

  }

  //printf("];\n\n");

  printf("\n\nROOFLINE(%d,%d,%d,:)  = [", p_N,p_Ne,p_Nb);
  for (int k=kMin; k<=kMax; k++){

    printf(" %16.17f ", roofline[k]);
  }

  printf("]\n\n");


  printf("\n\nResults(%d,%d,%d,:)  = [", p_N, p_Ne, p_Nb);
  for (int k=kMin; k<=kMax; k++){

    printf(" %16.17f ", results3D[k]);
  }

  printf("];\n\n");

  printf("\n\nDOFSvGFLOPSvGDOFS(%d,:)  = [", p_N);
  for (int k=kMin; k<=kMax; k++){
    printf(" %d %16.17f %16.17f ", p_Np*E, results3D[k], (p_Np*E/1.e9)/timeData[k]);
  }

  printf("];\n\n");


  printf("\n\ntimeData(:, %d) = [",p_N);
  for (int k=kMin; k<=kMax; k++){

    printf(" %16.17f ", timeData[k]);
  }

  printf("];\n\n");

  printf("\n\n giganodes(%d,%d,%d,:)  = [", p_N,p_Ne,p_Nb);
  for (int k=kMin; k<=kMax; k++){

    printf(" %16.17f ", (p_Np*E/1.e9)/timeData[k]);
  }

  printf("];\n\n");


  exit(0);
  return 0;

}
