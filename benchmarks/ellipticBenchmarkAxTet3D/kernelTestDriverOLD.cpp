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

  int NKernels = 9;

  // default to 512 elements if no arg is given
  int E = (argc>=2) ? atoi(argv[1]):512;
  int p_N = (argc>=3) ? atoi(argv[2]):5;
  int option = (argc>=4) ? atoi(argv[3]):1;  
  int p_Ne = (argc>=5) ? atoi(argv[4]):1;
  int p_Nb = (argc>=6) ? atoi(argv[5]):1;

  int kMin = (argc>=7) ? atoi(argv[6]):0;
  int kMax = (argc>=8) ? atoi(argv[7]):8;


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
  printf("Elements per block : %d\n", p_Ne);
  printf("Outputs per thread : %d\n", p_Nb);
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
  if (option == 0){
    gflops =  p_Np*24+p_Np*p_Np*8+p_NfacesNfp*37 + p_Np*p_NfacesNfp * 8;}
  else {
    gflops = p_Np*20*(1+p_Np);
  }

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


  if (option == 0){
    datafloat *vgeo, *sgeo, *gradq, *Aq, *LIFTT, *MM, *DsT, *DrT, *DtT ;
    int * elementList,* vmapM, *vmapP, *EToB; 
    occa::memory o_vgeo, o_sgeo, o_gradq, o_Aq, o_LIFTT, o_MM, o_DsT, o_DrT, o_DtT, o_vmapM, o_vmapP, o_EToB, o_elementList;

    int p_Nmax;
    if  (p_Np > p_Nfp*p_Nfaces){
      p_Nmax = p_Np;
    }
    else {p_Nmax = p_Nfp*p_Nfaces;}
    kernelInfo.addDefine("p_Nmax", p_Nmax);

    kernelInfo.addDefine("p_NXID", 0);
    kernelInfo.addDefine("p_NYID", 1);
    kernelInfo.addDefine("p_NZID", 2);
    kernelInfo.addDefine("p_SJID", 3);
    kernelInfo.addDefine("p_IJID", 4);
    kernelInfo.addDefine("p_IHID", 5);
    kernelInfo.addDefine("p_WSJID", 6);

    kernelInfo.addDefine("p_RXID", 0);
    kernelInfo.addDefine("p_RYID", 1);
    kernelInfo.addDefine("p_RZID", 2);
    kernelInfo.addDefine("p_SXID", 3);
    kernelInfo.addDefine("p_SYID", 4);
    kernelInfo.addDefine("p_SZID", 5);
    kernelInfo.addDefine("p_TXID", 6);
    kernelInfo.addDefine("p_TYID", 7);
    kernelInfo.addDefine("p_TZID", 8);
    kernelInfo.addDefine("p_JID" , 9);
    kernelInfo.addDefine("p_JWID", 10);

    char buf[200];
    for (int i =kMin; i<kMax+1; i++){
      printf("compiling kernel %d ...\n", i);
      sprintf(buf, "ellipticPartialAxIpdgTet3D_Ref%d", i); 
      Tet3Dkernel[i] = device.buildKernelFromSource("ellipticAxIpdgTet3D.okl", buf, kernelInfo);
    }

    // initialize with random numbers on Host and Device
    srand48(12345);

    randCalloc(device, E*BSIZE*11, &vgeo, o_vgeo);
    randCalloc(device, E*BSIZE*p_Nfaces*p_Nfp*7, &sgeo, o_sgeo);
    randCalloc(device, BSIZE*BSIZE, &DrT, o_DrT);
    randCalloc(device, BSIZE*BSIZE, &DsT, o_DsT);
    randCalloc(device, BSIZE*BSIZE, &DtT, o_DtT);
    randCalloc(device, BSIZE*p_Nfaces*p_Nfp, &LIFTT, o_LIFTT);
    randCalloc(device, BSIZE*BSIZE, &MM, o_MM);
    randCalloc(device, E*BSIZE*4, &gradq, o_gradq);
    randCalloc(device, E*BSIZE, &Aq, o_Aq);


    vmapM = (int *) calloc(E*p_Nfp*p_Nfaces, sizeof(int));
    vmapP = (int *) calloc(E*p_Nfp*p_Nfaces, sizeof(int));
    EToB = (int *) calloc(E*p_Nfaces, sizeof(int));
    elementList = (int *) calloc(E, sizeof(int));


    for(int j=0;j<E*p_Nfp*p_Nfaces;++j){
      vmapM[j] = rand()%BSIZE;
      vmapP[j] = rand()%BSIZE;
    }
    for(int j=0;j<E*p_Nfaces;++j){
      EToB[j] = rand()%2+1;
    }
    for(int j=0;j<E;++j){
      elementList[j] = j;
    }



    o_vmapM = device.malloc(E*p_Nfaces*p_Nfp*sizeof(int), vmapM);
    o_vmapP = device.malloc(E*p_Nfaces*p_Nfp*sizeof(int), vmapP);
    o_EToB = device.malloc(E*p_Nfaces*sizeof(int), EToB);
    o_elementList = device.malloc(E*sizeof(int), elementList);


    // initialize timer
    occa::initTimer(device);

    // queue Ax kernels
    for (int i =kMin;i<=kMax; i++){
      datafloat lambda = 1.;
      datafloat tau  = 0.5;    
      occa::streamTag startTag = device.tagStream();

      // launch kernel
      for(it=0;it<Niter;++it){  
        Tet3Dkernel[i](E, o_elementList,
		       o_vmapM,
		       o_vmapP,
		       lambda,
		       tau,
		       o_vgeo,
		       o_sgeo,
		       o_EToB,
		       o_DrT,
		       o_DsT,
		       o_DtT,
		       o_LIFTT,
		       o_MM,
		       o_gradq,
		       o_Aq);
      }

      occa::streamTag stopTag = device.tagStream();
      double elapsed = device.timeBetween(startTag, stopTag);

      printf("\n\nKERNEL %d  ================================================== \n\n", i);
      printf("OCCA elapsed time = %g\n", elapsed);
      printf("gflops = %f time = %f \n", gflops, elapsed);

      printf("GFL %17.17f \n",E*gflops/(elapsed*1.e9) );      

      results3D[i] = elapsed/Niter;

      printf("OCCA: estimated gflops = %17.15f\n", E*gflops/(elapsed*1.e9));

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

    }
  }
  else{
    datafloat  *ggeo, *SrrT,  *SrsT,  *SrtT,  *SssT,  *SsrT, *SstT,  *StsT, *StrT,  *SttT,  *MM, *q, *Aq;
    int  *elementList;
    occa::memory o_elementList, o_ggeo, o_SrrT, o_SrsT, o_SrtT, o_SsrT,o_SssT,o_SstT, o_StrT, o_StsT, o_SttT, o_MM, o_q, o_Aq; 

    srand(12345);

    kernelInfo.addDefine("p_G00ID", 0);
    kernelInfo.addDefine("p_G01ID", 1);
    kernelInfo.addDefine("p_G02ID", 2);
    kernelInfo.addDefine("p_G11ID", 3);
    kernelInfo.addDefine("p_G12ID", 4);
    kernelInfo.addDefine("p_G22ID", 5);
    kernelInfo.addDefine("p_GWJID", 6);
    randCalloc(device, E*7, &ggeo, o_ggeo);
    randCalloc(device, BSIZE*BSIZE, &SrrT, o_SrrT);
    randCalloc(device, BSIZE*BSIZE, &SrsT, o_SrsT);
    randCalloc(device, BSIZE*BSIZE, &SrtT, o_SrtT);
    randCalloc(device, BSIZE*BSIZE, &SsrT, o_SsrT);
    randCalloc(device, BSIZE*BSIZE, &SssT, o_SssT);
    randCalloc(device, BSIZE*BSIZE, &SstT, o_SstT);
    randCalloc(device, BSIZE*BSIZE, &StrT, o_StrT);
    randCalloc(device, BSIZE*BSIZE, &StsT, o_StsT);
    randCalloc(device, BSIZE*BSIZE, &SttT, o_SttT);
    randCalloc(device, BSIZE*E, &q, o_q);
    randCalloc(device, BSIZE*E, &Aq, o_Aq);
    randCalloc(device, BSIZE*BSIZE, &MM, o_MM);

    elementList = (int *) calloc(E, sizeof(int));

    for(int j=0;j<E;++j){
      elementList[j] = j;
    }
    o_elementList = device.malloc(E*sizeof(int), elementList);
    char buf[200];
    for (int i =kMin; i<kMax+1; i++){
      printf("compiling 3D kernel %d ...\n", i);
      sprintf(buf, "ellipticPartialAxTet3D_Ref%d", i); 
      Tet3Dkernel[i] = device.buildKernelFromSource("ellipticAxTet3D.okl", buf, kernelInfo);
    }
    int elementOffset = 0;
    occa::initTimer(device);

    // queue Ax kernels
    for (int i =kMin;i<=kMax; i++){
      datafloat lambda = 1.;
      occa::streamTag startTag = device.tagStream();
      if (i<9){
	// launch kernel
	for(it=0;it<Niter;++it){
	  Tet3Dkernel[i](E,
			 o_elementList,
			 o_ggeo,
			 o_SrrT,
			 o_SrsT,
			 o_SrtT,
			 o_SsrT,
			 o_SssT,
			 o_SstT,
			 o_StrT,
			 o_StsT,
			 o_SttT,
			 o_MM,
			 lambda,
			 o_q,       
			 o_Aq);
	}
      }
      else{
	for(it=0;it<Niter;++it){
	  Tet3Dkernel[i](E,
			 elementOffset,
			 o_ggeo,
			 o_SrrT,
			 o_SrsT,
			 o_SrtT,
			 o_SsrT,
			 o_SssT,
			 o_SstT,
			 o_StrT,
			 o_StsT,
			 o_SttT,
			 o_MM,
			 lambda,
			 o_q,       
			 o_Aq);
	}


      }

      //adjust GFLOPS for Ref2 and Ref3
      occa::streamTag stopTag = device.tagStream();
      double elapsed = device.timeBetween(startTag, stopTag);
<<<<<<< HEAD

 if (i>5){
=======
      if (i>5){
>>>>>>> 9969b2dadca2a1aeeb06960203d30f8d145d3aab
        // old: gflops = p_Np*20*(1+p_Np)
        gflops = p_Np*(p_Np*14 +14);
	gflops *=Niter;      


      } 
      printf("\n\nKERNEL %d  ================================================== \n\n", i);
      printf("OCCA elapsed time = %g\n", elapsed);
      timeData[i] = elapsed/Niter;
      printf("number of flops = %f time = %f \n", gflops, elapsed);
<<<<<<< HEAD
printf("DOFS %lu  %17.17f \n",p_Np*E,  elapsed);    

      results3D[i] =E*gflops/(elapsed*1000*1000*1000); 
//elapsed/Niter;
=======
      results3D[i] =E*gflops/(elapsed*1.e9); 
      //elapsed/Niter;
>>>>>>> 9969b2dadca2a1aeeb06960203d30f8d145d3aab
      //
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

    }//for




  }//else


  //printf("\n\nBWfromCopy%d = [", E);
  for (int k=kMin; k<=kMax; k++){
    printf("==== this is kernel %d \n", k);
    int p_Nq = k+1;
    int p_gjNq = k+2;
    int Niter = 10;
    double gflops;
    long long int Nbytes;
    if (option == 0){
      gflops =
        p_Np*24+p_Np*p_Np*8+p_NfacesNfp*37 + p_Np*p_NfacesNfp * 8;

      Nbytes = E*(p_Nfp*p_Nfaces*2*sizeof(int) +
		  10*sizeof(datafloat)+ 6*p_Nfaces*sizeof(datafloat)+
		  p_Nfaces*sizeof(int) +1*sizeof(int)+5*p_Np*sizeof(datafloat))+4*p_Np*p_Np*sizeof(datafloat)+p_Np*p_Nfp*p_Nfaces*sizeof(datafloat); 

      Nbytes /= 2;

      //printf("Building %d bytes\n", Nbytes*E);

    }
    else {
      if (k<=5)
	Nbytes = 10*p_Np*p_Np*sizeof(datafloat) + E*7*sizeof(datafloat)+E*sizeof(int)+E*p_Np*2*sizeof(datafloat);
      else {
	Nbytes =  7*p_Np*p_Np*sizeof(datafloat) + E*7*sizeof(datafloat)+E*sizeof(int)+E*p_Np*2*sizeof(datafloat);
      }

      Nbytes /= 2;
      gflops = p_Np*20*(1+p_Np); 
      if (k>5){
        // old: gflops = p_Np*20*(1+p_Np)
        gflops = p_Np*(p_Np*14 +14);
      }    

    }
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
    //  printf(" %16.17f ", copyBandwidth);

    //    Nbytes = E*(p_Nfp*p_Nfaces*2*sizeof(int) +
    //      10*sizeof(datafloat)+ 6*p_Nfaces*sizeof(datafloat)+
    //    p_Nfaces*sizeof(int) +1*sizeof(int)+5*p_Np*sizeof(datafloat))+4*p_Np*p_Np*sizeof(datafloat)+p_Np*p_Nfp*p_Nfaces*sizeof(datafloat); 
    printf("pNp = %d p_Nfp = %d Nbytes = %d \n", p_Np, p_Nfp, Nbytes);    
    printf("copy BW = %f gflops = %f bytes = %d \n", copyBandwidth, gflops, Nbytes);
    //    'roofline[k-1] = copyElapsed/Niter; 
    roofline[k] = (copyBandwidth*gflops*E)/(2.*Nbytes);

    //((E*gflops*(double)Niter))/(1e9*copyElapsed);
    o_foo.free();
    o_bah.free();

  }

  //printf("];\n\n");

  printf("\n\nROOFLINE(%d,%d,%d,:)  = [", p_N,p_Ne,p_Nb);
  for (int k=kMin; k<=kMax; k++){

    printf(" %16.17f ", roofline[k]);
  }

  printf("]\n\n");


  printf("\n\nResults(%d,%d,%d,:)  = [", p_N,p_Ne,p_Nb);
  for (int k=kMin; k<=kMax; k++){

    printf(" %16.17f ", results3D[k]);
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
