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

#include "test.hpp"
#include "mesh.hpp"

#define dlong int
#define dlongString "int"

#define dfloat double
#define ogs_dfloat ogs_double
#define MPI_DFLOAT MPI_DOUBLE
#define dfloatFormat "%lf"
#define dfloatString "double"

#define hlong long long int
#define ogs_hlong ogs_long_long
#define MPI_HLONG MPI_LONG_LONG_INT
#define hlongFormat "%lld"
#define hlongString "long long int"


/* offsets for geometric factors */
#define RXID 0
#define RYID 1
#define SXID 2
#define SYID 3
#define  JID 4
#define JWID 5
#define IJWID 6
#define RZID 7
#define SZID 8
#define TXID 9
#define TYID 10
#define TZID 11

/* offsets for second order geometric factors */
#define G00ID 0
#define G01ID 1
#define G11ID 2
#define GWJID 3
#define G12ID 4
#define G02ID 5
#define G22ID 6


/* offsets for nx, ny, sJ, 1/J */
#define NXID 0
#define NYID 1
#define SJID 2
#define IJID 3
#define IHID 4
#define WSJID 5
#define WIJID 6
#define NZID 7
#define STXID 8
#define STYID 9
#define STZID 10
#define SBXID 11
#define SBYID 12
#define SBZID 13
#define SURXID 14
#define SURYID 15
#define SURZID 16




void randDeviceCalloc(occa::device &device, int sz, dfloat **pt, occa::memory &o_pt){

  *pt = (dfloat*) calloc(sz, sizeof(dfloat));

  dfloat sum = 0;
  for(int n=0;n<sz;++n){
    (*pt)[n] = drand48()-0.5;
    sum += pow((*pt)[n],2);
  }

  o_pt = device.malloc(sz*sizeof(dfloat), *pt);

}


void setDeviceProperties(occa::device &device, occa::properties &props){

  props["defines/" "dfloat"]="double";
  props["defines/" "dfloat2"]="double2";
  props["defines/" "dfloat4"]="double4";
  props["defines/" "dfloat8"]="double8";

  props["defines/" "dlong"]="int";

  props["defines/" "p_NXID"]= NXID;
  props["defines/" "p_NYID"]= NYID;
  props["defines/" "p_NZID"]= NZID;
  props["defines/" "p_SJID"]= SJID;
  props["defines/" "p_IJID"]= IJID;
  props["defines/" "p_IHID"]= IHID;
  props["defines/" "p_WSJID"]= WSJID;
  props["defines/" "p_WIJID"]= WIJID;
  props["defines/" "p_STXID"]= STXID;
  props["defines/" "p_STYID"]= STYID;
  props["defines/" "p_STZID"]= STZID;
  props["defines/" "p_SBXID"]= SBXID;
  props["defines/" "p_SBYID"]= SBYID;
  props["defines/" "p_SBZID"]= SBZID;

  props["defines/" "p_G00ID"]= G00ID;
  props["defines/" "p_G01ID"]= G01ID;
  props["defines/" "p_G02ID"]= G02ID;
  props["defines/" "p_G11ID"]= G11ID;
  props["defines/" "p_G12ID"]= G12ID;
  props["defines/" "p_G22ID"]= G22ID;
  props["defines/" "p_GWJID"]= GWJID;


  props["defines/" "p_RXID"]= RXID;
  props["defines/" "p_SXID"]= SXID;
  props["defines/" "p_TXID"]= TXID;

  props["defines/" "p_RYID"]= RYID;
  props["defines/" "p_SYID"]= SYID;
  props["defines/" "p_TYID"]= TYID;

  props["defines/" "p_RZID"]= RZID;
  props["defines/" "p_SZID"]= SZID;
  props["defines/" "p_TZID"]= TZID;

  props["defines/" "p_JID"]= JID;
  props["defines/" "p_JWID"]= JWID;
  props["defines/" "p_IJWID"]= IJWID;

  props["compiler_flags"] += "--ftz=true ";
  props["compiler_flags"] += "--prec-div=false ";
  props["compiler_flags"] += "--prec-sqrt=false ";
  props["compiler_flags"] += "--use_fast_math ";
  props["compiler_flags"] += "--fmad=true "; // compiler option for cuda
  props["compiler_flags"] += "-Xptxas -dlcm=ca";
}

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  int NKernels = 1;

  // default to 512 elements if no arg is given
  dlong E     = (argc>=2) ? atoi(argv[1]):512;
  int N       = (argc>=3) ? atoi(argv[2]):5;
  int option  = (argc>=4) ? atoi(argv[3]):1;  
  int Ne      = (argc>=5) ? atoi(argv[4]):1;
  int Nb      = (argc>=6) ? atoi(argv[5]):1;

  int kMin = (argc>=7) ? atoi(argv[6]):0;
  int kMax = (argc>=8) ? atoi(argv[7]):2;

  int Np = (N+1)*(N+1);
  int Nfp = (N+1);
  int Nfaces = 4;
  int NfacesNfp = Nfaces*Nfp;
  int Nq = N+1;
  int Nggeo = 7;
  
  printf("==============================================================\n");
  printf("===================== BASIC INFO =============================\n");
  printf("==============================================================\n");
  printf("Number of elements : %d\n", E);
  printf("Polynomial degree  : %d\n", N);
  printf("Nodes per element  : %d\n", Np);
  printf("Elements per block : %d\n", Ne);
  printf("Outputs per thread : %d\n", Nb);
  printf("Running kernels    : %d to %d \n", kMin, kMax);
  printf("==============================================================\n");
  printf("==============================================================\n");
  printf("==============================================================\n");
  printf("\n\n");
  int BSIZE  = Np;

  // number of geometric factors
  int Nvgeo =  Nggeo;
  int Niter = 10, it;
  double gflops;
  
  gflops = Np*(4*Nq + 8 + 4*Nq); // FIX LATER
  gflops *= Niter;

  // build some dummy storage & parameters
  double results[15];
  double roofline[15];
  double timeData[15];

  occa::device device;
  occa::properties kernelInfo;
  // occa::kernel ellipticAxKernel;
  occa::kernel ellipticAxKernel[15];
  occa::kernel correctRes;

  // specify device configuration
   char deviceConfig[BUFSIZ];
   int device_id   = 0; 
   int platform_id = 0; 
   
   sprintf(deviceConfig, "mode: 'CUDA', device_id: %d", device_id);
   // sprintf(deviceConfig, "mode: 'OPENCL', platform_id = %d, device_id: %d", platform_id, device_id);
   device.setup( (std::string)deviceConfig);
   
   setDeviceProperties(device, kernelInfo); 
  // device.setup("mode = Serial");

  kernelInfo["defines/" "p_N"]          = N;
  kernelInfo["defines/" "p_Nq"]         = Nq;
  kernelInfo["defines/" "p_Np"]         = Np;
  kernelInfo["defines/" "p_Nfp"]        = Nfp;
  kernelInfo["defines/" "p_Nfaces"]     = Nfaces;
  kernelInfo["defines/" "p_NfacesNfp"]  = Nfaces*Nfp;
  kernelInfo["defines/" "BSIZE"]        = (int)BSIZE;
  kernelInfo["defines/" "p_Nggeo"]      = Nggeo;
  
  kernelInfo["defines/" "p_Ne"]         = Ne;
  kernelInfo["defines/" "p_Nb"]         = Nb;

  dlong *elementList = (dlong *)calloc(E, sizeof(dlong));
  for(dlong e = 0; e<E; e++) elementList[e] = e; 
  occa::memory o_elementList = device.malloc(E*sizeof(dlong), elementList);

  dfloat    *ggeo, *DT, *S, *MM,  *coeff, *q, *Aq;
  occa::memory o_ggeo, o_DT, o_S, o_MM, o_coeff, o_q, o_Aq; 
  
  srand(12345);
  
  randDeviceCalloc(device, E*Nggeo*BSIZE*BSIZE, &ggeo, o_ggeo);
  randDeviceCalloc(device, BSIZE*BSIZE, &DT, o_DT);
  randDeviceCalloc(device, BSIZE*BSIZE, &S, o_S);
  randDeviceCalloc(device, BSIZE*BSIZE, &MM, o_MM);
  randDeviceCalloc(device, BSIZE*E, &q, o_q);
  randDeviceCalloc(device, BSIZE*E, &Aq, o_Aq);
  randDeviceCalloc(device, BSIZE*E*2, &coeff, o_coeff);

  // char buf[200];
  char fileName[BUFSIZ], kernelName[BUFSIZ];
  // for (int k =kMin; k<1; k++){
  for (int k =kMin; k<kMax+1; k++){
    sprintf(fileName, DTEST "/okl/ellipticAxHex3D.okl");
    sprintf(kernelName, "ellipticPartialAxHex3D_%d", k);
    ellipticAxKernel[k] = buildKernel(device, fileName, kernelName, kernelInfo, comm); 
  }
  
  occa::streamTag startTag, stopTag; 
  // queue Ax kernels
  for (int k =kMin;k<=kMax; k++){
    if(k==0){
    dfloat lambda = 1.;
    startTag = device.tagStream();

    for(it=0;it<Niter;++it){
      ellipticAxKernel[k](E,
                          o_elementList,
                          o_ggeo,
                          o_DT,
                          o_S,
                          o_MM,
                          lambda,
                          o_q,       
                          o_Aq);
    }
    stopTag = device.tagStream();
  }else{

    startTag = device.tagStream();

    for(it=0;it<Niter;++it){
      ellipticAxKernel[k](E,
                          o_elementList,
                          o_ggeo,
                          o_DT,
                          o_S,
                          o_MM,
                          o_coeff,
                          o_q,       
                          o_Aq);
    }
    stopTag = device.tagStream();
  }
    

    //adjust GFLOPS for Ref2 and Ref3
    double elapsed = device.timeBetween(startTag, stopTag);

    gflops = Np*(4*Nq + 8 + 4*Nq); // FIX LATER
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
    int pNq = k+1;
    int pgjNq = k+2;
    int Niter = 10;
    double gflops;
    long long int Nbytes;

    Nbytes = Nq*Nq*sizeof(dfloat) + E*Np*Nggeo*sizeof(dfloat)+E*sizeof(int)+E*Np*2*sizeof(dfloat);

    
    Nbytes /= 2;
    gflops = Np*(4*Nq + 8 + 4*Nq); // FIX LATER

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
    printf("pNp = %d p_Nfp = %d Nbytes = %d \n", Np, Nfp, Nbytes);    
    printf("copy BW = %f gflops = %f bytes = %d \n", copyBandwidth, gflops, Nbytes);
    roofline[k] = (copyBandwidth*gflops*E)/(2*Nbytes);

    o_foo.free();
    o_bah.free();
  }

  // printf("];\n\n");
  
  // printf("\n\nROOFLINE = [");
  // for (int k=kMin; k<=kMax; k++){
    
  //   printf(" %16.17f ", roofline[k]);
  // }
  
  // printf("]\n\n");


  // printf("\n\nResults(:,%d)  = [", N);
  // for (int k=kMin; k<=kMax; k++){

  //   printf(" %16.17f ", results[k]);
  // }

  // printf("];\n\n");

  printf("\n\ntimeData(:, %d, %d) = [",E, N);
  for (int k=kMin; k<=kMax; k++){

    printf(" %16.17f ", timeData[k]);
  }

  printf("];\n\n");

  // printf("\n\ngigaNodes(:, %d) = [",N);
  // for (int k=kMin; k<=kMax; k++){

  //   printf(" %16.17f ", (Np*E/10e9)/timeData[k]);
  // }

  // printf("];\n\n");

// close down MPI
  MPI_Finalize();
  return 0;

}
