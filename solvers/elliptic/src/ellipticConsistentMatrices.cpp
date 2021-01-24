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

#include "elliptic.hpp"
#include "ellipticPrecon.hpp"

void ellipticBuildOperatorConsistentDiagonal(elliptic_t &elliptic, dfloat *diagA){

  mesh_t &mesh = elliptic.mesh;

  int useCubature = ((mesh.elementType==HEXAHEDRA||mesh.elementType==QUADRILATERALS) &&
			 (elliptic.settings.compareSetting("ELLIPTIC INTEGRATION", "CUBATURE") ? 1:0));
  
  int Np = mesh.Np;
  int Nelements = mesh.Nelements;

  //build kernels
  occa::properties kernelInfo = elliptic.platform.props;
  occa::kernel setLocalNodeKernel = 
    elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticSetLocalNode.okl", "ellipticSetLocalNodeKernel", kernelInfo);
  occa::kernel stridedCopyKernel = 
    elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticSetLocalNode.okl", "ellipticStridedCopyKernel", kernelInfo);
  
  dfloat *qL = (dfloat*) calloc(Np*Nelements, sizeof(dfloat));
  dfloat *AqL = (dfloat*) calloc(Np*Nelements, sizeof(dfloat));
  dfloat *diagAL = (dfloat*) calloc(Np*Nelements, sizeof(dfloat));

  occa::memory o_qL  = elliptic.platform.device.malloc(Np*Nelements*sizeof(dfloat));
  occa::memory o_AqL = elliptic.platform.device.malloc(Np*Nelements*sizeof(dfloat));
  occa::memory o_diagAL = elliptic.platform.device.malloc(Np*Nelements*sizeof(dfloat));  
  int useGlobalToLocal = 0;
  
  double tic = MPI_Wtime();
  
  for(dlong n=0;n<Np;++n){

    setLocalNodeKernel(Nelements, Np, n, (dfloat)1.0, o_qL);

    elliptic.platform.device.finish();
    //    printf("Stage %d\n", 1);
    
    if(mesh.NglobalGatherElements) {

      if(useCubature==0) { // GLL or non-hex
	printf("WARNING: using GLL kernel\n");
	elliptic.partialAxKernel(mesh.NglobalGatherElements,
				 mesh.o_globalGatherElementList,
				 elliptic.ogsMasked->o_GlobalToLocal,
				 useGlobalToLocal,
				 mesh.o_ggeo, mesh.o_D, mesh.o_S, mesh.o_MM, elliptic.lambda, o_qL, o_AqL);
      }else{
	elliptic.partialCubatureAxKernel(mesh.NglobalGatherElements,
					 mesh.o_globalGatherElementList,
					 elliptic.ogsMasked->o_GlobalToLocal,
					 useGlobalToLocal,
					 mesh.o_cubggeo,
					 mesh.o_cubD, // check layout
					 mesh.o_cubInterp, // check layout
					 elliptic.lambda, o_qL, o_AqL);
      }	
    }

    elliptic.platform.device.finish();
    //    printf("Stage %d\n", 2);
    
    if(mesh.NlocalGatherElements){
      if(useCubature==0) { // GLL or non-hex
	printf("WARNING: using GLL kernel\n");
	elliptic.partialAxKernel(mesh.NlocalGatherElements, 
				 mesh.o_localGatherElementList,
				 elliptic.ogsMasked->o_GlobalToLocal,
				 useGlobalToLocal,
				 mesh.o_ggeo, mesh.o_D, mesh.o_S, mesh.o_MM, elliptic.lambda, o_qL, o_AqL);
      }else{
	elliptic.partialCubatureAxKernel(mesh.NlocalGatherElements,					 
					 mesh.o_localGatherElementList,
					 elliptic.ogsMasked->o_GlobalToLocal,
					 useGlobalToLocal,
					 mesh.o_cubggeo,					 
					 mesh.o_cubD, //right layout
					 mesh.o_cubInterp, // dropped T ?
					 elliptic.lambda,
					 o_qL,
					 o_AqL);
      }
    }

    elliptic.platform.device.finish();
    //    printf("Stage %d\n", 3);
    
    // strided by Np, offset n
    stridedCopyKernel(Nelements, Np, n, o_AqL, o_diagAL);

    elliptic.platform.device.finish();
    //    printf("Stage %d\n", 4);

    //    elliptic.ogsMasked->Gather(diagA, o_diagAL, ogs_dfloat, ogs_add, ogs_trans);
    
  }
  
  o_diagAL.copyTo(diagAL);

  for(dlong e=0;e<Nelements;++e){
    for(dlong n=0;n<Np;++n){
      dlong id = e*Np+n;
      if(diagAL[id] == 0 || elliptic.mapB[id])
	diagAL[id] = 1;
      
      if (elliptic.allNeumann) {
	if (elliptic.mapB[id]!=1) { //dont fill rows for masked nodes
	  diagAL[id] += elliptic.allNeumannPenalty*elliptic.allNeumannScale*elliptic.allNeumannScale;
	}
      }
    }
  }
  
  elliptic.ogsMasked->Gather(diagA, diagAL, ogs_dfloat, ogs_add, ogs_trans);
  double toc = MPI_Wtime();

  printf("Diagonal build (sloppy) took %g seconds\n", toc-tic);
  
  free(qL);
  free(AqL);
  free(diagAL);
}

#define nonZero_t parAlmond::parCOO::nonZero_t

dlong ellipticBuildOperatorConsistentMatrix(elliptic_t &elliptic, nonZero_t *A){

  mesh_t &mesh = elliptic.mesh;

  int useCubature = ((mesh.elementType==HEXAHEDRA||mesh.elementType==QUADRILATERALS) &&
			 elliptic.settings.compareSetting("ELLIPTIC INTEGRATION", "CUBATURE")) ? 1:0;

  //  useCubature = 0;
  
  int Np = mesh.Np;
  int Nelements = mesh.Nelements;

  //build kernels
  occa::properties kernelInfo = elliptic.platform.props;

  //  printf("CONSISTENT MATRIX INTEGRATION TYPE: %d\n",  useCubature);
  //  std::cout << kernelInfo << std::endl;
  
  occa::kernel setLocalNodeKernel = 
    elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticSetLocalNode.okl", "ellipticSetLocalNodeKernel", kernelInfo);
  occa::kernel stridedCopyKernel = 
    elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticSetLocalNode.okl", "ellipticStridedCopyKernel", kernelInfo);
  
  dfloat *q = (dfloat*) calloc(Np*Nelements, sizeof(dfloat));

  occa::memory o_q  = elliptic.platform.device.malloc(Np*Nelements*sizeof(dfloat));
  occa::memory o_Aq = elliptic.platform.device.malloc(Np*Np*Nelements*sizeof(dfloat));

  double tic = MPI_Wtime();

  int useGlobalToLocal = 0;
  
#if 0
  dfloat *tmpCubInterp = (dfloat*) calloc(mesh.Nq*mesh.cubNq, sizeof(dfloat));
  dfloat *tmpCubD      = (dfloat*) calloc(mesh.cubNq*mesh.cubNq, sizeof(dfloat));
  mesh.o_cubInterp.copyTo(tmpCubInterp);
  mesh.o_cubD.copyTo(tmpCubD);
  printf("cubInterp=\n");
  for(int n=0;n<mesh.cubNq;++n){
    for(int m=0;m<mesh.Nq;++m){
      printf("%g, ", tmpCubInterp[n+ m*mesh.cubNq]);
    }
    printf("\n");
  }
  printf("cubD=\n");
  for(int n=0;n<mesh.cubNq;++n){
    for(int m=0;m<mesh.cubNq;++m){
      printf("%g, ", tmpCubD[n*mesh.cubNq+m]);
    }
    printf("\n");
  }
#endif
  
  for(dlong n=0;n<Np;++n){

    size_t offset = n*(size_t)(mesh.Nelements*mesh.Np*sizeof(dfloat));
    
    occa::memory o_Aqn = o_Aq+offset;
    
    setLocalNodeKernel(Nelements, Np, n, (dfloat)1.0, o_q);
    
    if(mesh.NglobalGatherElements) {

      if(useCubature==0) { // GLL or non-hex
	printf("WARNING: using partial Ax\n");
	elliptic.partialAxKernel(mesh.NglobalGatherElements,
				 mesh.o_globalGatherElementList,
				 elliptic.ogsMasked->o_GlobalToLocal,
				 useGlobalToLocal,
				 mesh.o_ggeo,
				 mesh.o_D,
				 mesh.o_S,
				 mesh.o_MM,
				 elliptic.lambda,
				 o_q,
				 o_Aqn);
      }else{
	printf("WARNING: using cubature partial Ax\n");
	elliptic.partialCubatureAxKernel(mesh.NglobalGatherElements,
					 mesh.o_globalGatherElementList,
					 elliptic.ogsMasked->o_GlobalToLocal,
					 useGlobalToLocal,
					 mesh.o_cubggeo,
					 mesh.o_cubD, // check layout
					 mesh.o_cubInterp, // check layout
					 elliptic.lambda,
					 o_q,
					 o_Aqn);
      }	
    }
    
    if(mesh.NlocalGatherElements){
      if(useCubature==0) { // GLL or non-hex
	printf("WARNING: using partial Ax\n");
	elliptic.partialAxKernel(mesh.NlocalGatherElements, mesh.o_localGatherElementList,
				 mesh.o_ggeo, mesh.o_D, mesh.o_S, mesh.o_MM, elliptic.lambda, o_q, o_Aqn);
      }else{
	printf("WARNING: using cubature partial Ax\n");
	elliptic.partialCubatureAxKernel(mesh.NlocalGatherElements,
					 mesh.o_localGatherElementList,
					 mesh.o_cubggeo,
					 mesh.o_cubD, //right layout
					 mesh.o_cubInterp, // dropped T ?
					 elliptic.lambda,
					 o_q,
					 o_Aqn);
      }
    }
  }

  double toc = MPI_Wtime();
  
  printf("BUILDING CONSISTENT MATRIX took %g seconds on device\n", toc-tic);
  
  dlong allEntries = mesh.Np*mesh.Np*mesh.Nelements;
  dfloat *tmpA = (dfloat*) calloc(allEntries, sizeof(dfloat));
  o_Aq.copyTo(tmpA);
  
  dlong cnt = 0;
  for(dlong e=0;e<mesh.Nelements;++e){
    for(int n=0;n<mesh.Np;++n){
      dlong row = elliptic.maskedGlobalNumbering[e*mesh.Np+n];
      for(int m=0;m<mesh.Np;++m){
	dlong col = elliptic.maskedGlobalNumbering[e*mesh.Np+m];
	cnt = e*mesh.Np*mesh.Np + n*mesh.Np + m;
	if(row>=0 && col>=0){
	  A[cnt].row = row;
	  A[cnt].col = col;
	  A[cnt].val = tmpA[m*mesh.Np*mesh.Nelements + e*mesh.Np + n];
	}
      }
    }
  }

  cnt = mesh.Nelements*mesh.Np*mesh.Np;
  
  free(tmpA);
  free(q);

  return cnt;
}
