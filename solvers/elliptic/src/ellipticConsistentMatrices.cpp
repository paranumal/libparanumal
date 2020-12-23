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

  int integrationType = ((mesh.elementType==HEXAHEDRA||mesh.elementType==QUADRILATERALS) &&
			 elliptic.settings.compareSetting("ELLIPTIC INTEGRATION", "CUBATURE")) ? 1:0;

  int Np = mesh.Np;
  int Nelements = mesh.Nelements;

  //build kernels
  occa::properties kernelInfo = elliptic.platform.props;
  occa::kernel setLocalNodeKernel = 
    elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticSetLocalNode.okl", "ellipticSetLocalNodeKernel", kernelInfo);
  occa::kernel stridedCopyKernel = 
    elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticSetLocalNode.okl", "ellipticStridedCopyKernel", kernelInfo);
  
  dfloat *q = (dfloat*) calloc(Np*Nelements, sizeof(dfloat));
  dfloat *Aq = (dfloat*) calloc(Np*Nelements, sizeof(dfloat));

  occa::memory o_q  = elliptic.platform.device.malloc(Np*Nelements*sizeof(dfloat));
  occa::memory o_Aq = elliptic.platform.device.malloc(Np*Nelements*sizeof(dfloat));
  occa::memory o_diagA = elliptic.platform.device.malloc(Np*Nelements*sizeof(dfloat));  

  double tic = MPI_Wtime();
  
  for(dlong n=0;n<Np;++n){

    setLocalNodeKernel(Nelements, Np, n, (dfloat)1.0, o_q);

    if(mesh.NglobalGatherElements) {

      if(integrationType==0) { // GLL or non-hex
	elliptic.partialAxKernel(mesh.NglobalGatherElements, mesh.o_globalGatherElementList,
				 mesh.o_ggeo, mesh.o_D, mesh.o_S, mesh.o_MM, elliptic.lambda, o_q, o_Aq);
      }else{
	elliptic.partialCubatureAxKernel(mesh.NglobalGatherElements,
					 mesh.o_globalGatherElementList,
					 mesh.o_cubggeo,
					 mesh.o_cubD, // check layout
					 mesh.o_cubInterp, // check layout
					 elliptic.lambda, o_q, o_Aq);
      }	
    }
    
    if(mesh.NlocalGatherElements){
      if(integrationType==0) { // GLL or non-hex
	elliptic.partialAxKernel(mesh.NlocalGatherElements, mesh.o_localGatherElementList,
				 mesh.o_ggeo, mesh.o_D, mesh.o_S, mesh.o_MM, elliptic.lambda, o_q, o_Aq);
      }else{
	elliptic.partialCubatureAxKernel(mesh.NlocalGatherElements,
					 mesh.o_localGatherElementList,
					 mesh.o_cubggeo,
					 mesh.o_cubD, //right layout
					 mesh.o_cubInterp, // dropped T ?
					 elliptic.lambda,
					 o_q,
					 o_Aq);
      }
    }

    // strided by Np, offset n
    stridedCopyKernel(Nelements, Np, n, o_Aq, o_diagA);
  }
  
  o_diagA.copyTo(diagA);

  for(dlong e=0;e<Nelements;++e){
    for(dlong n=0;n<Np;++n){
      dlong id = e*Np+n;
      if(diagA[id] == 0)
	diagA[id] = 1;
    }
  }
  
  elliptic.ogsMasked->GatherScatter(diagA, ogs_dfloat, ogs_add, ogs_sym);
  double toc = MPI_Wtime();

  printf("Diagonal build (sloppy) took %g seconds\n", toc-tic);
  
  free(q);
  free(Aq);
}

#define nonZero_t parAlmond::parCOO::nonZero_t

dlong ellipticBuildOperatorConsistentMatrix(elliptic_t &elliptic, nonZero_t *A){

  mesh_t &mesh = elliptic.mesh;

  int integrationType = ((mesh.elementType==HEXAHEDRA||mesh.elementType==QUADRILATERALS) &&
			 elliptic.settings.compareSetting("ELLIPTIC INTEGRATION", "CUBATURE")) ? 1:0;

  int Np = mesh.Np;
  int Nelements = mesh.Nelements;

  //build kernels
  occa::properties kernelInfo = elliptic.platform.props;
  occa::kernel setLocalNodeKernel = 
    elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticSetLocalNode.okl", "ellipticSetLocalNodeKernel", kernelInfo);
  occa::kernel stridedCopyKernel = 
    elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticSetLocalNode.okl", "ellipticStridedCopyKernel", kernelInfo);
  
  dfloat *q = (dfloat*) calloc(Np*Nelements, sizeof(dfloat));
  dfloat *Aq = (dfloat*) calloc(Np*Nelements, sizeof(dfloat));

  occa::memory o_q  = elliptic.platform.device.malloc(Np*Nelements*sizeof(dfloat));
  occa::memory o_Aq = elliptic.platform.device.malloc(Np*Np*Nelements*sizeof(dfloat));

  double tic = MPI_Wtime();
  
  for(dlong n=0;n<Np;++n){

    size_t offset = n*mesh.Np*mesh.Nelements*sizeof(dfloat);
    
    setLocalNodeKernel(Nelements, Np, n, (dfloat)1.0, o_q);

    
    if(mesh.NglobalGatherElements) {

      if(integrationType==0) { // GLL or non-hex
	elliptic.partialAxKernel(mesh.NglobalGatherElements, mesh.o_globalGatherElementList,
				 mesh.o_ggeo, mesh.o_D, mesh.o_S, mesh.o_MM, elliptic.lambda, o_q, o_Aq+offset);
      }else{
	elliptic.partialCubatureAxKernel(mesh.NglobalGatherElements,
					 mesh.o_globalGatherElementList,
					 mesh.o_cubggeo,
					 mesh.o_cubD, // check layout
					 mesh.o_cubInterp, // check layout
					 elliptic.lambda, o_q, o_Aq+offset);
      }	
    }
    
    if(mesh.NlocalGatherElements){
      if(integrationType==0) { // GLL or non-hex
	elliptic.partialAxKernel(mesh.NlocalGatherElements, mesh.o_localGatherElementList,
				 mesh.o_ggeo, mesh.o_D, mesh.o_S, mesh.o_MM, elliptic.lambda, o_q, o_Aq+offset);
      }else{
	elliptic.partialCubatureAxKernel(mesh.NlocalGatherElements,
					 mesh.o_localGatherElementList,
					 mesh.o_cubggeo,
					 mesh.o_cubD, //right layout
					 mesh.o_cubInterp, // dropped T ?
					 elliptic.lambda,
					 o_q,
					 o_Aq+offset);
      }
    }
  }

  double toc = MPI_Wtime();
  
  printf("BUILDING CONSISTENT MATRIX took %g seconds on device\n", toc-tic);
  
  dlong allEntries = mesh.Np*mesh.Np*mesh.Nelements;
  dfloat *tmpA = (dfloat*) calloc(allEntries, sizeof(dfloat));
  o_Aq.copyTo(tmpA);
  
  dlong cnt = 0;
  for(int m=0;m<mesh.Np;++m){
    for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.Np;++n){
	dlong id = e*mesh.Np + n;
	dlong row = elliptic.maskedGlobalNumbering[id];
	dlong col = elliptic.maskedGlobalNumbering[e*mesh.Np+m];
	if(row>=0 && col>=0){
	  A[cnt].row = row;
	  A[cnt].col = col;
	  A[cnt].val = tmpA[m*mesh.Np*mesh.Nelements + id];
	  ++cnt;
	}
      }
    }
  }

  free(tmpA);
  free(q);
  free(Aq);

  return cnt;
}
