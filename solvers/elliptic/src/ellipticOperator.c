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

#include "elliptic.h"
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "ogsInterface.h"

#include "omp.h"

#if 0
// 2 thread version
void ellipticOperator(elliptic_t *elliptic, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *precision){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  const cgOptions_t &cgOptions = elliptic->cgOptions;

  dfloat *sendBuffer = elliptic->sendBuffer;
  dfloat *recvBuffer = elliptic->recvBuffer;
  dfloat *gradSendBuffer = elliptic->gradSendBuffer;
  dfloat *gradRecvBuffer = elliptic->gradRecvBuffer;

  dfloat alpha = 0., alphaG = 0.;
  dlong Nblock = elliptic->Nblock;
  dfloat *tmp = elliptic->tmp;
  occa::memory &o_tmp = elliptic->o_tmp;

  int Nbytes = sizeof(dfloat);
	
  if(cgOptions.continuous){
    ogs_t *ellipticOgs = elliptic->ogs;

#pragma omp parallel num_threads(2)
    {
      int tid = omp_get_thread_num();
      
      if(tid==0){
	
	if(mesh->NglobalGatherElements) {
	  if(!cgOptions.serial || elliptic->elementType!=HEXAHEDRA)
	    elliptic->partialAxKernel(mesh->NglobalGatherElements, mesh->o_globalGatherElementList,
				      mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq);
	  else
	    ellipticSerialPartialAxHexKernel3D(mesh->Nq, mesh->NglobalGatherElements, mesh->o_globalGatherElementList,
					       mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq);
	}
	
	// Gather on device halo nodes
	occaGather(ellipticOgs->NhaloGather, ellipticOgs->o_haloGatherOffsets, ellipticOgs->o_haloGatherIds, dfloatString, "add", o_Aq, ogs::o_haloBuf);
	
	// copy device halo to host
	if(!cgOptions.serial){
	  ogs::o_haloBuf.copyTo(ogs::haloBuf, ellipticOgs->NhaloGather*Nbytes, 0, "async: true");
	  
	  // host gather scatter using libgs 
	  ogsHostGatherScatter(ogs::haloBuf, dfloatString, "add", ellipticOgs->haloGshSym);
	  
	  // copy totally gather halo data back from HOST to DEVICE
	  ogs::o_haloBuf.copyFrom(ogs::haloBuf, ellipticOgs->NhaloGather*Nbytes, 0, "async: true");
	}
	else{
	  // host gather scatter using libgs (directly on serial halo buffer)
	  ogsHostGatherScatter(ogs::o_haloBuf.ptr(), dfloatString, "add", ellipticOgs->haloGshSym);
	}

	// do scatter back to local nodes
	occaScatter(ellipticOgs->NhaloGather, ellipticOgs->o_haloGatherOffsets, ellipticOgs->o_haloGatherIds, dfloatString, "add", ogs::o_haloBuf, o_Aq);
      }

      if(tid==1){
      	if(mesh->NlocalGatherElements){

	  if(!cgOptions.serial || elliptic->elementType!=HEXAHEDRA){
	    elliptic->partialAxKernel2(mesh->NlocalGatherElements, mesh->o_localGatherElementList,
				       mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq);
	  }
	  else
	    ellipticSerialPartialAxHexKernel3D(mesh->Nq, mesh->NlocalGatherElements, mesh->o_localGatherElementList,
					       mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq);
	}
      }
    }    
     
    // have to do this after data from gslib gather scatter because owned nodes
    if(ellipticOgs->NlocalGather) {
      
      if(!cgOptions.serial)
	occaGatherScatter(ellipticOgs->NlocalGather, ellipticOgs->o_localGatherOffsets, ellipticOgs->o_localGatherIds, dfloatString, "add", o_Aq);
      else{
	const int Ngather = ellipticOgs->NlocalGather;
	const int * __restrict__ cpu_gatherStarts = (int*) ellipticOgs->o_localGatherOffsets.ptr();
	const int * __restrict__ cpu_gatherIds    = (int*) ellipticOgs->o_localGatherIds.ptr();
	dfloat * __restrict__ cpu_Aq = (dfloat*) o_Aq.ptr();

	cpu_Aq = (dfloat*)__builtin_assume_aligned(cpu_Aq, USE_OCCA_MEM_BYTE_ALIGN) ;
	cpu_gatherStarts = (int*)__builtin_assume_aligned(cpu_gatherStarts, USE_OCCA_MEM_BYTE_ALIGN) ;
	cpu_gatherIds = (int*)__builtin_assume_aligned(cpu_gatherIds, USE_OCCA_MEM_BYTE_ALIGN) ;
	
	for(dlong g=0;g<Ngather;++g){		
	  
	  const dlong start = cpu_gatherStarts[g];
	  const dlong end = cpu_gatherStarts[g+1];
	  if((start+1)!=end){
	    
	    double gAq = 0;
	    
	    for(dlong n=start;n<end;++n){
	      const dlong id = cpu_gatherIds[n];
	      gAq += cpu_Aq[id];
	    }
	  
	    for(dlong n=start;n<end;++n){
	      const dlong id = cpu_gatherIds[n];
	      cpu_Aq[id] = gAq;
	    }
	  }
	}	 
      }
    }
    
    if(elliptic->allNeumann) {
      
      elliptic->innerProductKernel(mesh->Nelements*mesh->Np, elliptic->o_invDegree, o_q, o_tmp);
      o_tmp.copyTo(tmp);
      
      for(dlong n=0;n<Nblock;++n)
	alpha += tmp[n];
      
      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
      alphaG *= elliptic->allNeumannPenalty*elliptic->allNeumannScale*elliptic->allNeumannScale;
      
      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);
    }
    
    //post-mask
    if (elliptic->Nmasked) 
      mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_Aq);
    
  } else {
    dlong offset = 0;
    dfloat alpha = 0., alphaG =0.;
    dlong Nblock = elliptic->Nblock;
    dfloat *tmp = elliptic->tmp;
    occa::memory &o_tmp = elliptic->o_tmp;

    ellipticStartHaloExchange(elliptic, o_q, mesh->Np, sendBuffer, recvBuffer);

    if(options.compareArgs("BASIS", "NODAL")) {
      elliptic->partialGradientKernel(mesh->Nelements,
				      offset,
				      mesh->o_vgeo,
				      mesh->o_Dmatrices,
				      o_q,
				      elliptic->o_grad);
    } else if(options.compareArgs("BASIS", "BERN")) {
      elliptic->partialGradientKernel(mesh->Nelements,
				      offset,
				      mesh->o_vgeo,
				      mesh->o_D1ids,
				      mesh->o_D2ids,
				      mesh->o_D3ids,
				      mesh->o_Dvals,
				      o_q,
				      elliptic->o_grad);
    }

    ellipticInterimHaloExchange(elliptic, o_q, mesh->Np, sendBuffer, recvBuffer);

    //Start the rank 1 augmentation if all BCs are Neumann
    //TODO this could probably be moved inside the Ax kernel for better performance
    if(elliptic->allNeumann)
      mesh->sumKernel(mesh->Nelements*mesh->Np, o_q, o_tmp);
    
    if(mesh->NinternalElements) {
      if(options.compareArgs("BASIS", "NODAL")) {
	elliptic->partialIpdgKernel(mesh->NinternalElements,
				    mesh->o_internalElementIds,
				    mesh->o_vmapM,
				    mesh->o_vmapP,
				    lambda,
				    elliptic->tau,
				    mesh->o_vgeo,
				    mesh->o_sgeo,
				    elliptic->o_EToB,
				    mesh->o_Dmatrices,
				    mesh->o_LIFTT,
				    mesh->o_MM,
				    elliptic->o_grad,
				    o_Aq);
      } else if(options.compareArgs("BASIS", "BERN")) {
	elliptic->partialIpdgKernel(mesh->NinternalElements,
				    mesh->o_internalElementIds,
				    mesh->o_vmapM,
				    mesh->o_vmapP,
				    lambda,
				    elliptic->tau,
				    mesh->o_vgeo,
				    mesh->o_sgeo,
				    elliptic->o_EToB,
				    mesh->o_D1ids,
				    mesh->o_D2ids,
				    mesh->o_D3ids,
				    mesh->o_Dvals,
				    mesh->o_L0vals,
				    mesh->o_ELids,
				    mesh->o_ELvals,
				    mesh->o_BBMM,
				    elliptic->o_grad,
				    o_Aq);
      }
    }
    
    if(elliptic->allNeumann) {
      o_tmp.copyTo(tmp);
      
      for(dlong n=0;n<Nblock;++n)
	alpha += tmp[n];
      
      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
      alphaG *= elliptic->allNeumannPenalty*elliptic->allNeumannScale*elliptic->allNeumannScale;
    }
    
    ellipticEndHaloExchange(elliptic, o_q, mesh->Np, recvBuffer);
    
    if(mesh->totalHaloPairs){
      offset = mesh->Nelements;
      if(options.compareArgs("BASIS", "NODAL")) {
	elliptic->partialGradientKernel(mesh->totalHaloPairs,
					offset,
					mesh->o_vgeo,
					mesh->o_Dmatrices,
					o_q,
					elliptic->o_grad);
      } else if(options.compareArgs("BASIS", "BERN")) {
	elliptic->partialGradientKernel(mesh->totalHaloPairs,
					offset,
					mesh->o_vgeo,
					mesh->o_D1ids,
					mesh->o_D2ids,
					mesh->o_D3ids,
					mesh->o_Dvals,
					o_q,
					elliptic->o_grad);
      }
    }
    
    if(mesh->NnotInternalElements) {
      if(options.compareArgs("BASIS", "NODAL")) {
	elliptic->partialIpdgKernel(mesh->NnotInternalElements,
				    mesh->o_notInternalElementIds,
				    mesh->o_vmapM,
				    mesh->o_vmapP,
				    lambda,
				    elliptic->tau,
				    mesh->o_vgeo,
				    mesh->o_sgeo,
				    elliptic->o_EToB,
				    mesh->o_Dmatrices,
				    mesh->o_LIFTT,
				    mesh->o_MM,
				    elliptic->o_grad,
				    o_Aq);
      } else if(options.compareArgs("BASIS", "BERN")) {
	elliptic->partialIpdgKernel(mesh->NnotInternalElements,
				    mesh->o_notInternalElementIds,
				    mesh->o_vmapM,
				    mesh->o_vmapP,
				    lambda,
				    elliptic->tau,
				    mesh->o_vgeo,
				    mesh->o_sgeo,
				    elliptic->o_EToB,
				    mesh->o_D1ids,
				    mesh->o_D2ids,
				    mesh->o_D3ids,
				    mesh->o_Dvals,
				    mesh->o_L0vals,
				    mesh->o_ELids,
				    mesh->o_ELvals,
				    mesh->o_BBMM,
				    elliptic->o_grad,
				    o_Aq);
      }
    }
    
    if(elliptic->allNeumann)
      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);
  } 
	
}


#endif

#if 1
// simplified version
void ellipticOperator(elliptic_t *elliptic, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *precision){
  
  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;
  const cgOptions_t &cgOptions = elliptic->cgOptions;

  dfloat *sendBuffer = elliptic->sendBuffer;
  dfloat *recvBuffer = elliptic->recvBuffer;
  dfloat *gradSendBuffer = elliptic->gradSendBuffer;
  dfloat *gradRecvBuffer = elliptic->gradRecvBuffer;

  dfloat alpha = 0., alphaG = 0.;

  if(cgOptions.continuous){
    if(cgOptions.serial && (elliptic->elementType == HEXAHEDRA))
      ellipticSerialAxHexKernel3D(mesh->Nq,  mesh->Nelements, mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq, elliptic->o_ggeoNoJW);
    else
      elliptic->AxKernel(mesh->Nelements, mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq);

    ogs_t *ellipticOgs = elliptic->ogs;
    int Nbytes = sizeof(dfloat);
#if 0    
    if(cgOptions.enableGatherScatters)
      ogsGatherScatter(o_Aq, ogsDfloat, ogsAdd, ellipticOgs);
#endif

    // Gather on device halo nodes
    occaGather(ellipticOgs->NhaloGather, ellipticOgs->o_haloGatherOffsets, ellipticOgs->o_haloGatherIds, dfloatString, "add", o_Aq, ogs::o_haloBuf);
    
    // copy device halo to host
    if(!cgOptions.serial){
      ogs::o_haloBuf.copyTo(ogs::haloBuf, ellipticOgs->NhaloGather*Nbytes, 0, "async: true");
      
      // host gather scatter using libgs 
      ogsHostGatherScatter(ogs::haloBuf, dfloatString, "add", ellipticOgs->haloGshSym);
      
      // copy totally gather halo data back from HOST to DEVICE
      ogs::o_haloBuf.copyFrom(ogs::haloBuf, ellipticOgs->NhaloGather*Nbytes, 0, "async: true");
    }
    else{
      // host gather scatter using libgs (directly on serial halo buffer)
      ogsHostGatherScatter(ogs::o_haloBuf.ptr(), dfloatString, "add", ellipticOgs->haloGshSym);
    }
    
    // do scatter back to local nodes
occaScatter(ellipticOgs->NhaloGather, ellipticOgs->o_haloGatherOffsets, ellipticOgs->o_haloGatherIds, dfloatString, "add", ogs::o_haloBuf, o_Aq);
    
    // have to do this after data from gslib gather scatter because owned nodes
    if(ellipticOgs->NlocalGather) {
      
      if(!cgOptions.serial)
	occaGatherScatter(ellipticOgs->NlocalGather, ellipticOgs->o_localGatherOffsets, ellipticOgs->o_localGatherIds, dfloatString, "add", o_Aq);
      else{
	const int Ngather = ellipticOgs->NlocalGather;
	const int * __restrict__ cpu_gatherStarts = (int*) ellipticOgs->o_localGatherOffsets.ptr();
	const int * __restrict__ cpu_gatherIds    = (int*) ellipticOgs->o_localGatherIds.ptr();
	dfloat * __restrict__ cpu_Aq = (dfloat*) o_Aq.ptr();

	cpu_Aq = (dfloat*)__builtin_assume_aligned(cpu_Aq, USE_OCCA_MEM_BYTE_ALIGN) ;
	cpu_gatherStarts = (int*)__builtin_assume_aligned(cpu_gatherStarts, USE_OCCA_MEM_BYTE_ALIGN) ;
	cpu_gatherIds = (int*)__builtin_assume_aligned(cpu_gatherIds, USE_OCCA_MEM_BYTE_ALIGN) ;
	
	for(dlong g=0;g<Ngather;++g){		
	  
	  const dlong start = cpu_gatherStarts[g];
	  const dlong end = cpu_gatherStarts[g+1];
	  if((start+1)!=end){
	    
	    double gAq = 0;
	    
	    for(dlong n=start;n<end;++n){
	      const dlong id = cpu_gatherIds[n];
	      gAq += cpu_Aq[id];
	    }
	  
	    for(dlong n=start;n<end;++n){
	      const dlong id = cpu_gatherIds[n];
	      cpu_Aq[id] = gAq;
	    }
	  }
	}	 
      }
    }

    if(elliptic->allNeumann) { // inspect this later
      
      dlong Nblock = elliptic->Nblock;
      dfloat *tmp = elliptic->tmp;
      occa::memory &o_tmp = elliptic->o_tmp;

      elliptic->innerProductKernel(mesh->Nelements*mesh->Np, elliptic->o_invDegree, o_q, o_tmp);
      o_tmp.copyTo(tmp);
      
      for(dlong n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
      alphaG *= elliptic->allNeumannPenalty*elliptic->allNeumannScale*elliptic->allNeumannScale;

      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);
    }

    //post-mask
    if (elliptic->Nmasked) 
      mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_Aq);

  } else if(options.compareArgs("DISCRETIZATION", "IPDG")) {
    printf("WARNING: DEBUGGING C0\n");
    MPI_Finalize();
    exit(-1);
  } 
  
}
#endif

#if 0
void ellipticOperator(elliptic_t *elliptic, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *precision){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  dfloat *sendBuffer = elliptic->sendBuffer;
  dfloat *recvBuffer = elliptic->recvBuffer;
  dfloat *gradSendBuffer = elliptic->gradSendBuffer;
  dfloat *gradRecvBuffer = elliptic->gradRecvBuffer;

  dfloat alpha = 0., alphaG = 0.;
  dlong Nblock = elliptic->Nblock;
  dfloat *tmp = elliptic->tmp;
  occa::memory &o_tmp = elliptic->o_tmp;

  int DEBUG_ENABLE_OGS = 1;
  options.getArgs("DEBUG ENABLE OGS", DEBUG_ENABLE_OGS);


  if(options.compareArgs("DISCRETIZATION", "CONTINUOUS")){
    ogs_t *ogs = elliptic->ogs;

#if 1
    int mapType = (elliptic->elementType==HEXAHEDRA &&
                   options.compareArgs("ELEMENT MAP", "TRILINEAR")) ? 1:0;

    int integrationType = (elliptic->elementType==HEXAHEDRA &&
			   options.compareArgs("ELLIPTIC INTEGRATION", "CUBATURE")) ? 1:0;
    
    occa::kernel &partialAxKernel = (strstr(precision, "float")) ? elliptic->partialFloatAxKernel : elliptic->partialAxKernel;
    
    if(mesh->NglobalGatherElements) {
      
      if(integrationType==0) { // GLL or non-hex
	if(mapType==0){
	  if(!options.compareArgs("THREAD MODEL", "Serial") || elliptic->elementType!=HEXAHEDRA)
	    partialAxKernel(mesh->NglobalGatherElements, mesh->o_globalGatherElementList,
			    mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq);
	  else
	    ellipticSerialPartialAxHexKernel3D(mesh->Nq, mesh->NglobalGatherElements, mesh->o_globalGatherElementList,
					       mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq);
	}
	else
	  partialAxKernel(mesh->NglobalGatherElements, mesh->o_globalGatherElementList,
			  elliptic->o_EXYZ, elliptic->o_gllzw, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq);
      }
      else{
	elliptic->partialCubatureAxKernel(mesh->NglobalGatherElements,
					  mesh->o_globalGatherElementList,
					  mesh->o_cubggeo,
					  mesh->o_cubD,
					  mesh->o_cubInterpT,
					  lambda, o_q, o_Aq);
      }
    }
#else

    elliptic->AxKernel(mesh->Nelements, mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq);
#endif
    if(DEBUG_ENABLE_OGS==1)
      ogsGatherScatterStart(o_Aq, ogsDfloat, ogsAdd, ogs);

#if 1
    if(mesh->NlocalGatherElements){
      if(integrationType==0) { // GLL or non-hex
	if(mapType==0){
	  if(!options.compareArgs("THREAD MODEL", "Serial") || elliptic->elementType!=HEXAHEDRA)
	    partialAxKernel(mesh->NlocalGatherElements, mesh->o_localGatherElementList,
			    mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq);
	  else
	    ellipticSerialPartialAxHexKernel3D(mesh->Nq, mesh->NlocalGatherElements, mesh->o_localGatherElementList,
					       mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq);
	}
	else
	  partialAxKernel(mesh->NlocalGatherElements, mesh->o_localGatherElementList,
			  elliptic->o_EXYZ, elliptic->o_gllzw, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq);
      }
      else{

	elliptic->partialCubatureAxKernel(mesh->NlocalGatherElements,
					  mesh->o_localGatherElementList,
					  mesh->o_cubggeo,
					  mesh->o_cubD,
					  mesh->o_cubInterpT,
					  lambda,
					  o_q,
					  o_Aq);
      }
    }
#endif

    // finalize gather using local and global contributions
    if(DEBUG_ENABLE_OGS==1)
      ogsGatherScatterFinish(o_Aq, ogsDfloat, ogsAdd, ogs);

    if(elliptic->allNeumann) {
      // mesh->sumKernel(mesh->Nelements*mesh->Np, o_q, o_tmp);
      elliptic->innerProductKernel(mesh->Nelements*mesh->Np, elliptic->o_invDegree, o_q, o_tmp);
      o_tmp.copyTo(tmp);

      for(dlong n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
      alphaG *= elliptic->allNeumannPenalty*elliptic->allNeumannScale*elliptic->allNeumannScale;

      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);
    }

    //post-mask
    if (elliptic->Nmasked) 
      mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_Aq);

  } else if(options.compareArgs("DISCRETIZATION", "IPDG")) {
    dlong offset = 0;
    dfloat alpha = 0., alphaG =0.;
    dlong Nblock = elliptic->Nblock;
    dfloat *tmp = elliptic->tmp;
    occa::memory &o_tmp = elliptic->o_tmp;

    ellipticStartHaloExchange(elliptic, o_q, mesh->Np, sendBuffer, recvBuffer);

    if(options.compareArgs("BASIS", "NODAL")) {
      elliptic->partialGradientKernel(mesh->Nelements,
				      offset,
				      mesh->o_vgeo,
				      mesh->o_Dmatrices,
				      o_q,
				      elliptic->o_grad);
    } else if(options.compareArgs("BASIS", "BERN")) {
      elliptic->partialGradientKernel(mesh->Nelements,
				      offset,
				      mesh->o_vgeo,
				      mesh->o_D1ids,
				      mesh->o_D2ids,
				      mesh->o_D3ids,
				      mesh->o_Dvals,
				      o_q,
				      elliptic->o_grad);
    }

    ellipticInterimHaloExchange(elliptic, o_q, mesh->Np, sendBuffer, recvBuffer);

    //Start the rank 1 augmentation if all BCs are Neumann
    //TODO this could probably be moved inside the Ax kernel for better performance
    if(elliptic->allNeumann)
      mesh->sumKernel(mesh->Nelements*mesh->Np, o_q, o_tmp);

    if(mesh->NinternalElements) {
      if(options.compareArgs("BASIS", "NODAL")) {
        elliptic->partialIpdgKernel(mesh->NinternalElements,
				    mesh->o_internalElementIds,
				    mesh->o_vmapM,
				    mesh->o_vmapP,
				    lambda,
				    elliptic->tau,
				    mesh->o_vgeo,
				    mesh->o_sgeo,
				    elliptic->o_EToB,
				    mesh->o_Dmatrices,
				    mesh->o_LIFTT,
				    mesh->o_MM,
				    elliptic->o_grad,
				    o_Aq);
      } else if(options.compareArgs("BASIS", "BERN")) {
        elliptic->partialIpdgKernel(mesh->NinternalElements,
				    mesh->o_internalElementIds,
				    mesh->o_vmapM,
				    mesh->o_vmapP,
				    lambda,
				    elliptic->tau,
				    mesh->o_vgeo,
				    mesh->o_sgeo,
				    elliptic->o_EToB,
				    mesh->o_D1ids,
				    mesh->o_D2ids,
				    mesh->o_D3ids,
				    mesh->o_Dvals,
				    mesh->o_L0vals,
				    mesh->o_ELids,
				    mesh->o_ELvals,
				    mesh->o_BBMM,
				    elliptic->o_grad,
				    o_Aq);
      }
    }

    if(elliptic->allNeumann) {
      o_tmp.copyTo(tmp);

      for(dlong n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
      alphaG *= elliptic->allNeumannPenalty*elliptic->allNeumannScale*elliptic->allNeumannScale;
    }

    ellipticEndHaloExchange(elliptic, o_q, mesh->Np, recvBuffer);

    if(mesh->totalHaloPairs){
      offset = mesh->Nelements;
      if(options.compareArgs("BASIS", "NODAL")) {
        elliptic->partialGradientKernel(mesh->totalHaloPairs,
					offset,
					mesh->o_vgeo,
					mesh->o_Dmatrices,
					o_q,
					elliptic->o_grad);
      } else if(options.compareArgs("BASIS", "BERN")) {
        elliptic->partialGradientKernel(mesh->totalHaloPairs,
					offset,
					mesh->o_vgeo,
					mesh->o_D1ids,
					mesh->o_D2ids,
					mesh->o_D3ids,
					mesh->o_Dvals,
					o_q,
					elliptic->o_grad);
      }
    }

    if(mesh->NnotInternalElements) {
      if(options.compareArgs("BASIS", "NODAL")) {
        elliptic->partialIpdgKernel(mesh->NnotInternalElements,
				    mesh->o_notInternalElementIds,
				    mesh->o_vmapM,
				    mesh->o_vmapP,
				    lambda,
				    elliptic->tau,
				    mesh->o_vgeo,
				    mesh->o_sgeo,
				    elliptic->o_EToB,
				    mesh->o_Dmatrices,
				    mesh->o_LIFTT,
				    mesh->o_MM,
				    elliptic->o_grad,
				    o_Aq);
      } else if(options.compareArgs("BASIS", "BERN")) {
        elliptic->partialIpdgKernel(mesh->NnotInternalElements,
				    mesh->o_notInternalElementIds,
				    mesh->o_vmapM,
				    mesh->o_vmapP,
				    lambda,
				    elliptic->tau,
				    mesh->o_vgeo,
				    mesh->o_sgeo,
				    elliptic->o_EToB,
				    mesh->o_D1ids,
				    mesh->o_D2ids,
				    mesh->o_D3ids,
				    mesh->o_Dvals,
				    mesh->o_L0vals,
				    mesh->o_ELids,
				    mesh->o_ELvals,
				    mesh->o_BBMM,
				    elliptic->o_grad,
				    o_Aq);
      }
    }

    if(elliptic->allNeumann)
      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);
  } 

}
#endif
