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

#include "adaptive.h"
#include "ogsInterface.h"

#include "omp.h"


//host gs
void adaptiveSerialOperator(adaptive_t *adaptive, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *precision){
  
  mesh_t *mesh = adaptive->mesh;
  setupAide &options = adaptive->options;

  int enableGatherScatters = 1;
  int enableReductions = 1;
  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int serial = options.compareArgs("THREAD MODEL", "Serial");
  int ipdg = options.compareArgs("DISCRETIZATION", "IPDG");
  
  options.getArgs("DEBUG ENABLE REDUCTIONS", enableReductions);
  options.getArgs("DEBUG ENABLE OGS", enableGatherScatters);
  
  //  printf("serialOperator: gathers = %d, reductions = %d, cts = %d, serial = %d, ipdg = %d\n",
  //	 enableGatherScatters, enableReductions, continuous, serial, ipdg);
  
  dfloat *sendBuffer = adaptive->sendBuffer;
  dfloat *recvBuffer = adaptive->recvBuffer;
  dfloat *gradSendBuffer = adaptive->gradSendBuffer;
  dfloat *gradRecvBuffer = adaptive->gradRecvBuffer;

  dfloat alpha = 0., alphaG = 0.;

  if(continuous && serial){

    adaptiveSerialAxHexKernel3D(mesh->Nq,  mesh->Nelements, mesh->o_ggeo, 
				mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq, adaptive->o_ggeoNoJW);
    
    ogs_t *adaptiveOgs = adaptive->ogs;
    
    ogsHostGatherScatter(o_Aq.ptr(), dfloatString, "add", adaptiveOgs->hostGsh);

#if USE_NULL_BOOST==1
    if(adaptive->allNeumann) { // inspect this later
      
      dlong Nblock = adaptive->Nblock;
      dfloat *tmp = adaptive->tmp;
      occa::memory &o_tmp = adaptive->o_tmp;

      adaptive->innerProductKernel(mesh->Nelements*mesh->Np, adaptive->o_invDegree, o_q, o_tmp);
      o_tmp.copyTo(tmp);
      
      for(dlong n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
      alphaG *= adaptive->allNeumannPenalty*adaptive->allNeumannScale*adaptive->allNeumannScale;

      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);
    }
#endif
    
    //post-mask
    if (adaptive->Nmasked) 
      mesh->maskKernel(adaptive->Nmasked, adaptive->o_maskIds, o_Aq);

  } else if(ipdg){
    printf("WARNING: DEBUGGING C0\n");
    MPI_Finalize();
    exit(-1);
  } 
  
}

void adaptiveOperator(adaptive_t *adaptive, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *precision){

  mesh_t *mesh = adaptive->mesh;
  setupAide &options = adaptive->options;

  int enableGatherScatters = 1;
  int enableReductions = 1;
  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int serial = options.compareArgs("THREAD MODEL", "Serial");
  int ipdg = options.compareArgs("DISCRETIZATION", "IPDG");
  
  options.getArgs("DEBUG ENABLE REDUCTIONS", enableReductions);
  options.getArgs("DEBUG ENABLE OGS", enableGatherScatters);

  //  printf("generalOperator: gathers = %d, reductions = %d, cts = %d, serial = %d, ipdg = %d\n",
  //	 enableGatherScatters, enableReductions, continuous, serial, ipdg);
  
  dfloat *sendBuffer = adaptive->sendBuffer;
  dfloat *recvBuffer = adaptive->recvBuffer;
  dfloat *gradSendBuffer = adaptive->gradSendBuffer;
  dfloat *gradRecvBuffer = adaptive->gradRecvBuffer;

  dfloat alpha = 0., alphaG = 0.;
  dlong Nblock = adaptive->Nblock;
  dfloat *tmp = adaptive->tmp;
  occa::memory &o_tmp = adaptive->o_tmp;

  if(continuous){

    // if in serial mode bypass occa
    if(serial && adaptive->elementType==HEXAHEDRA){
      adaptiveSerialOperator(adaptive, lambda, o_q, o_Aq, precision);
      return;
    }

    ogs_t *ogs = adaptive->ogs;

    if(mesh->NglobalGatherElements) {
      // do elements that touch partition boundary
      adaptive->partialAxKernel(mesh->NglobalGatherElements, mesh->o_globalGatherElementList,
				mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq);

    }
    
    if(enableGatherScatters)
      ogsGatherScatterStart(o_Aq, ogsDfloat, ogsAdd, ogs);
    

    if(mesh->NlocalGatherElements){
      // do elements that do not touch partition boundary
      adaptive->partialAxKernel(mesh->NlocalGatherElements, mesh->o_localGatherElementList,
				mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq);
    }

    // finalize gather using local and global contributions
    if(enableGatherScatters==1)
      ogsGatherScatterFinish(o_Aq, ogsDfloat, ogsAdd, ogs);

#if USE_NULL_BOOST==1
    if(adaptive->allNeumann) {

      adaptive->innerProductKernel(mesh->Nelements*mesh->Np, adaptive->o_invDegree, o_q, o_tmp);
      o_tmp.copyTo(tmp);
      
      for(dlong n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
      alphaG *= adaptive->allNeumannPenalty*adaptive->allNeumannScale*adaptive->allNeumannScale;

      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);
    }
#endif
    
    //post-mask
    if (adaptive->Nmasked) 
      mesh->maskKernel(adaptive->Nmasked, adaptive->o_maskIds, o_Aq);

  } else if(ipdg){
    dlong offset = 0;
    dfloat alpha = 0., alphaG =0.;
    dlong Nblock = adaptive->Nblock;
    dfloat *tmp = adaptive->tmp;
    occa::memory &o_tmp = adaptive->o_tmp;

    adaptiveStartHaloExchange(adaptive, o_q, mesh->Np, sendBuffer, recvBuffer);

    adaptive->partialGradientKernel(mesh->Nelements,
				    offset,
				    mesh->o_vgeo,
				    mesh->o_Dmatrices,
				    o_q,
				    adaptive->o_grad);

    adaptiveInterimHaloExchange(adaptive, o_q, mesh->Np, sendBuffer, recvBuffer);

    //Start the rank 1 augmentation if all BCs are Neumann
    //TODO this could probably be moved inside the Ax kernel for better performance
#if USE_NULL_BOOST==1
    if(adaptive->allNeumann)
      mesh->sumKernel(mesh->Nelements*mesh->Np, o_q, o_tmp);
#endif
    
    if(mesh->NinternalElements) {
      adaptive->partialIpdgKernel(mesh->NinternalElements,
				  mesh->o_internalElementIds,
				  mesh->o_vmapM,
				  mesh->o_vmapP,
				  lambda,
				  adaptive->tau,
				  mesh->o_vgeo,
				  mesh->o_sgeo,
				  adaptive->o_EToB,
				  mesh->o_Dmatrices,
				  mesh->o_LIFTT,
				  mesh->o_MM,
				  adaptive->o_grad,
				  o_Aq);
    }

#if USE_NULL_BOOST==1
    if(adaptive->allNeumann) {
      o_tmp.copyTo(tmp);

      for(dlong n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
      alphaG *= adaptive->allNeumannPenalty*adaptive->allNeumannScale*adaptive->allNeumannScale;
    }
#endif
    
    adaptiveEndHaloExchange(adaptive, o_q, mesh->Np, recvBuffer);

    if(mesh->totalHaloPairs){
      offset = mesh->Nelements;
      adaptive->partialGradientKernel(mesh->totalHaloPairs,
				      offset,
				      mesh->o_vgeo,
				      mesh->o_Dmatrices,
				      o_q,
				      adaptive->o_grad);
    }

    if(mesh->NnotInternalElements) {
      adaptive->partialIpdgKernel(mesh->NnotInternalElements,
				  mesh->o_notInternalElementIds,
				  mesh->o_vmapM,
				  mesh->o_vmapP,
				  lambda,
				  adaptive->tau,
				  mesh->o_vgeo,
				  mesh->o_sgeo,
				  adaptive->o_EToB,
				  mesh->o_Dmatrices,
				  mesh->o_LIFTT,
				  mesh->o_MM,
				  adaptive->o_grad,
				  o_Aq);
    }

#if USE_NULL_BOOST==1
    if(adaptive->allNeumann)
      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);
#endif
  } 

}

