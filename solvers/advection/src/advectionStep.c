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

#include "advection.h"

void advectionDopriStep(advection_t *advection, setupAide &newOptions, const dfloat time){

  mesh_t *mesh = advection->mesh;
  
  //RK step
  for(int rk=0;rk<advection->Nrk;++rk){
    
    // t_rk = t + C_rk*dt
    dfloat currentTime = time + advection->rkC[rk]*mesh->dt;
    
    //compute RK stage 
    // rkq = q + dt sum_{i=0}^{rk-1} a_{rk,i}*rhsq_i
    advection->rkStageKernel(mesh->Nelements,
			     rk,
			     mesh->dt,
			     advection->o_rkA,
			     advection->o_q,
			     advection->o_rkrhsq,
			     advection->o_rkq);
    
    //compute RHS
    // rhsq = F(currentTIme, rkq)
    // extract q halo on DEVICE
    if(mesh->totalHaloPairs>0){
      int Nentries = mesh->Np*advection->Nfields;          
      mesh->haloExtractKernel(mesh->totalHaloPairs,
			      Nentries,
			      mesh->o_haloElementList,
			      advection->o_rkq,
			      advection->o_haloBuffer);
      

      mesh->device.finish();
      mesh->device.setStream(mesh->dataStream);
      
      // copy extracted halo to HOST 
      advection->o_haloBuffer.copyTo(advection->sendBuffer);      

      mesh->device.setStream(mesh->defaultStream);
      
    }

    if(newOptions.compareArgs("ADVECTION FORMULATION", "NODAL")){
      advection->volumeKernel(mesh->Nelements, 
			      mesh->o_vgeo, 
			      mesh->o_Dmatrices,
			      advection->o_advectionVelocityJW,
			      advection->o_rkq, 
			      advection->o_rhsq);
    }
    
    if(newOptions.compareArgs("ADVECTION FORMULATION", "CUBATURE")){
      advection->volumeKernel(mesh->Nelements, 
			      mesh->o_vgeo,
			      mesh->o_cubvgeo, 
			      mesh->o_cubDWmatrices,
			      mesh->o_cubInterpT,
			      mesh->o_cubProjectT,
			      advection->o_cubAdvectionVelocityJW,
			      advection->o_rkq, 
			      advection->o_rhsq);
    }

    
    // wait for q halo data to arrive
    if(mesh->totalHaloPairs>0){
      mesh->device.setStream(mesh->dataStream);
      mesh->device.finish();
      
      // start halo exchange
      meshHaloExchangeStart(mesh, mesh->Np*advection->Nfields*sizeof(dfloat), advection->sendBuffer, advection->recvBuffer);

      meshHaloExchangeFinish(mesh);
      
      // copy halo data to DEVICE
      size_t offset = mesh->Np*advection->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      advection->o_rkq.copyFrom(advection->recvBuffer, advection->haloBytes, offset, "async: true");

      mesh->device.finish(); // finish copy to device
      mesh->device.setStream(mesh->defaultStream);
    }
    
    advection->surfaceKernel(mesh->Nelements, 
			     mesh->o_sgeo, 
			     mesh->o_LIFTT, 
			     mesh->o_vmapM, 
			     mesh->o_vmapP, 
			     mesh->o_EToB,
			     currentTime, 
			     mesh->o_x, 
			     mesh->o_y,
			     mesh->o_z,
			     advection->o_advectionVelocityM,
			     advection->o_advectionVelocityP,
			     advection->o_rkq, 
			     advection->o_rhsq);
    
    // update solution using Runge-Kutta
    // rkrhsq_rk = rhsq
    // if rk==6 
    //   q = q + dt*sum_{i=0}^{rk} rkA_{rk,i}*rkrhs_i
    //   rkerr = dt*sum_{i=0}^{rk} rkE_{rk,i}*rkrhs_i
    advection->rkUpdateKernel(mesh->Nelements, 
			      rk,
			      mesh->dt, 
			      advection->o_rkA, 
			      advection->o_rkE, 
			      advection->o_q,
			      advection->o_rhsq, 
			      advection->o_rkrhsq, 
			      advection->o_rkq,
			      advection->o_rkerr);
  }
}


void advectionLserkStep(advection_t *advection, setupAide &newOptions, const dfloat time){

  mesh_t *mesh = advection->mesh;

  static dfloat *hostSendBuffer = NULL;

  if(!hostSendBuffer)
    hostSendBuffer = (dfloat*) calloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields, sizeof(dfloat));

  static dfloat *hostRecvBuffer = NULL;

  if(!hostRecvBuffer)
    hostRecvBuffer = (dfloat*) calloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields, sizeof(dfloat));
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(int rk=0;rk<mesh->Nrk;++rk){
      
    dfloat currentTime = time + mesh->rkc[rk]*mesh->dt;

    if(!newOptions.compareArgs("ADVECTION FORMULATION", "COMBINED")){

      // extract q halo on DEVICE
      if(mesh->totalHaloPairs>0){
	int Nentries = mesh->Np*advection->Nfields;
        
	mesh->haloExtractKernel(mesh->totalHaloPairs,
				Nentries,
				mesh->o_haloElementList,
				advection->o_q,
				advection->o_haloBuffer);

	mesh->device.finish();
	mesh->device.setStream(mesh->dataStream);

        // copy extracted halo to HOST
        advection->o_haloBuffer.copyTo(advection->sendBuffer,"async: true");
        mesh->device.setStream(mesh->defaultStream);
	
      }

      if(newOptions.compareArgs("ADVECTION FORMULATION", "NODAL")){
	advection->volumeKernel(mesh->Nelements, 
				mesh->o_vgeo, 
				mesh->o_Dmatrices,
				advection->o_advectionVelocityJW,
				advection->o_q, 
				advection->o_rhsq);
      }

      if(newOptions.compareArgs("ADVECTION FORMULATION", "CUBATURE")){
	advection->volumeKernel(mesh->Nelements, 
				mesh->o_vgeo,
				mesh->o_cubvgeo, 
				mesh->o_cubDWmatrices,
				mesh->o_cubInterpT,
				mesh->o_cubProjectT,
				advection->o_cubAdvectionVelocityJW,
				advection->o_q, 
				advection->o_rhsq);
      }
      
      // wait for q halo data to arrive
      if(mesh->totalHaloPairs>0){

	mesh->device.setStream(mesh->dataStream);
	mesh->device.finish(); //finish copy to host

#if 0
	memcpy(hostSendBuffer, advection->sendBuffer,
	       mesh->totalHaloPairs*mesh->Np*advection->Nfields*sizeof(dfloat));
	
	// halo exchange
        meshHaloExchangeStart(mesh, mesh->Np*advection->Nfields*sizeof(dfloat),
			      hostSendBuffer,
			      hostRecvBuffer);
	meshHaloExchangeFinish(mesh);

	memcpy(advection->recvBuffer, hostRecvBuffer,
	       mesh->totalHaloPairs*mesh->Np*advection->Nfields*sizeof(dfloat));
#else
	// halo exchange
        meshHaloExchangeStart(mesh, mesh->Np*advection->Nfields*sizeof(dfloat),
			      advection->sendBuffer,
			      advection->recvBuffer);
	meshHaloExchangeFinish(mesh);
#endif
	
        // copy halo data to DEVICE
        size_t offset = mesh->Np*advection->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
	advection->o_q.copyFrom(advection->recvBuffer, advection->haloBytes, offset, "async: true");
      }
      
      advection->surfaceKernel(mesh->Nelements, 
			       mesh->o_sgeo, 
			       mesh->o_LIFTT, 
			       mesh->o_vmapM, 
			       mesh->o_vmapP, 
			       mesh->o_EToB,
			       currentTime, 
			       mesh->o_x, 
			       mesh->o_y,
			       mesh->o_z,
			       advection->o_advectionVelocityM,
			       advection->o_advectionVelocityP,
			       advection->o_q, 
			       advection->o_rhsq);
      
      // update solution using Runge-Kutta
      advection->updateKernel(mesh->Nelements, 
			      mesh->dt, 
			      mesh->rka[rk], 
			      mesh->rkb[rk], 
			      advection->o_rhsq, 
			      advection->o_resq, 
			      advection->o_q);
    }
    else{

      occa::memory &o_sourceq = (!(rk%2)) ? advection->o_q: advection->o_rhsq;
      occa::memory &o_destq   = (!(rk%2)) ? advection->o_rhsq: advection->o_q;

      //      occa::streamTag haloCopyToHostEvent;

      //      mesh->device.setStream(mesh->computeStream);         

      // extract q halo on DEVICE
      if(mesh->totalHaloPairs>0){
        mesh->device.setStream(mesh->dataStream);

	// NFP 
        int Nentries = mesh->Nfp*advection->Nfields;

        mesh->haloExtractKernel(mesh->totalHaloPairs,
                                Nentries,
                                mesh->o_haloElementList,
                                o_sourceq,
                                advection->o_haloBuffer);

        // copy extracted halo to HOST using async copy
        advection->o_haloBuffer.copyTo(advection->sendBuffer, "async: true");
	
	mesh->device.setStream(mesh->defaultStream);
      }

      if(mesh->NinternalElements>0){
	
	// NODAL only
	advection->volumeKernel(mesh->NinternalElements,
				mesh->o_internalElementIds,
				mesh->dt,
				mesh->rka[rk], 
				mesh->rkb[rk], 
				mesh->o_vgeo, 
				mesh->o_Dmatrices,
				advection->o_advectionVelocityJW,
				mesh->o_vmapM,
				mesh->o_vmapP,
				advection->o_advectionVelocityM,
				advection->o_advectionVelocityP,
				advection->o_resq,
				o_sourceq, 
				o_destq);
      }

      
      if(mesh->totalHaloPairs>0){

	mesh->device.setStream(mesh->dataStream);

	mesh->device.finish(); //finish halo extract

#if 0
	memcpy(hostSendBuffer, advection->sendBuffer,
	       mesh->totalHaloPairs*mesh->Nfp*advection->Nfields*sizeof(dfloat));
	
	// halo exchange (NFP)
        meshHaloExchangeStart(mesh, mesh->Nfp*advection->Nfields*sizeof(dfloat),
			      hostSendBuffer,
			      hostRecvBuffer);
	meshHaloExchangeFinish(mesh);

	memcpy(advection->recvBuffer, hostRecvBuffer,
	       mesh->totalHaloPairs*mesh->Nfp*advection->Nfields*sizeof(dfloat));
#else
        // start MPI halo exchange
        meshHaloExchangeStart(mesh, mesh->Nfp*advection->Nfields*sizeof(dfloat), advection->sendBuffer, advection->recvBuffer);
        meshHaloExchangeFinish(mesh);
#endif
        // copy halo data to DEVICE
        size_t offset = mesh->Np*advection->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data

	hlong faceBytes = mesh->totalHaloPairs*mesh->Nfp*sizeof(dfloat); // NFP !!
        o_sourceq.copyFrom(advection->recvBuffer, faceBytes, offset, "async: true");

        mesh->device.finish(); //finish halo extract
        mesh->device.setStream(mesh->defaultStream);

        if(mesh->NnotInternalElements){
          advection->volumeKernel(mesh->NnotInternalElements,
                                  mesh->o_notInternalElementIds,
                                  mesh->dt,
                                  mesh->rka[rk],
                                  mesh->rkb[rk],
                                  mesh->o_vgeo,
                                  mesh->o_Dmatrices,
                                  advection->o_advectionVelocityJW,
                                  mesh->o_vmapM,
                                  mesh->o_vmapP,
                                  advection->o_advectionVelocityM,
                                  advection->o_advectionVelocityP,
                                  advection->o_resq,
                                  o_sourceq,
                                  o_destq);
	  mesh->device.finish(); //finish halo extract
	}
      }


      
      if(rk==4){
	//  TW ?
	//	mesh->device.finish();
	o_destq.copyTo(o_sourceq);
      }
    }
  }

  
  
}
