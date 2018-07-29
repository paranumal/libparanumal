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

#include "acoustics.h"

void acousticsDopriStep(acoustics_t *acoustics, setupAide &newOptions, const dfloat time){

  mesh_t *mesh = acoustics->mesh;
  
  //RK step
  for(int rk=0;rk<acoustics->Nrk;++rk){
    
    // t_rk = t + C_rk*dt
    dfloat currentTime = time + acoustics->rkC[rk]*mesh->dt;
    
    //compute RK stage 
    // rkq = q + dt sum_{i=0}^{rk-1} a_{rk,i}*rhsq_i
    acoustics->rkStageKernel(mesh->Nelements,
		       rk,
		       mesh->dt,
		       acoustics->o_rkA,
		       acoustics->o_q,
		       acoustics->o_rkrhsq,
		       acoustics->o_rkq);
    
    //compute RHS
    // rhsq = F(currentTIme, rkq)

    // extract q halo on DEVICE
    if(mesh->totalHaloPairs>0){
      int Nentries = mesh->Np*acoustics->Nfields;          
      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, acoustics->o_rkq, acoustics->o_haloBuffer);
      
      // copy extracted halo to HOST 
      acoustics->o_haloBuffer.copyTo(acoustics->sendBuffer);      
      
      // start halo exchange
      meshHaloExchangeStart(mesh, mesh->Np*acoustics->Nfields*sizeof(dfloat), acoustics->sendBuffer, acoustics->recvBuffer);
    }

    acoustics->volumeKernel(mesh->Nelements, 
		      mesh->o_vgeo, 
		      mesh->o_Dmatrices,
		      acoustics->o_rkq, 
		      acoustics->o_rhsq);

    // wait for q halo data to arrive
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);
          
      // copy halo data to DEVICE
      size_t offset = mesh->Np*acoustics->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      acoustics->o_rkq.copyFrom(acoustics->recvBuffer, acoustics->haloBytes, offset);
    }

    acoustics->surfaceKernel(mesh->Nelements, 
			     mesh->o_sgeo, 
			     mesh->o_LIFTT, 
			     mesh->o_vmapM, 
			     mesh->o_vmapP, 
			     mesh->o_EToB,
			     currentTime, 
			     mesh->o_x, 
			     mesh->o_y,
			     mesh->o_z, 
			     acoustics->o_rkq, 
			     acoustics->o_rhsq);
    
    // update solution using Runge-Kutta
    // rkrhsq_rk = rhsq
    // if rk==6 
    //   q = q + dt*sum_{i=0}^{rk} rkA_{rk,i}*rkrhs_i
    //   rkerr = dt*sum_{i=0}^{rk} rkE_{rk,i}*rkrhs_i
    acoustics->rkUpdateKernel(mesh->Nelements, 
			rk,
			mesh->dt, 
			acoustics->o_rkA, 
			acoustics->o_rkE, 
			acoustics->o_q,
			acoustics->o_rhsq, 
			acoustics->o_rkrhsq, 
			acoustics->o_rkq,
			acoustics->o_rkerr);
  }
}


void acousticsLserkStep(acoustics_t *acoustics, setupAide &newOptions, const dfloat time){

  mesh_t *mesh = acoustics->mesh;
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(int rk=0;rk<mesh->Nrk;++rk){
      
    dfloat currentTime = time + mesh->rkc[rk]*mesh->dt;
      
    // extract q halo on DEVICE
    if(mesh->totalHaloPairs>0){
      int Nentries = mesh->Np*acoustics->Nfields;
        
      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, acoustics->o_q, acoustics->o_haloBuffer);
        
      // copy extracted halo to HOST 
      acoustics->o_haloBuffer.copyTo(acoustics->sendBuffer);      
        
      // start halo exchange
      meshHaloExchangeStart(mesh, mesh->Np*acoustics->Nfields*sizeof(dfloat), acoustics->sendBuffer, acoustics->recvBuffer);
    }

    acoustics->volumeKernel(mesh->Nelements, 
		      mesh->o_vgeo, 
		      mesh->o_Dmatrices,
		      acoustics->o_q, 
		      acoustics->o_rhsq);
    
    
    // wait for q halo data to arrive
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);
      
      // copy halo data to DEVICE
      size_t offset = mesh->Np*acoustics->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      acoustics->o_q.copyFrom(acoustics->recvBuffer, acoustics->haloBytes, offset);
    }

    acoustics->surfaceKernel(mesh->Nelements, 
		       mesh->o_sgeo, 
		       mesh->o_LIFTT, 
		       mesh->o_vmapM, 
		       mesh->o_vmapP, 
		       mesh->o_EToB,
		       currentTime, 
		       mesh->o_x, 
		       mesh->o_y,
		       mesh->o_z, 
		       acoustics->o_q, 
		       acoustics->o_rhsq);
        
    // update solution using Runge-Kutta
    acoustics->updateKernel(mesh->Nelements, 
		      mesh->dt, 
		      mesh->rka[rk], 
		      mesh->rkb[rk], 
		      acoustics->o_rhsq, 
		      acoustics->o_resq, 
		      acoustics->o_q);
  }
}
