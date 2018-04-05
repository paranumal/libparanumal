#include "cnsQuad2D.h"

void cnsRunQuad2D(cns_t *cns, char *options){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *mesh = cns->mesh;
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(int tstep=0;tstep<mesh->NtimeSteps;++tstep){

    int advSwitch = 1;//(tstep>100);
    
    for(int rk=0;rk<mesh->Nrk;++rk){

      dfloat currentTime = tstep*mesh->dt + mesh->rkc[rk]*mesh->dt;
      
      // extract q halo on DEVICE
      if(mesh->totalHaloPairs>0){
        int Nentries = mesh->Np*cns->Nfields;
        
        mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_q, cns->o_haloBuffer);
        
        // copy extracted halo to HOST 
        mesh->o_haloBuffer.copyTo(cns->sendBuffer);      
        
        // start halo exchange
        meshHaloExchangeStart(mesh, mesh->Np*cns->Nfields*sizeof(dfloat), cns->sendBuffer, cns->recvBuffer);
      }

      // now compute viscous stresses
      cns->stressesVolumeKernel(mesh->Nelements, 
                                mesh->o_vgeo, 
                                mesh->o_D, 
                                cns->mu, 
                                cns->o_q, 
                                cns->o_viscousStresses);

      // wait for q halo data to arrive
      if(mesh->totalHaloPairs>0){
        meshHaloExchangeFinish(mesh);
        
        // copy halo data to DEVICE
        size_t offset = mesh->Np*cns->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
        cns->o_q.copyFrom(cns->recvBuffer, cns->haloBytes, offset);
      }
      
      cns->stressesSurfaceKernel(mesh->Nelements, 
                                 mesh->o_sgeo, 
                                 cns->o_LIFTT,
                                 mesh->o_vmapM, 
                                 mesh->o_vmapP, 
                                 mesh->o_EToB, 
                                 currentTime,
                                 mesh->o_x, 
                                 mesh->o_y, 
                                 cns->mu, 
                                 cns->o_q, 
                                 cns->o_viscousStresses);
      
      // extract stresses halo on DEVICE
      if(mesh->totalHaloPairs>0){
        int Nentries = mesh->Np*cns->Nstresses;
        
        mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_viscousStresses, cns->o_haloStressesBuffer);
        
        // copy extracted halo to HOST 
        mesh->o_haloBuffer.copyTo(cns->sendStressesBuffer);      
        
        // start halo exchange
        meshHaloExchangeStart(mesh, mesh->Np*cns->Nstresses*sizeof(dfloat), cns->sendStressesBuffer, cns->recvStressesBuffer);
      }
      
      // compute volume contribution to DG cns RHS
      if (strstr(options,"CUBATURE")) {
        cns->cubatureVolumeKernel(mesh->Nelements, 
                                  advSwitch, 
                                  mesh->o_vgeo, 
                                  mesh->o_cubvgeo,
                                  mesh->o_cubDWT,
                                  mesh->o_cubInterpT,
                                  mesh->o_cubProjectT,
                                  cns->o_viscousStresses, 
                                  cns->o_q, 
                                  cns->o_rhsq);
      } else {
        cns->volumeKernel(mesh->Nelements, 
                          advSwitch, 
                          mesh->o_vgeo, 
                          mesh->o_D, 
                          cns->o_viscousStresses, 
                          cns->o_q, 
                          cns->o_rhsq);
      }

      // wait for halo stresses data to arrive
      if(mesh->totalHaloPairs>0){
        meshHaloExchangeFinish(mesh);
        
        // copy halo data to DEVICE
        size_t offset = mesh->Np*cns->Nstresses*mesh->Nelements*sizeof(dfloat); // offset for halo data
        cns->o_viscousStresses.copyFrom(cns->recvStressesBuffer, cns->haloStressesBytes, offset);
      }
      
      // compute surface contribution to DG cns RHS (LIFTT ?)
      if (strstr(options,"CUBATURE")) {
        cns->cubatureSurfaceKernel(mesh->Nelements, 
                                   advSwitch, 
                                   mesh->o_vgeo, 
                                   mesh->o_cubsgeo, 
                                   mesh->o_vmapM, 
                                   mesh->o_vmapP, 
                                   mesh->o_EToB,
                                   mesh->o_cubInterpT,
                                   mesh->o_cubProjectT,
                                   currentTime, 
                                   mesh->o_intx, 
                                   mesh->o_inty, 
                                   cns->mu, 
                                   cns->o_q, 
                                   cns->o_viscousStresses, 
                                   cns->o_rhsq);
      } else {
        cns->surfaceKernel(mesh->Nelements, 
                           advSwitch, 
                           mesh->o_sgeo, 
                           cns->o_LIFTT, 
                           mesh->o_vmapM, 
                           mesh->o_vmapP, 
                           mesh->o_EToB,
                           currentTime, 
                           mesh->o_x, 
                           mesh->o_y, 
                           cns->mu, 
                           cns->o_q, 
                           cns->o_viscousStresses, 
                           cns->o_rhsq);
      }
      
      // update solution using Runge-Kutta
      cns->updateKernel(mesh->Nelements, mesh->dt, mesh->rka[rk], mesh->rkb[rk], cns->o_rhsq, cns->o_resq, cns->o_q);
      
    }
    
    if(((tstep+1)%mesh->errorStep)==0){
      cnsReportQuad2D(cns, tstep+1, options);
    }
  }
}
