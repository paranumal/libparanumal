#include "cnsTri2D.h"

void cnsRunTri2D(cns_t *cns, char *options){

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
        cns->o_haloBuffer.copyTo(cns->sendBuffer);      
        
        // start halo exchange
        meshHaloExchangeStart(mesh, mesh->Np*cns->Nfields*sizeof(dfloat), cns->sendBuffer, cns->recvBuffer);
      }

      // now compute viscous stresses
      cns->stressesVolumeKernel(mesh->Nelements, 
                                mesh->o_vgeo, 
                                mesh->o_DrT, 
                                mesh->o_DsT, 
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
                                 mesh->o_LIFTT,
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
        cns->o_haloStressesBuffer.copyTo(cns->sendStressesBuffer);      
        
        // start halo exchange
        meshHaloExchangeStart(mesh, mesh->Np*cns->Nstresses*sizeof(dfloat), cns->sendStressesBuffer, cns->recvStressesBuffer);
      }
      
      // compute volume contribution to DG cns RHS
      if (strstr(options,"CUBATURE")) {
        cns->cubatureVolumeKernel(mesh->Nelements, 
                                  advSwitch, 
                                  mesh->o_vgeo, 
                                  mesh->o_cubDrWT,
                                  mesh->o_cubDsWT,
                                  mesh->o_cubInterpT,
                                  cns->o_viscousStresses, 
                                  cns->o_q, 
                                  cns->o_rhsq);
      } else {
        cns->volumeKernel(mesh->Nelements, 
                          advSwitch, 
                          mesh->o_vgeo, 
                          mesh->o_DrT,
                          mesh->o_DsT,
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
                                   mesh->o_sgeo, 
                                   mesh->o_intInterpT,
                                   mesh->o_intLIFTT, 
                                   mesh->o_vmapM, 
                                   mesh->o_vmapP, 
                                   mesh->o_EToB,
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
                           mesh->o_LIFTT, 
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
      cns->updateKernel(mesh->Nelements, 
                        mesh->dt, 
                        mesh->rka[rk], 
                        mesh->rkb[rk], 
                        cns->o_rhsq, 
                        cns->o_resq, 
                        cns->o_q);
    }
    
    if(((tstep+1)%mesh->errorStep)==0){
      cnsReportTri2D(cns, tstep+1, options);
    }
  }
}
