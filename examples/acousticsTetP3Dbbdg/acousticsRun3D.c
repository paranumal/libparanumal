#include "acoustics3D.h"


void acousticsRun3Dbbdg(mesh3D *mesh){

  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->NpMax*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);
  
  //temp buffer to store the solution while plotting is happening.
  //Should probably do this better, but we have a Vandermond matrix to 
  // go from degree N modal space to degree NMax nodal space, but no inverse yet.
  dfloat *qtmp = (dfloat*) calloc(mesh->Nelements*mesh->NpMax*mesh->Nfields,sizeof(dfloat));


  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(iint rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat time = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

      if(mesh->totalHaloPairs>0){
        // extract halo node data
        meshHaloExtract(mesh, mesh->NpMax*mesh->Nfields*sizeof(dfloat),
              mesh->q, sendBuffer);
        
        // start halo exchange
        meshHaloExchangeStart(mesh,
              mesh->NpMax*mesh->Nfields*sizeof(dfloat),
              sendBuffer,
              recvBuffer);
      } 
      
      // compute volume contribution to DG acoustics RHS
      acousticsVolume3Dbbdg(mesh);

      if(mesh->totalHaloPairs>0){
        // wait for halo data to arrive
        meshHaloExchangeFinish(mesh);
    
        // copy data to halo zone of q
        memcpy(mesh->q+mesh->NpMax*mesh->Nfields*mesh->Nelements, recvBuffer, haloBytes);
      }

      // compute surface contribution to DG acoustics RHS
      acousticsSurface3Dbbdg(mesh, time);

      // update solution using LSERK4
      acousticsUpdate3D(mesh, mesh->rka[rk], mesh->rkb[rk]);
    }
    
    // estimate maximum error
    if((tstep%mesh->errorStep)==0) {
      //Save and transform to nodal basis
      for (iint e =0;e<mesh->Nelements;e++){
        iint id = e*mesh->NpMax*mesh->Nfields;
        iint N = mesh->N[e];
        
        for (iint n=0; n<mesh->Np[N]; n++){
          qtmp[id+n*mesh->Nfields+0] = mesh->q[id+n*mesh->Nfields+0];
          qtmp[id+n*mesh->Nfields+1] = mesh->q[id+n*mesh->Nfields+1];
          qtmp[id+n*mesh->Nfields+2] = mesh->q[id+n*mesh->Nfields+2];
          qtmp[id+n*mesh->Nfields+3] = mesh->q[id+n*mesh->Nfields+3];
        }
        for (iint n=0;n<mesh->NpMax;n++) {
          mesh->q[id+n*mesh->Nfields+0] = 0.0;
          mesh->q[id+n*mesh->Nfields+1] = 0.0;
          mesh->q[id+n*mesh->Nfields+2] = 0.0;
          mesh->q[id+n*mesh->Nfields+3] = 0.0;
        }
        for (iint n=0;n<mesh->NpMax;n++){
          for (iint m=0; m<mesh->Np[N]; m++){
            mesh->q[id+n*mesh->Nfields + 0] += mesh->VBplot[N][n*mesh->Np[N]+m]*qtmp[id+m*mesh->Nfields+0];
            mesh->q[id+n*mesh->Nfields + 1] += mesh->VBplot[N][n*mesh->Np[N]+m]*qtmp[id+m*mesh->Nfields+1];
            mesh->q[id+n*mesh->Nfields + 2] += mesh->VBplot[N][n*mesh->Np[N]+m]*qtmp[id+m*mesh->Nfields+2];
            mesh->q[id+n*mesh->Nfields + 3] += mesh->VBplot[N][n*mesh->Np[N]+m]*qtmp[id+m*mesh->Nfields+3];
          }
        }
      }
      
      // do error stuff on host
      acousticsError3D(mesh, mesh->dt*(tstep+1));

      // output field files
      iint fld = 3;
      meshPlotVTU3DP(mesh, "foo", fld);
      
      //Recover saved q
      for (iint e =0;e<mesh->Nelements;e++){
        iint id = e*mesh->NpMax*mesh->Nfields;
        iint N = mesh->N[e];

        for (iint n=0; n<mesh->Np[N]; n++){
          mesh->q[id+n*mesh->Nfields+0] = qtmp[id+n*mesh->Nfields+0];
          mesh->q[id+n*mesh->Nfields+1] = qtmp[id+n*mesh->Nfields+1];
          mesh->q[id+n*mesh->Nfields+2] = qtmp[id+n*mesh->Nfields+2];
          mesh->q[id+n*mesh->Nfields+3] = qtmp[id+n*mesh->Nfields+3];
        }
      }
    }
  }

  free(recvBuffer);
  free(sendBuffer);
  free(qtmp);
}

void acousticsOccaRun3Dbbdg(mesh3D *mesh){

  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->NpMax*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){

    for(iint rk=0;rk<mesh->Nrk;++rk){
      // intermediate stage time
      dfloat time = tstep*mesh->dt + mesh->dt*mesh->rkc[rk];

      if(mesh->totalHaloPairs>0){
        // extract halo on DEVICE
        iint Nentries = mesh->NpMax*mesh->Nfields;
        
        mesh->haloExtractKernel(mesh->totalHaloPairs,
              Nentries,
              mesh->o_haloElementList,
              mesh->o_q,
              mesh->o_haloBuffer);
        
        // copy extracted halo to HOST 
        mesh->o_haloBuffer.copyTo(sendBuffer);      
        
        // start halo exchange
        meshHaloExchangeStart(mesh,
              mesh->NpMax*mesh->Nfields*sizeof(dfloat),
              sendBuffer,
              recvBuffer);
      }
          
      for (iint p=1;p<=mesh->NMax;p++) {
        // compute volume contribution to DG acoustics RHS
        if (mesh->NelOrder[p])
          mesh->volumeKernel[p](mesh->NelOrder[p],
                mesh->o_NelList[p],
                mesh->o_vgeo,
                mesh->o_D0ids[p],
                mesh->o_D1ids[p],
                mesh->o_D2ids[p],
                mesh->o_D3ids[p],
                mesh->o_Dvals[p],
                mesh->o_q,
                mesh->o_rhsq);
      }
          
      if(mesh->totalHaloPairs>0){
      // wait for halo data to arrive
      meshHaloExchangeFinish(mesh);
      
      // copy halo data to DEVICE
      size_t offset = mesh->NpMax*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      mesh->o_q.copyFrom(recvBuffer, haloBytes, offset);
      }
          
      for (iint p=1;p<=mesh->NMax;p++) {
        // compute surface contribution to DG acoustics RHS
        if (mesh->NelOrder[p])
          mesh->surfaceKernel[p](mesh->NelOrder[p],
                  mesh->o_NelList[p],
                  mesh->o_N,
                  mesh->o_sgeo,
                  mesh->o_L0ids[p],
                  mesh->o_L0vals[p],
                  mesh->o_ELids[p],
                  mesh->o_ELvals[p],
                  mesh->o_BBLower[p],
                  mesh->o_BBRaiseids[p],
                  mesh->o_BBRaiseVals[p],
                  mesh->o_vmapM,
                  mesh->o_vmapP,
                  mesh->o_EToE,
                  mesh->o_EToF,
                  mesh->o_EToB,
                  time,
                  mesh->o_x,
                  mesh->o_y,
                  mesh->o_z,
                  mesh->o_q,
                  mesh->o_rhsq);
      }
      
      // update solution using Runge-Kutta
      mesh->updateKernel(mesh->Nelements*mesh->NpMax*mesh->Nfields,
       mesh->dt,
       mesh->rka[rk],
       mesh->rkb[rk],
       mesh->o_rhsq,
       mesh->o_resq,
       mesh->o_q);
      
    }
    
    // estimate maximum error
    if((tstep%mesh->errorStep)==0){

      // copy data back to host
      mesh->o_q.copyTo(mesh->q);

      //Transform to nodal basis
      dfloat qtmp[mesh->Nfields*mesh->NpMax];
      for (iint e =0;e<mesh->Nelements;e++){
        iint id = e*mesh->NpMax*mesh->Nfields;
        iint N = mesh->N[e];

        for (iint n=0; n<mesh->Np[N]; n++){
          qtmp[n*mesh->Nfields + 0] = mesh->q[id+n*mesh->Nfields+0];
          qtmp[n*mesh->Nfields + 1] = mesh->q[id+n*mesh->Nfields+1];
          qtmp[n*mesh->Nfields + 2] = mesh->q[id+n*mesh->Nfields+2];
          qtmp[n*mesh->Nfields + 3] = mesh->q[id+n*mesh->Nfields+3];
        }
        for (iint n=0;n<mesh->NpMax;n++) {
          mesh->q[id+n*mesh->Nfields+0] = 0.0;
          mesh->q[id+n*mesh->Nfields+1] = 0.0;
          mesh->q[id+n*mesh->Nfields+2] = 0.0;
          mesh->q[id+n*mesh->Nfields+3] = 0.0;
        }
        for (iint n=0;n<mesh->NpMax;n++){
          for (iint m=0; m<mesh->Np[N]; m++){
            mesh->q[id+n*mesh->Nfields + 0] += mesh->VBplot[N][n*mesh->Np[N]+m]*qtmp[m*mesh->Nfields+0];
            mesh->q[id+n*mesh->Nfields + 1] += mesh->VBplot[N][n*mesh->Np[N]+m]*qtmp[m*mesh->Nfields+1];
            mesh->q[id+n*mesh->Nfields + 2] += mesh->VBplot[N][n*mesh->Np[N]+m]*qtmp[m*mesh->Nfields+2];
            mesh->q[id+n*mesh->Nfields + 3] += mesh->VBplot[N][n*mesh->Np[N]+m]*qtmp[m*mesh->Nfields+3];
          }
        }
      }

      // do error stuff on host
      acousticsError3D(mesh, mesh->dt*(tstep+1));
      
      // output field files
      iint fld = 3;
      meshPlotVTU3DP(mesh, "foo", fld);
    }
  }

  free(recvBuffer);
  free(sendBuffer);
}
