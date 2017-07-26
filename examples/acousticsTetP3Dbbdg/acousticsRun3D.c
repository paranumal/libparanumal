#include "acoustics3D.h"

void addSourceField(mesh3D *mesh, dfloat *q, dfloat t);

void acousticsRun3Dbbdg(mesh3D *mesh){

  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->NfpMax*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);
  
  //temp buffer to store the solution while plotting is happening.
  //Should probably do this better, but we have a Vandermond matrix to 
  // go from degree N modal space to degree NMax nodal space, but no inverse yet.
  dfloat *qtmp = (dfloat*) calloc(mesh->Nelements*mesh->NpMax*mesh->Nfields,sizeof(dfloat));

  int Nframe =0;

  //populate the trace buffer fQ
  for (iint l=0;l<mesh->MRABNlevels;l++) {
    acousticsMRABpmlUpdate3D(mesh, 0., 0., 0., l, 0.);
    acousticsMRABUpdate3D(mesh, 0., 0., 0., l, 0., 0.);
  }


  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){ 
    for (iint Ntick=0; Ntick < pow(2,mesh->MRABNlevels-1);Ntick++) {

      // intermediate stage time
      dfloat t = mesh->dt*(tstep*pow(2,mesh->MRABNlevels-1) + Ntick);

      iint lev;
      for (lev=0;lev<mesh->MRABNlevels;lev++) 
        if (Ntick % (1<<lev) != 0) break; //find the max lev to compute rhs

      if(mesh->totalHaloPairs>0){
        // extract halo node data
        meshHaloExtract(mesh, mesh->Nfields*mesh->NfpMax*mesh->Nfaces*sizeof(dfloat),
              mesh->fQP, sendBuffer);
        
        // start halo exchange
        meshHaloExchangeStart(mesh,
              mesh->Nfields*mesh->NfpMax*mesh->Nfaces*sizeof(dfloat),
              sendBuffer,
              recvBuffer);
      }     

      // compute volume contribution to DG acoustics RHS
      for (iint l=0;l<lev;l++) {
        acousticsPmlVolume3Dbbdg(mesh,l);
        acousticsVolume3Dbbdg(mesh,l);
      }

      if(mesh->totalHaloPairs>0){
        // wait for halo data to arrive
        meshHaloExchangeFinish(mesh);
    
        // copy data to halo zone of q
        memcpy(mesh->fQP+mesh->NfpMax*mesh->Nfields*mesh->Nelements*mesh->Nfaces, recvBuffer, haloBytes);
      }
      
      // compute surface contribution to DG acoustics RHS
      for (iint l=0;l<lev;l++) {
        acousticsPmlSurface3Dbbdg(mesh,l,t);
        acousticsSurface3Dbbdg(mesh,l,t);
      }
      
      dfloat a1, a2, a3;
      dfloat b1, b2, b3;
      if (tstep==0) {
        a1 = 1.0;
        a2 = 0.0;
        a3 = 0.0;
        b1 = 0.5;
        b2 = 0.0;
        b3 = 0.0; 
      } else if (tstep ==1) {
        a1 =  3.0/2.0;
        a2 = -1.0/2.0;
        a3 =  0.0;
        b1 =  5.0/8.0;
        b2 = -1.0/8.0;
        b3 =  0.0;
      } else {
        a1 =  23./12.;
        a2 = -16./12.;
        a3 =   5./12.;
        b1 =  17./24.;
        b2 =  -7./24.;
        b3 =   2./24.;
      }

      for (lev=0;lev<mesh->MRABNlevels;lev++) 
        if ((Ntick+1) % (1<<lev) !=0) break; //find the max lev to update
      
      #if WADG
        for (iint l=0; l<lev; l++) {
          acousticsMRABpmlUpdate3D_wadg(mesh, a1, a2, a3, l, mesh->dt*pow(2,l));
          acousticsMRABUpdate3D_wadg(mesh, a1, a2, a3, l, t, mesh->dt*pow(2,l));
        }
        if (lev<mesh->MRABNlevels) {
          acousticsMRABpmlUpdateTrace3D_wadg(mesh, b1, b2, b3, lev, mesh->dt*pow(2,lev-1));
          acousticsMRABUpdateTrace3D_wadg(mesh, b1, b2, b3, lev, t, mesh->dt*pow(2,lev-1));
        }
      #else     
        for (iint l=0; l<lev; l++) {
          acousticsMRABpmlUpdate3D(mesh, a1, a2, a3, l, mesh->dt*pow(2,l));
          acousticsMRABUpdate3D(mesh, a1, a2, a3, l, t, mesh->dt*pow(2,l));
        }
        if (lev<mesh->MRABNlevels) {
          acousticsMRABpmlUpdateTrace3D(mesh, b1, b2, b3, lev, mesh->dt*pow(2,lev-1));
          acousticsMRABUpdateTrace3D(mesh, b1, b2, b3, lev, t, mesh->dt*pow(2,lev-1));
        }
      #endif
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
      acousticsError3D(mesh, (mesh->dt)*(tstep+1)*pow(2,mesh->MRABNlevels-1));

      // output field files
      iint fld = 3;
      char fileName[BUFSIZ];

      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      sprintf(fileName, "foo_%04d_%04d.vtu", rank, Nframe++);
      meshPlotVTU3DP(mesh, fileName, fld);
      
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
  iint haloBytes = mesh->totalHaloPairs*mesh->NfpMax*mesh->Nfaces*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  int Nframe=0;
  
  //populate the trace buffer fQ
  dfloat zero = 0.0;
  for (iint l=0; l<mesh->MRABNlevels; l++) {
    for (iint p=1;p<=mesh->NMax;p++) {
      #if WADG
      if (mesh->MRABNelP[l][p]) 
        mesh->updateKernel[p](mesh->MRABNelP[l][p],
                              mesh->o_MRABelIdsP[l][p],
                              mesh->o_N,
                              zero,
                              zero,zero,zero,
                              mesh->o_cubInterpT[p],
                              mesh->o_cubProjectT[p],
                              mesh->o_c2,
                              zero,
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              mesh->o_invVB2DT[p],
                              mesh->o_EToB,
                              mesh->o_EToE,
                              mesh->o_BBLower[p],
                              mesh->o_BBRaiseids[p],
                              mesh->o_BBRaiseVals[p],
                              mesh->o_vmapM,
                              mesh->o_rhsq,
                              mesh->o_q,
                              mesh->o_fQM,
                              mesh->o_fQP,
                              mesh->MRABshiftIndex[l]);
      if (mesh->MRABpmlNelP[l][p])
        mesh->pmlUpdateKernel[p](mesh->MRABpmlNelP[l][p],
                              mesh->o_MRABpmlElIdsP[l][p],
                              mesh->o_MRABpmlIdsP[l][p],
                              mesh->o_N,
                              zero,
                              zero,zero,zero,
                              mesh->o_cubInterpT[p],
                              mesh->o_cubProjectT[p],
                              mesh->o_c2,
                              mesh->o_EToE,
                              mesh->o_BBLower[p],
                              mesh->o_BBRaiseids[p],
                              mesh->o_BBRaiseVals[p],
                              mesh->o_vmapM,
                              mesh->o_rhsq,
                              mesh->o_pmlrhsq,
                              mesh->o_q,
                              mesh->o_pmlq,
                              mesh->o_fQM,
                              mesh->o_fQP,
                              mesh->MRABshiftIndex[l]);
      #else     
      if (mesh->MRABNelP[l][p]) 
        mesh->updateKernel[p](mesh->MRABNelP[l][p],
                              mesh->o_MRABelIdsP[l][p],
                              mesh->o_N,
                              zero,
                              zero,zero,zero,
                              zero,
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              mesh->o_invVB2DT[p],
                              mesh->o_EToB,
                              mesh->o_EToE,
                              mesh->o_BBLower[p],
                              mesh->o_BBRaiseids[p],
                              mesh->o_BBRaiseVals[p],
                              mesh->o_vmapM,
                              mesh->o_rhsq,
                              mesh->o_q,
                              mesh->o_fQM,
                              mesh->o_fQP,
                              mesh->MRABshiftIndex[l]);    
      if (mesh->MRABpmlNelP[l][p])
        mesh->pmlUpdateKernel[p](mesh->MRABpmlNelP[l][p],
                              mesh->o_MRABpmlElIdsP[l][p],
                              mesh->o_MRABpmlIdsP[l][p],
                              mesh->o_N,
                              zero,
                              zero,zero,zero,
                              mesh->o_EToE,
                              mesh->o_BBLower[p],
                              mesh->o_BBRaiseids[p],
                              mesh->o_BBRaiseVals[p],
                              mesh->o_vmapM,
                              mesh->o_rhsq,
                              mesh->o_pmlrhsq,
                              mesh->o_q,
                              mesh->o_pmlq,
                              mesh->o_fQM,
                              mesh->o_fQP,
                              mesh->MRABshiftIndex[l]);
      #endif
    }
  }

  for(iint tstep=0;tstep<mesh->NtimeSteps;++tstep){ 
    for (iint Ntick=0; Ntick < pow(2,mesh->MRABNlevels-1);Ntick++) {

      // intermediate stage time
      dfloat t = mesh->dt*(tstep*pow(2,mesh->MRABNlevels-1) + Ntick);

      iint lev;
      for (lev=0;lev<mesh->MRABNlevels;lev++) 
        if (Ntick % (1<<lev) != 0) break; //find the max lev to compute rhs

      if(mesh->totalHaloPairs>0){
        // extract halo on DEVICE
        iint Nentries = mesh->NfpMax*mesh->Nfields*mesh->Nfaces;

        mesh->haloExtractKernel(mesh->totalHaloPairs,
                    Nentries,
                    mesh->o_haloElementList,
                    mesh->o_fQP,
                    mesh->o_haloBuffer);

        // copy extracted halo to HOST 
        mesh->o_haloBuffer.copyTo(sendBuffer);      

        // start halo exchange
        meshHaloExchangeStart(mesh,
              mesh->Nfields*mesh->NfpMax*mesh->Nfaces*sizeof(dfloat),
              sendBuffer,
              recvBuffer);
      }     

      // compute volume contribution to DG acoustics RHS
      for (iint l=0;l<lev;l++) {
        for (iint p=1;p<=mesh->NMax;p++) {
          if (mesh->MRABNelP[l][p]) 
            mesh->volumeKernel[p](mesh->MRABNelP[l][p],
                                  mesh->o_MRABelIdsP[l][p],
                                  mesh->o_vgeo,
                                  mesh->o_D0ids[p],
                                  mesh->o_D1ids[p],
                                  mesh->o_D2ids[p],
                                  mesh->o_D3ids[p],
                                  mesh->o_Dvals[p],
                                  mesh->o_q,
                                  mesh->o_rhsq,
                                  mesh->MRABshiftIndex[l]);
          if (mesh->MRABpmlNelP[l][p])
            mesh->pmlVolumeKernel[p](mesh->MRABpmlNelP[l][p],
                                    mesh->o_MRABpmlElIdsP[l][p],
                                    mesh->o_MRABpmlIdsP[l][p],
                                    mesh->o_vgeo,
                                    mesh->o_D0ids[p],
                                    mesh->o_D1ids[p],
                                    mesh->o_D2ids[p],
                                    mesh->o_D3ids[p],
                                    mesh->o_Dvals[p],
                                    mesh->o_pmlSigmaX,
                                    mesh->o_pmlSigmaY,
                                    mesh->o_pmlSigmaZ,
                                    mesh->o_cubInterpT[p],
                                    mesh->o_cubProjectT[p],
                                    mesh->o_q,
                                    mesh->o_pmlq,
                                    mesh->o_rhsq,
                                    mesh->o_pmlrhsq,
                                    mesh->MRABshiftIndex[l]);
        }
      }

      if(mesh->totalHaloPairs>0){
        // wait for halo data to arrive
        meshHaloExchangeFinish(mesh);

        // copy halo data to DEVICE
        size_t offset = mesh->Nfaces*mesh->NfpMax*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
        mesh->o_fQP.copyFrom(recvBuffer, haloBytes, offset);
      }
      
      // compute surface contribution to DG acoustics RHS
      for (iint l=0;l<lev;l++) {
        for (iint p=1;p<=mesh->NMax;p++) {
          if (mesh->MRABNelP[l][p]) 
            mesh->surfaceKernel[p](mesh->MRABNelP[l][p],
                                  mesh->o_MRABelIdsP[l][p],
                                  mesh->o_sgeo,
                                  mesh->o_L0ids[p],
                                  mesh->o_L0vals[p],
                                  mesh->o_ELids[p],
                                  mesh->o_ELvals[p],
                                  mesh->o_vmapM,
                                  mesh->o_mapP,
                                  mesh->o_EToB,
                                  t,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mesh->o_q,
                                  mesh->o_fQM,
                                  mesh->o_fQP,
                                  mesh->o_rhsq,
                                  mesh->MRABshiftIndex[l]);
          if (mesh->MRABpmlNelP[l][p])
            mesh->pmlSurfaceKernel[p](mesh->MRABpmlNelP[l][p],
                                mesh->o_MRABpmlElIdsP[l][p],
                                mesh->o_MRABpmlIdsP[l][p],
                                mesh->o_sgeo,
                                mesh->o_L0ids[p],
                                mesh->o_L0vals[p],
                                mesh->o_ELids[p],
                                mesh->o_ELvals[p],
                                mesh->o_vmapM,
                                mesh->o_mapP,
                                mesh->o_EToB,
                                t,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                mesh->o_q,
                                mesh->o_pmlq,
                                mesh->o_fQM,
                                mesh->o_fQP,
                                mesh->o_rhsq,
                                mesh->o_pmlrhsq,
                                mesh->MRABshiftIndex[l]);
        }
      }

      dfloat a1, a2, a3;
      dfloat b1, b2, b3;
      if (tstep==0) {
        a1 = 1.0;
        a2 = 0.0;
        a3 = 0.0;
        b1 = 0.5;
        b2 = 0.0;
        b3 = 0.0; 
      } else if (tstep ==1) {
        a1 =  3.0/2.0;
        a2 = -1.0/2.0;
        a3 =  0.0;
        b1 =  5.0/8.0;
        b2 = -1.0/8.0;
        b3 =  0.0;
      } else {
        a1 =  23./12.;
        a2 = -16./12.;
        a3 =   5./12.;
        b1 =  17./24.;
        b2 =  -7./24.;
        b3 =   2./24.;
      }

      for (lev=0;lev<mesh->MRABNlevels;lev++) 
        if ((Ntick+1) % (1<<lev) !=0) break; //find the max lev to update
      
      #if WADG
      for (iint l=0; l<lev; l++) {
        for (iint p=1;p<=mesh->NMax;p++) {
          if (mesh->MRABNelP[l][p]) 
            mesh->updateKernel[p](mesh->MRABNelP[l][p],
                                  mesh->o_MRABelIdsP[l][p],
                                  mesh->o_N,
                                  mesh->dt*pow(2,l),
                                  a1,a2,a3,
                                  mesh->o_cubInterpT[p],
                                  mesh->o_cubProjectT[p],
                                  mesh->o_c2,
                                  t,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mesh->o_invVB2DT[p],
                                  mesh->o_EToB,
                                  mesh->o_EToE,
                                  mesh->o_BBLower[p],
                                  mesh->o_BBRaiseids[p],
                                  mesh->o_BBRaiseVals[p],
                                  mesh->o_vmapM,
                                  mesh->o_rhsq,
                                  mesh->o_q,
                                  mesh->o_fQM,
                                  mesh->o_fQP,
                                  mesh->MRABshiftIndex[l]);
          if (mesh->MRABpmlNelP[l][p])
            mesh->pmlUpdateKernel[p](mesh->MRABpmlNelP[l][p],
                                  mesh->o_MRABpmlElIdsP[l][p],
                                  mesh->o_MRABpmlIdsP[l][p],
                                  mesh->o_N,
                                  mesh->dt*pow(2,l),
                                  a1,a2,a3,
                                  mesh->o_cubInterpT[p],
                                  mesh->o_cubProjectT[p],
                                  mesh->o_c2,
                                  mesh->o_EToE,
                                  mesh->o_BBLower[p],
                                  mesh->o_BBRaiseids[p],
                                  mesh->o_BBRaiseVals[p],
                                  mesh->o_vmapM,
                                  mesh->o_rhsq,
                                  mesh->o_pmlrhsq,
                                  mesh->o_q,
                                  mesh->o_pmlq,
                                  mesh->o_fQM,
                                  mesh->o_fQP,
                                  mesh->MRABshiftIndex[l]);
        }
        //rotate index
        mesh->MRABshiftIndex[l] = (mesh->MRABshiftIndex[l]+1)%3;
      }
      if (lev<mesh->MRABNlevels) {
        for (iint p=1;p<=mesh->NMax;p++) {
          if (mesh->MRABNhaloEleP[lev][p])
            mesh->traceUpdateKernel[p](mesh->MRABNhaloEleP[lev][p],
                                      mesh->o_MRABhaloIdsP[lev][p],
                                      mesh->o_N,
                                      mesh->dt*pow(2,lev-1),
                                      b1,b2,b3,
                                      mesh->o_cubInterpT[p],
                                      mesh->o_cubProjectT[p],
                                      mesh->o_c2,
                                      t,
                                      mesh->o_x,
                                      mesh->o_y,
                                      mesh->o_z,
                                      mesh->o_invVB2DT[p],
                                      mesh->o_EToB,
                                      mesh->o_EToE,
                                      mesh->o_BBLower[p],
                                      mesh->o_BBRaiseids[p],
                                      mesh->o_BBRaiseVals[p],
                                      mesh->o_vmapM,
                                      mesh->o_rhsq,
                                      mesh->o_q,
                                      mesh->o_fQM,
                                      mesh->o_fQP,
                                      mesh->MRABshiftIndex[lev]);
          if (mesh->MRABpmlNhaloEleP[lev][p])
            mesh->pmlTraceUpdateKernel[p](mesh->MRABpmlNhaloEleP[lev][p],
                                      mesh->o_MRABpmlHaloEleIdsP[lev][p],
                                      mesh->o_MRABpmlHaloIdsP[lev][p],
                                      mesh->o_N,
                                      mesh->dt*pow(2,lev-1),
                                      b1,b2,b3,
                                      mesh->o_cubInterpT[p],
                                      mesh->o_cubProjectT[p],
                                      mesh->o_c2,
                                      mesh->o_EToE,
                                      mesh->o_BBLower[p],
                                      mesh->o_BBRaiseids[p],
                                      mesh->o_BBRaiseVals[p],
                                      mesh->o_vmapM,
                                      mesh->o_rhsq,
                                      mesh->o_pmlrhsq,
                                      mesh->o_q,
                                      mesh->o_pmlq,
                                      mesh->o_fQM,
                                      mesh->o_fQP,
                                      mesh->MRABshiftIndex[lev]);
        }
      }
      #else     
      for (iint l=0; l<lev; l++) {
        for (iint p=1;p<=mesh->NMax;p++) {
          if (mesh->MRABNelP[l][p]) 
            mesh->updateKernel[p](mesh->MRABNelP[l][p],
                                  mesh->o_MRABelIdsP[l][p],
                                  mesh->o_N,
                                  mesh->dt*pow(2,l),
                                  a1,a2,a3,
                                  t,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mesh->o_invVB2DT[p],
                                  mesh->o_EToB,
                                  mesh->o_EToE,
                                  mesh->o_BBLower[p],
                                  mesh->o_BBRaiseids[p],
                                  mesh->o_BBRaiseVals[p],
                                  mesh->o_vmapM,
                                  mesh->o_rhsq,
                                  mesh->o_q,
                                  mesh->o_fQM,
                                  mesh->o_fQP,
                                  mesh->MRABshiftIndex[l]);    
          if (mesh->MRABpmlNelP[l][p])
            mesh->pmlUpdateKernel[p](mesh->MRABpmlNelP[l][p],
                                  mesh->o_MRABpmlElIdsP[l][p],
                                  mesh->o_MRABpmlIdsP[l][p],
                                  mesh->o_N,
                                  mesh->dt*pow(2,l),
                                  a1,a2,a3,
                                  mesh->o_EToE,
                                  mesh->o_BBLower[p],
                                  mesh->o_BBRaiseids[p],
                                  mesh->o_BBRaiseVals[p],
                                  mesh->o_vmapM,
                                  mesh->o_rhsq,
                                  mesh->o_pmlrhsq,
                                  mesh->o_q,
                                  mesh->o_pmlq,
                                  mesh->o_fQM,
                                  mesh->o_fQP,
                                  mesh->MRABshiftIndex[l]);
        }
        //rotate index
        mesh->MRABshiftIndex[l] = (mesh->MRABshiftIndex[l]+1)%3;
      }

      if (lev<mesh->MRABNlevels) {
        for (iint p=1;p<=mesh->NMax;p++) {
          if (mesh->MRABNhaloEleP[lev][p]) 
            mesh->traceUpdateKernel[p](mesh->MRABNhaloEleP[lev][p],
                                      mesh->o_MRABhaloIdsP[lev][p],
                                      mesh->o_N,
                                      mesh->dt*pow(2,lev-1),
                                      b1,b2,b3,
                                      t,
                                      mesh->o_x,
                                      mesh->o_y,
                                      mesh->o_z,
                                      mesh->o_invVB2DT[p],
                                      mesh->o_EToB,
                                      mesh->o_EToE,
                                      mesh->o_BBLower[p],
                                      mesh->o_BBRaiseids[p],
                                      mesh->o_BBRaiseVals[p],
                                      mesh->o_vmapM,
                                      mesh->o_rhsq,
                                      mesh->o_q,
                                      mesh->o_fQM,
                                      mesh->o_fQP,
                                      mesh->MRABshiftIndex[lev]);
          if (mesh->MRABpmlNhaloEleP[lev][p])
            mesh->pmlTraceUpdateKernel[p](mesh->MRABpmlNhaloEleP[lev][p],
                                      mesh->o_MRABpmlHaloEleIdsP[lev][p],
                                      mesh->o_MRABpmlHaloIdsP[lev][p],
                                      mesh->o_N,
                                      mesh->dt*pow(2,lev-1),
                                      b1,b2,b3,
                                      mesh->o_EToE,
                                      mesh->o_BBLower[p],
                                      mesh->o_BBRaiseids[p],
                                      mesh->o_BBRaiseVals[p],
                                      mesh->o_vmapM,
                                      mesh->o_rhsq,
                                      mesh->o_pmlrhsq,
                                      mesh->o_q,
                                      mesh->o_pmlq,
                                      mesh->o_fQM,
                                      mesh->o_fQP,
                                      mesh->MRABshiftIndex[lev]);
        }
      }
      #endif
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

      addSourceField(mesh, mesh->q,mesh->dt*(tstep+1)*pow(2,mesh->MRABNlevels-1));

      // do error stuff on host
      acousticsError3D(mesh, (mesh->dt)*(tstep+1)*pow(2,mesh->MRABNlevels-1));

      // output field files
      iint fld = 3;
      char fileName[BUFSIZ];

      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      sprintf(fileName, "foo_%04d_%04d.vtu", rank, Nframe++);
      meshPlotVTU3DP(mesh, fileName, fld);
    }
  }

  free(recvBuffer);
  free(sendBuffer);
}


//Ricker pulse
void acousticsRickerPulse3D(dfloat x, dfloat y, dfloat z, dfloat t, dfloat f, dfloat c,
                           dfloat *u, dfloat *v, dfloat *w, dfloat *p);

void addSourceField(mesh3D *mesh, dfloat *q, dfloat t) {

  for (iint m=0;m<mesh->sourceNelements;m++) {
    iint e = mesh->sourceElements[m];
    
    iint vid = e*mesh->Nverts;
    dfloat xe1 = mesh->EX[vid+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[vid+1];
    dfloat xe3 = mesh->EX[vid+2];
    dfloat xe4 = mesh->EX[vid+3];
    
    dfloat ye1 = mesh->EY[vid+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[vid+1];
    dfloat ye3 = mesh->EY[vid+2];
    dfloat ye4 = mesh->EY[vid+3];

    dfloat ze1 = mesh->EZ[vid+0]; /* y-coordinates of vertices */
    dfloat ze2 = mesh->EZ[vid+1];
    dfloat ze3 = mesh->EZ[vid+2];
    dfloat ze4 = mesh->EZ[vid+3];
    
    for (iint n=0;n<mesh->NpMax;n++) {
      iint id = n + e*mesh->NpMax;

      dfloat rn = mesh->r[mesh->NMax][n];
      dfloat sn = mesh->s[mesh->NMax][n];
      dfloat tn = mesh->t[mesh->NMax][n];
  
      /* physical coordinate of interpolation node */
        dfloat x = -0.5*(rn+sn+tn+1.)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3 + 0.5*(1+tn)*xe4 ;
        dfloat y = -0.5*(rn+sn+tn+1.)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3 + 0.5*(1+tn)*ye4 ;
        dfloat z = -0.5*(rn+sn+tn+1.)*ze1 + 0.5*(1+rn)*ze2 + 0.5*(1+sn)*ze3 + 0.5*(1+tn)*ze4 ;

      dfloat x0 = mesh->sourceX0;
      dfloat y0 = mesh->sourceY0;
      dfloat z0 = mesh->sourceZ0;
      dfloat t0 = mesh->sourceT0;
      dfloat freq = mesh->sourceFreq;

      dfloat c = sqrt(mesh->sourceC2);

      dfloat u, v, w, p;
      acousticsRickerPulse3D(x-x0, y-y0, z-z0, t+t0, freq,c, &u, &v, &w, &p);

      iint qid = mesh->Nfields*id;
      q[qid+0] += u;
      q[qid+1] += v;
      q[qid+2] += w;
      q[qid+3] += p;
    }
  }
}
