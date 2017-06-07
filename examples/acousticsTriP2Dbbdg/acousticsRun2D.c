#include "acoustics2D.h"

void acousticsRun2Dbbdg(mesh2D *mesh){

  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->NfpMax*mesh->Nfields*mesh->Nfaces*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  //temp buffer to store the solution while plotting is happening.
  //Should probably do this better, but we have a Vandermond matrix to 
  // go from degree N modal space to degree NMax nodal space, but no inverse yet.
  dfloat *qtmp = (dfloat*) calloc(mesh->Nelements*mesh->NpMax*mesh->Nfields,sizeof(dfloat));

  //populate the trace buffer fQ
  for (iint l=0;l<mesh->MRABNlevels;l++)
    acousticsMRABUpdate2D(mesh, 0., 0., 0., l, 0.);


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
      for (iint l=0;l<lev;l++)
        acousticsVolume2Dbbdg(mesh,l);

      if(mesh->totalHaloPairs>0){
        // wait for halo data to arrive
        meshHaloExchangeFinish(mesh);
    
        // copy data to halo zone of q
        memcpy(mesh->fQP+mesh->NfpMax*mesh->Nfields*mesh->Nelements*mesh->Nfaces, recvBuffer, haloBytes);
      }
      
      // compute surface contribution to DG acoustics RHS
      for (iint l=0;l<lev;l++) {
        acousticsSurface2Dbbdg(mesh,l,t);
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
          acousticsMRABUpdate2D_wadg(mesh, a1, a2, a3, l, mesh->dt*pow(2,l));
        }
        if (lev<mesh->MRABNlevels) {
          acousticsMRABUpdateTrace2D_wadg(mesh, b1, b2, b3, lev, mesh->dt*pow(2,lev-1));
        }
      #else     
        for (iint l=0; l<lev; l++) {
          acousticsMRABUpdate2D(mesh, a1, a2, a3, l, mesh->dt*pow(2,l));
        }
        if (lev<mesh->MRABNlevels) {
          acousticsMRABUpdateTrace2D(mesh, b1, b2, b3, lev, mesh->dt*pow(2,lev-1));
        }
      #endif
    }
    
    // estimate maximum error
    if((tstep%mesh->errorStep)==0){

      //Save and transform to nodal basis
      for (iint e =0;e<mesh->Nelements;e++){
        iint id = e*mesh->NpMax*mesh->Nfields;
        iint N = mesh->N[e];
        
        for (iint n=0; n<mesh->Np[N]; n++){
          qtmp[id+n*mesh->Nfields+0] = mesh->q[id+n*mesh->Nfields+0];
          qtmp[id+n*mesh->Nfields+1] = mesh->q[id+n*mesh->Nfields+1];
          qtmp[id+n*mesh->Nfields+2] = mesh->q[id+n*mesh->Nfields+2];
        }
        for (iint n=0;n<mesh->NpMax;n++) {
          mesh->q[id+n*mesh->Nfields+0] = 0.0;
          mesh->q[id+n*mesh->Nfields+1] = 0.0;
          mesh->q[id+n*mesh->Nfields+2] = 0.0;
        }
        for (iint n=0;n<mesh->NpMax;n++){
          for (iint m=0; m<mesh->Np[N]; m++){
            mesh->q[id+n*mesh->Nfields + 0] += mesh->VBplot[N][n*mesh->Np[N]+m]*qtmp[id+m*mesh->Nfields+0];
            mesh->q[id+n*mesh->Nfields + 1] += mesh->VBplot[N][n*mesh->Np[N]+m]*qtmp[id+m*mesh->Nfields+1];
            mesh->q[id+n*mesh->Nfields + 2] += mesh->VBplot[N][n*mesh->Np[N]+m]*qtmp[id+m*mesh->Nfields+2];
          }
        }
      }

      // do error stuff on host
      acousticsError2D(mesh, (mesh->dt)*(tstep+1)*pow(2,mesh->MRABNlevels-1));

      // output field files
      iint fld = 2;
      meshPlotVTU2DP(mesh, "foo", fld);

      //Recover saved q
      for (iint e =0;e<mesh->Nelements;e++){
        iint id = e*mesh->NpMax*mesh->Nfields;
        iint N = mesh->N[e];

        for (iint n=0; n<mesh->Np[N]; n++){
          mesh->q[id+n*mesh->Nfields+0] = qtmp[id+n*mesh->Nfields+0];
          mesh->q[id+n*mesh->Nfields+1] = qtmp[id+n*mesh->Nfields+1];
          mesh->q[id+n*mesh->Nfields+2] = qtmp[id+n*mesh->Nfields+2];
        }
      }
    }
  }

  free(recvBuffer);
  free(sendBuffer);
  free(qtmp);
}

void acousticsOccaRun2Dbbdg(mesh2D *mesh){

  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->NfpMax*mesh->Nfields*mesh->Nfaces*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  //populate the trace buffer fQ
  for (iint l=0; l<mesh->MRABNlevels; l++) {
    for (iint p=1;p<=mesh->NMax;p++) {
      if (mesh->MRABNelP[l][p]) {
      #if WADG
        mesh->updateKernel[p](mesh->MRABNelP[l][p],
              mesh->o_MRABelIdsP[l][p],
              mesh->o_N,
              0.0,
              0.,0.,0.,
              mesh->o_cubInterpT[p],
              mesh->o_cubProjectT[p],
              mesh->o_c2,
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
      #else     
        mesh->updateKernel[p](mesh->MRABNelP[l][p],
              mesh->o_MRABelIdsP[l][p],
              mesh->o_N,
              0.0,
              0.,0.,0.,
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
      #endif
      }
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
          if (mesh->MRABNelP[l][p]) {
            mesh->volumeKernel[p](mesh->MRABNelP[l][p],
                  mesh->o_MRABelIdsP[l][p],
                  mesh->o_vgeo,
                  mesh->o_D1ids[p],
                  mesh->o_D2ids[p],
                  mesh->o_D3ids[p],
                  mesh->o_Dvals[p],
                  mesh->o_q,
                  mesh->o_rhsq,
                  mesh->MRABshiftIndex[l]);
          }
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
          if (mesh->MRABNelP[l][p]) {
            mesh->surfaceKernel[p](mesh->MRABNelP[l][p],
                      mesh->o_MRABelIdsP[l][p],
                      mesh->o_sgeo,
                      mesh->o_L0vals[p],
                      mesh->o_ELids[p],
                      mesh->o_ELvals[p],
                      mesh->o_vmapM,
                      mesh->o_mapP,
                      mesh->o_EToB,
                      t,
                      mesh->o_x,
                      mesh->o_y,
                      mesh->o_q,
                      mesh->o_fQM,
                      mesh->o_fQP,
                      mesh->o_rhsq,
                      mesh->MRABshiftIndex[l]);
          }
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
          if (mesh->MRABNelP[l][p]) {
            mesh->updateKernel[p](mesh->MRABNelP[l][p],
                  mesh->o_MRABelIdsP[l][p],
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
                  mesh->o_q,
                  mesh->o_fQM,
                  mesh->o_fQP,
                  mesh->MRABshiftIndex[l]);
          }
        }
        //rotate index
        mesh->MRABshiftIndex[l] = (mesh->MRABshiftIndex[l]+1)%3;
      }
      if (lev<mesh->MRABNlevels) {
        for (iint p=1;p<=mesh->NMax;p++) {
          if (mesh->MRABNhaloEleP[lev][p]) {
            mesh->traceUpdateKernel[p](mesh->MRABNhaloEleP[lev][p],
                  mesh->o_MRABhaloIdsP[lev][p],
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
                  mesh->o_q,
                  mesh->o_fQM,
                  mesh->o_fQP,
                  mesh->MRABshiftIndex[lev]);
          }
        }
      }
      #else     
      for (iint l=0; l<lev; l++) {
        for (iint p=1;p<=mesh->NMax;p++) {
          if (mesh->MRABNelP[l][p]) {
            mesh->updateKernel[p](mesh->MRABNelP[l][p],
                  mesh->o_MRABelIdsP[l][p],
                  mesh->o_N,
                  mesh->dt*pow(2,l),
                  a1,a2,a3,
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
          }
        }
        //rotate index
        mesh->MRABshiftIndex[l] = (mesh->MRABshiftIndex[l]+1)%3;
      }

      if (lev<mesh->MRABNlevels) {
        for (iint p=1;p<=mesh->NMax;p++) {
          if (mesh->MRABNhaloEleP[lev][p]) {
            mesh->traceUpdateKernel[p](mesh->MRABNhaloEleP[lev][p],
                  mesh->o_MRABhaloIdsP[lev][p],
                  mesh->o_N,
                  mesh->dt*pow(2,lev-1),
                  b1,b2,b3,
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
          }
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
        }
        for (iint n=0;n<mesh->NpMax;n++) {
          mesh->q[id+n*mesh->Nfields+0] = 0.0;
          mesh->q[id+n*mesh->Nfields+1] = 0.0;
          mesh->q[id+n*mesh->Nfields+2] = 0.0;
        }
        for (iint n=0;n<mesh->NpMax;n++){
          for (iint m=0; m<mesh->Np[N]; m++){
            mesh->q[id+n*mesh->Nfields + 0] += mesh->VBplot[N][n*mesh->Np[N]+m]*qtmp[m*mesh->Nfields+0];
            mesh->q[id+n*mesh->Nfields + 1] += mesh->VBplot[N][n*mesh->Np[N]+m]*qtmp[m*mesh->Nfields+1];
            mesh->q[id+n*mesh->Nfields + 2] += mesh->VBplot[N][n*mesh->Np[N]+m]*qtmp[m*mesh->Nfields+2];
          }
        }
      }

      // do error stuff on host
      acousticsError2D(mesh, mesh->dt*(tstep+1)*pow(2,mesh->MRABNlevels-1));

      // output field files
      iint fld = 2;
      meshPlotVTU2DP(mesh, "foo", fld);
    }
  }

  free(recvBuffer);
  free(sendBuffer);
}




