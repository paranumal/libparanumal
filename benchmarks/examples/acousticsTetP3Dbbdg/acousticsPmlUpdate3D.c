#include "acoustics3D.h"


void acousticsMRABpmlUpdate3D(mesh3D *mesh,
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){

  dfloat *un = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *wn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *unp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *wnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  for(iint et=0;et<mesh->MRABpmlNelements[lev];et++){
    iint e = mesh->MRABpmlElementIds[lev][et];
    int N = mesh->N[e];
    iint pmlId = mesh->MRABpmlIds[lev][et];

    for(iint n=0;n<mesh->Np[N];++n){
      iint id = mesh->Nfields*(e*mesh->NpMax + n);
      iint pid = mesh->pmlNfields*(pmlId*mesh->NpMax + n);

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      iint pmlrhsId1 = 3*pid + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->pmlNfields;
      iint pmlrhsId2 = 3*pid + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->pmlNfields;
      iint pmlrhsId3 = 3*pid + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->pmlNfields;

      mesh->pmlq[pid+0] += dt*(a1*mesh->pmlrhsq[pmlrhsId1+0] + a2*mesh->pmlrhsq[pmlrhsId2+0] + a3*mesh->pmlrhsq[pmlrhsId3+0]);
      mesh->pmlq[pid+1] += dt*(a1*mesh->pmlrhsq[pmlrhsId1+1] + a2*mesh->pmlrhsq[pmlrhsId2+1] + a3*mesh->pmlrhsq[pmlrhsId3+1]);
      mesh->pmlq[pid+2] += dt*(a1*mesh->pmlrhsq[pmlrhsId1+2] + a2*mesh->pmlrhsq[pmlrhsId2+2] + a3*mesh->pmlrhsq[pmlrhsId3+2]);
      mesh->q[id+0] += dt*(a1*mesh->rhsq[rhsId1+0] + a2*mesh->rhsq[rhsId2+0] + a3*mesh->rhsq[rhsId3+0]);
      mesh->q[id+1] += dt*(a1*mesh->rhsq[rhsId1+1] + a2*mesh->rhsq[rhsId2+1] + a3*mesh->rhsq[rhsId3+1]);
      mesh->q[id+2] += dt*(a1*mesh->rhsq[rhsId1+2] + a2*mesh->rhsq[rhsId2+2] + a3*mesh->rhsq[rhsId3+2]);
      mesh->q[id+3] = mesh->pmlq[pid+0]+mesh->pmlq[pid+1]+mesh->pmlq[pid+2];
    }

    //project traces to proper order for neighbour
    for (iint f =0;f<mesh->Nfaces;f++) {
      //load local traces
      for (iint n=0;n<mesh->Nfp[N];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qidM = mesh->Nfields*mesh->vmapM[id];
        iint qid = mesh->Nfields*id;

        un[n] = mesh->q[qidM+0];
        vn[n] = mesh->q[qidM+1];
        wn[n] = mesh->q[qidM+2];
        pn[n] = mesh->q[qidM+3];

        mesh->fQM[qid+0] = un[n];
        mesh->fQM[qid+1] = vn[n];
        mesh->fQM[qid+2] = wn[n];
        mesh->fQM[qid+3] = pn[n];
      }

      // load element neighbour
      iint eP = mesh->EToE[e*mesh->Nfaces+f];
      if (eP<0) eP = e; //boundary
      iint NP = mesh->N[eP]; 

      if (NP > N) { 
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = 0.0;
          vnp[n] = 0.0;
          wnp[n] = 0.0;
          pnp[n] = 0.0;
          for (iint m=0;m<2;m++){ //apply raise operator sparsly
            dfloat BBRaiseVal = mesh->BBRaiseVals[N][2*n+m];
            iint BBRaiseid = mesh->BBRaiseids[N][2*n+m];
            unp[n] += BBRaiseVal*un[BBRaiseid];
            vnp[n] += BBRaiseVal*vn[BBRaiseid];
            wnp[n] += BBRaiseVal*wn[BBRaiseid];
            pnp[n] += BBRaiseVal*pn[BBRaiseid];
          }
        }
      } else if (NP < N) { 
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = 0.0;
          vnp[n] = 0.0;
          wnp[n] = 0.0;
          pnp[n] = 0.0;
          for (iint m=0;m<mesh->Nfp[N];m++){
            iint id = n*mesh->Nfp[N] + m;
            unp[n] += mesh->BBLower[N][id]*un[m];
            vnp[n] += mesh->BBLower[N][id]*vn[m];
            wnp[n] += mesh->BBLower[N][id]*wn[m];
            pnp[n] += mesh->BBLower[N][id]*pn[m];
          }
        }
      } else { //equal order neighbor
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = un[n];
          vnp[n] = vn[n];
          wnp[n] = wn[n];
          pnp[n] = pn[n];
        }
      }

      //write new traces to fQ
      for (iint n=0;n<mesh->Nfp[NP];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qid = mesh->Nfields*id;
        
        mesh->fQP[qid+0] = unp[n];
        mesh->fQP[qid+1] = vnp[n];
        mesh->fQP[qid+2] = wnp[n];
        mesh->fQP[qid+3] = pnp[n];
      }
    }
  }
  free(un); free(vn); free(wn); free(pn);
  free(unp); free(vnp); free(wnp); free(pnp);
}

void acousticsMRABpmlUpdateTrace3D(mesh3D *mesh,
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){

  dfloat *un = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *wn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *unp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *wnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *s_q = (dfloat*) calloc(mesh->NpMax*mesh->Nfields,sizeof(dfloat));

  for(iint et=0;et<mesh->MRABpmlNhaloElements[lev];et++){
    iint e = mesh->MRABpmlHaloElementIds[lev][et];
    int N = mesh->N[e];
    iint pmlId = mesh->MRABpmlHaloIds[lev][et];

    for(iint n=0;n<mesh->Np[N];++n){
      iint id = mesh->Nfields*(e*mesh->NpMax + n);
      iint pid = mesh->pmlNfields*(pmlId*mesh->NpMax + n);

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      iint pmlrhsId1 = 3*pid + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->pmlNfields;
      iint pmlrhsId2 = 3*pid + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->pmlNfields;
      iint pmlrhsId3 = 3*pid + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->pmlNfields;

      // Increment solutions
      dfloat px = mesh->pmlq[pid+0] + dt*(a1*mesh->pmlrhsq[pmlrhsId1+0] + a2*mesh->pmlrhsq[pmlrhsId2+0] + a3*mesh->pmlrhsq[pmlrhsId3+0]);
      dfloat py = mesh->pmlq[pid+1] + dt*(a1*mesh->pmlrhsq[pmlrhsId1+1] + a2*mesh->pmlrhsq[pmlrhsId2+1] + a3*mesh->pmlrhsq[pmlrhsId3+1]);
      dfloat pz = mesh->pmlq[pid+2] + dt*(a1*mesh->pmlrhsq[pmlrhsId1+2] + a2*mesh->pmlrhsq[pmlrhsId2+2] + a3*mesh->pmlrhsq[pmlrhsId3+2]);
      s_q[n*mesh->Nfields+0] = mesh->q[id+0] + dt*(a1*mesh->rhsq[rhsId1+0] + a2*mesh->rhsq[rhsId2+0] + a3*mesh->rhsq[rhsId3+0]);
      s_q[n*mesh->Nfields+1] = mesh->q[id+1] + dt*(a1*mesh->rhsq[rhsId1+1] + a2*mesh->rhsq[rhsId2+1] + a3*mesh->rhsq[rhsId3+1]);
      s_q[n*mesh->Nfields+2] = mesh->q[id+2] + dt*(a1*mesh->rhsq[rhsId1+2] + a2*mesh->rhsq[rhsId2+2] + a3*mesh->rhsq[rhsId3+2]);
      s_q[n*mesh->Nfields+3] = px+py+pz;
    }

    //write new traces to fQ
    for (iint f =0;f<mesh->Nfaces;f++) {
      //load local traces
      for (iint n=0;n<mesh->Nfp[N];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->NpMax);
        iint qid = mesh->Nfields*id;

        un[n] = s_q[qidM+0];
        vn[n] = s_q[qidM+1];
        wn[n] = s_q[qidM+2];
        pn[n] = s_q[qidM+3];

        mesh->fQM[qid+0] = un[n];
        mesh->fQM[qid+1] = vn[n];
        mesh->fQM[qid+2] = wn[n];
        mesh->fQM[qid+3] = pn[n];
      }

      // load element neighbour
      iint eP = mesh->EToE[e*mesh->Nfaces+f];
      if (eP<0) eP = e; //boundary
      iint NP = mesh->N[eP]; 

      if (NP > N) { 
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = 0.0;
          vnp[n] = 0.0;
          wnp[n] = 0.0;
          pnp[n] = 0.0;
          for (iint m=0;m<2;m++){ //apply raise operator sparsly
            dfloat BBRaiseVal = mesh->BBRaiseVals[N][2*n+m];
            iint BBRaiseid = mesh->BBRaiseids[N][2*n+m];
            unp[n] += BBRaiseVal*un[BBRaiseid];
            vnp[n] += BBRaiseVal*vn[BBRaiseid];
            wnp[n] += BBRaiseVal*wn[BBRaiseid];
            pnp[n] += BBRaiseVal*pn[BBRaiseid];
          }
        }
      } else if (NP < N) { 
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = 0.0;
          vnp[n] = 0.0;
          wnp[n] = 0.0;
          pnp[n] = 0.0;
          for (iint m=0;m<mesh->Nfp[N];m++){
            iint id = n*mesh->Nfp[N] + m;
            unp[n] += mesh->BBLower[N][id]*un[m];
            vnp[n] += mesh->BBLower[N][id]*vn[m];
            wnp[n] += mesh->BBLower[N][id]*wn[m];
            pnp[n] += mesh->BBLower[N][id]*pn[m];
          }
        }
      } else { //equal order neighbor
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = un[n];
          vnp[n] = vn[n];
          wnp[n] = wn[n];
          pnp[n] = pn[n];
        }
      }

      //write new traces to fQ
      for (iint n=0;n<mesh->Nfp[NP];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qid = mesh->Nfields*id;
        
        mesh->fQP[qid+0] = unp[n];
        mesh->fQP[qid+1] = vnp[n];
        mesh->fQP[qid+2] = wnp[n];
        mesh->fQP[qid+3] = pnp[n];
      }
    }
  }
  free(un); free(vn); free(wn); free(pn);
  free(unp); free(vnp); free(wnp); free(pnp);
  free(s_q);
}


void acousticsMRABpmlUpdate3D_wadg(mesh3D *mesh,
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){

  dfloat *un = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *wn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *unp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *wnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *p = (dfloat*) calloc(mesh->NpMax,sizeof(dfloat));
  dfloat *cubp = (dfloat*) calloc(mesh->cubNpMax,sizeof(dfloat));

  for(iint et=0;et<mesh->MRABpmlNelements[lev];et++){
    iint e = mesh->MRABpmlElementIds[lev][et];
    int N = mesh->N[e];
    iint pmlId = mesh->MRABpmlIds[lev][et];

    for(iint n=0;n<mesh->Np[N];++n){
      iint id = mesh->Nfields*(e*mesh->NpMax + n);
      iint pid = mesh->pmlNfields*(pmlId*mesh->NpMax + n);

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      iint pmlrhsId1 = 3*pid + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->pmlNfields;
      iint pmlrhsId2 = 3*pid + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->pmlNfields;
      iint pmlrhsId3 = 3*pid + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->pmlNfields;

      mesh->pmlq[pid+0] += dt*(a1*mesh->pmlrhsq[pmlrhsId1+0] + a2*mesh->pmlrhsq[pmlrhsId2+0] + a3*mesh->pmlrhsq[pmlrhsId3+0]);
      mesh->pmlq[pid+1] += dt*(a1*mesh->pmlrhsq[pmlrhsId1+1] + a2*mesh->pmlrhsq[pmlrhsId2+1] + a3*mesh->pmlrhsq[pmlrhsId3+1]);
      mesh->pmlq[pid+2] += dt*(a1*mesh->pmlrhsq[pmlrhsId1+2] + a2*mesh->pmlrhsq[pmlrhsId2+2] + a3*mesh->pmlrhsq[pmlrhsId3+2]);
      mesh->q[id+0] += dt*(a1*mesh->rhsq[rhsId1+0] + a2*mesh->rhsq[rhsId2+0] + a3*mesh->rhsq[rhsId3+0]);
      mesh->q[id+1] += dt*(a1*mesh->rhsq[rhsId1+1] + a2*mesh->rhsq[rhsId2+1] + a3*mesh->rhsq[rhsId3+1]);
      mesh->q[id+2] += dt*(a1*mesh->rhsq[rhsId1+2] + a2*mesh->rhsq[rhsId2+2] + a3*mesh->rhsq[rhsId3+2]);
      p[n] = mesh->pmlq[pid+0]+mesh->pmlq[pid+1]+mesh->pmlq[pid+2];
    }

    // Interpolate rhs to cubature nodes
    for(iint n=0;n<mesh->cubNp[N];++n){
      cubp[n] = 0.f;
      for (iint i=0;i<mesh->Np[N];++i){
        cubp[n] += mesh->cubInterp[N][n*mesh->Np[N] + i] * p[i];
      }
      // Multiply result by wavespeed c2 at cubature node
      cubp[n] *= mesh->c2[n + e*mesh->cubNpMax];
    }

    // Increment solution, project result back down
    for(iint n=0;n<mesh->Np[N];++n){
      // Project scaled rhs down
      dfloat c2p = 0.f;
      for (iint i=0;i<mesh->cubNp[N];++i){
        c2p += mesh->cubProject[N][n*mesh->cubNp[N] + i] * cubp[i];
      }
      iint id = mesh->Nfields*(e*mesh->NpMax + n);
      mesh->q[id+3] = c2p;
    }

    //project traces to proper order for neighbour
    for (iint f =0;f<mesh->Nfaces;f++) {
      //load local traces
      for (iint n=0;n<mesh->Nfp[N];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qidM = mesh->Nfields*mesh->vmapM[id];
        iint qid = mesh->Nfields*id;

        un[n] = mesh->q[qidM+0];
        vn[n] = mesh->q[qidM+1];
        wn[n] = mesh->q[qidM+2];
        pn[n] = mesh->q[qidM+3];

        mesh->fQM[qid+0] = un[n];
        mesh->fQM[qid+1] = vn[n];
        mesh->fQM[qid+2] = wn[n];
        mesh->fQM[qid+3] = pn[n];
      }

      // load element neighbour
      iint eP = mesh->EToE[e*mesh->Nfaces+f];
      if (eP<0) eP = e; //boundary
      iint NP = mesh->N[eP]; 

      if (NP > N) { 
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = 0.0;
          vnp[n] = 0.0;
          wnp[n] = 0.0;
          pnp[n] = 0.0;
          for (iint m=0;m<2;m++){ //apply raise operator sparsly
            dfloat BBRaiseVal = mesh->BBRaiseVals[N][2*n+m];
            iint BBRaiseid = mesh->BBRaiseids[N][2*n+m];
            unp[n] += BBRaiseVal*un[BBRaiseid];
            vnp[n] += BBRaiseVal*vn[BBRaiseid];
            wnp[n] += BBRaiseVal*wn[BBRaiseid];
            pnp[n] += BBRaiseVal*pn[BBRaiseid];
          }
        }
      } else if (NP < N) { 
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = 0.0;
          vnp[n] = 0.0;
          wnp[n] = 0.0;
          pnp[n] = 0.0;
          for (iint m=0;m<mesh->Nfp[N];m++){
            iint id = n*mesh->Nfp[N] + m;
            unp[n] += mesh->BBLower[N][id]*un[m];
            vnp[n] += mesh->BBLower[N][id]*vn[m];
            wnp[n] += mesh->BBLower[N][id]*wn[m];
            pnp[n] += mesh->BBLower[N][id]*pn[m];
          }
        }
      } else { //equal order neighbor
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = un[n];
          vnp[n] = vn[n];
          wnp[n] = wn[n];
          pnp[n] = pn[n];
        }
      }

      //write new traces to fQ
      for (iint n=0;n<mesh->Nfp[NP];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qid = mesh->Nfields*id;
        
        mesh->fQP[qid+0] = unp[n];
        mesh->fQP[qid+1] = vnp[n];
        mesh->fQP[qid+2] = wnp[n];
        mesh->fQP[qid+3] = pnp[n];
      }
    }
  }
  free(un); free(vn); free(wn); free(pn);
  free(unp); free(vnp); free(wnp); free(pnp);
  free(p); free(cubp);
}

void acousticsMRABpmlUpdateTrace3D_wadg(mesh3D *mesh,
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){

  dfloat *un = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *wn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *unp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *wnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *s_q = (dfloat*) calloc(mesh->NpMax*mesh->Nfields,sizeof(dfloat));
  dfloat *cubp = (dfloat*) calloc(mesh->cubNpMax,sizeof(dfloat));

  for(iint et=0;et<mesh->MRABpmlNhaloElements[lev];et++){
    iint e = mesh->MRABpmlHaloElementIds[lev][et];
    int N = mesh->N[e];
    iint pmlId = mesh->MRABpmlHaloIds[lev][et];

    for(iint n=0;n<mesh->Np[N];++n){
      iint id = mesh->Nfields*(e*mesh->NpMax + n);
      iint pid = mesh->pmlNfields*(pmlId*mesh->NpMax + n);

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      iint pmlrhsId1 = 3*pid + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->pmlNfields;
      iint pmlrhsId2 = 3*pid + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->pmlNfields;
      iint pmlrhsId3 = 3*pid + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->pmlNfields;

      // Increment solutions
      dfloat px = mesh->pmlq[pid+0] + dt*(a1*mesh->pmlrhsq[pmlrhsId1+0] + a2*mesh->pmlrhsq[pmlrhsId2+0] + a3*mesh->pmlrhsq[pmlrhsId3+0]);
      dfloat py = mesh->pmlq[pid+1] + dt*(a1*mesh->pmlrhsq[pmlrhsId1+1] + a2*mesh->pmlrhsq[pmlrhsId2+1] + a3*mesh->pmlrhsq[pmlrhsId3+1]);
      dfloat pz = mesh->pmlq[pid+2] + dt*(a1*mesh->pmlrhsq[pmlrhsId1+2] + a2*mesh->pmlrhsq[pmlrhsId2+2] + a3*mesh->pmlrhsq[pmlrhsId3+2]);
      s_q[n*mesh->Nfields+0] = mesh->q[id+0] + dt*(a1*mesh->rhsq[rhsId1+0] + a2*mesh->rhsq[rhsId2+0] + a3*mesh->rhsq[rhsId3+0]);
      s_q[n*mesh->Nfields+1] = mesh->q[id+1] + dt*(a1*mesh->rhsq[rhsId1+1] + a2*mesh->rhsq[rhsId2+1] + a3*mesh->rhsq[rhsId3+1]);
      s_q[n*mesh->Nfields+2] = mesh->q[id+2] + dt*(a1*mesh->rhsq[rhsId1+2] + a2*mesh->rhsq[rhsId2+2] + a3*mesh->rhsq[rhsId3+2]);
      s_q[n*mesh->Nfields+2] = px+py+pz;
    }

    // Interpolate rhs to cubature nodes
    for(iint n=0;n<mesh->cubNp[N];++n){
      cubp[n] = 0.f;
      for (iint i=0;i<mesh->Np[N];++i){
        cubp[n] += mesh->cubInterp[N][n*mesh->Np[N] + i] * s_q[i*mesh->Nfields+2];
      }
      // Multiply result by wavespeed c2 at cubature node
      cubp[n] *= mesh->c2[n + e*mesh->cubNpMax];
    }

    // Increment solution, project result back down
    for(iint n=0;n<mesh->Np[N];++n){
      // Project scaled rhs down
      s_q[n*mesh->Nfields+3] = 0.f;
      for (iint i=0;i<mesh->cubNp[N];++i){
        s_q[n*mesh->Nfields+3] += mesh->cubProject[N][n*mesh->cubNp[N] + i] * cubp[i];
      }
    }

    //write new traces to fQ
    for (iint f =0;f<mesh->Nfaces;f++) {
      //load local traces
      for (iint n=0;n<mesh->Nfp[N];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->NpMax);
        iint qid = mesh->Nfields*id;

        un[n] = s_q[qidM+0];
        vn[n] = s_q[qidM+1];
        wn[n] = s_q[qidM+2];
        pn[n] = s_q[qidM+3];

        mesh->fQM[qid+0] = un[n];
        mesh->fQM[qid+1] = vn[n];
        mesh->fQM[qid+2] = wn[n];
        mesh->fQM[qid+3] = pn[n];
      }

      // load element neighbour
      iint eP = mesh->EToE[e*mesh->Nfaces+f];
      if (eP<0) eP = e; //boundary
      iint NP = mesh->N[eP]; 

      if (NP > N) { 
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = 0.0;
          vnp[n] = 0.0;
          wnp[n] = 0.0;
          pnp[n] = 0.0;
          for (iint m=0;m<2;m++){ //apply raise operator sparsly
            dfloat BBRaiseVal = mesh->BBRaiseVals[N][2*n+m];
            iint BBRaiseid = mesh->BBRaiseids[N][2*n+m];
            unp[n] += BBRaiseVal*un[BBRaiseid];
            vnp[n] += BBRaiseVal*vn[BBRaiseid];
            wnp[n] += BBRaiseVal*wn[BBRaiseid];
            pnp[n] += BBRaiseVal*pn[BBRaiseid];
          }
        }
      } else if (NP < N) { 
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = 0.0;
          vnp[n] = 0.0;
          wnp[n] = 0.0;
          pnp[n] = 0.0;
          for (iint m=0;m<mesh->Nfp[N];m++){
            iint id = n*mesh->Nfp[N] + m;
            unp[n] += mesh->BBLower[N][id]*un[m];
            vnp[n] += mesh->BBLower[N][id]*vn[m];
            wnp[n] += mesh->BBLower[N][id]*wn[m];
            pnp[n] += mesh->BBLower[N][id]*pn[m];
          }
        }
      } else { //equal order neighbor
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = un[n];
          vnp[n] = vn[n];
          wnp[n] = wn[n];
          pnp[n] = pn[n];
        }
      }

      //write new traces to fQ
      for (iint n=0;n<mesh->Nfp[NP];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qid = mesh->Nfields*id;
        
        mesh->fQP[qid+0] = unp[n];
        mesh->fQP[qid+1] = vnp[n];
        mesh->fQP[qid+2] = wnp[n];
        mesh->fQP[qid+3] = pnp[n];
      }
    }
  }
  free(un); free(vn); free(wn); free(pn);
  free(unp); free(vnp); free(wnp); free(pnp);
  free(s_q); free(cubp);
}



