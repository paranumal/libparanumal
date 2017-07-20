#include "acoustics2D.h"


void acousticsMRABpmlUpdate2D(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){

  dfloat *un = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *unp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
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
      mesh->q[id+0] += dt*(a1*mesh->rhsq[rhsId1+0] + a2*mesh->rhsq[rhsId2+0] + a3*mesh->rhsq[rhsId3+0]);
      mesh->q[id+1] += dt*(a1*mesh->rhsq[rhsId1+1] + a2*mesh->rhsq[rhsId2+1] + a3*mesh->rhsq[rhsId3+1]);
      mesh->q[id+2] = mesh->pmlq[pid+0]+mesh->pmlq[pid+1];
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
        pn[n] = mesh->q[qidM+2];

        mesh->fQM[qid+0] = un[n];
        mesh->fQM[qid+1] = vn[n];
        mesh->fQM[qid+2] = pn[n];
      }

      // load element neighbour
      iint eP = mesh->EToE[e*mesh->Nfaces+f];
      if (eP<0) eP = e; //boundary
      iint NP = mesh->N[eP]; 

      if (NP > N) { 
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = 0.0;
          vnp[n] = 0.0;
          pnp[n] = 0.0;
          for (iint m=0;m<2;m++){ //apply raise operator sparsly
            dfloat BBRaiseVal = mesh->BBRaiseVals[N][2*n+m];
            iint BBRaiseid = mesh->BBRaiseids[N][2*n+m];
            unp[n] += BBRaiseVal*un[BBRaiseid];
            vnp[n] += BBRaiseVal*vn[BBRaiseid];
            pnp[n] += BBRaiseVal*pn[BBRaiseid];
          }
        }
      } else if (NP < N) { 
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = 0.0;
          vnp[n] = 0.0;
          pnp[n] = 0.0;
          for (iint m=0;m<mesh->Nfp[N];m++){
            iint id = n*mesh->Nfp[N] + m;
            unp[n] += mesh->BBLower[N][id]*un[m];
            vnp[n] += mesh->BBLower[N][id]*vn[m];
            pnp[n] += mesh->BBLower[N][id]*pn[m];
          }
        }
      } else { //equal order neighbor
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = un[n];
          vnp[n] = vn[n];
          pnp[n] = pn[n];
        }
      }
    
      //write new traces to fQ
      for (iint n=0;n<mesh->Nfp[NP];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qid = mesh->Nfields*id;

        mesh->fQP[qid+0] = unp[n];
        mesh->fQP[qid+1] = vnp[n];
        mesh->fQP[qid+2] = pnp[n];
      }
    }
  }
  free(un); free(vn); free(pn);
  free(unp); free(vnp); free(pnp);
}

void acousticsMRABpmlUpdateTrace2D(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){

  dfloat *un = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *unp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
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
      s_q[n*mesh->Nfields+0] = mesh->q[id+0] + dt*(a1*mesh->rhsq[rhsId1+0] + a2*mesh->rhsq[rhsId2+0] + a3*mesh->rhsq[rhsId3+0]);
      s_q[n*mesh->Nfields+1] = mesh->q[id+1] + dt*(a1*mesh->rhsq[rhsId1+1] + a2*mesh->rhsq[rhsId2+1] + a3*mesh->rhsq[rhsId3+1]);
      s_q[n*mesh->Nfields+2] = px+py;
    }

    //project traces to proper order for neighbour
    for (iint f =0;f<mesh->Nfaces;f++) {
      //load local traces
      for (iint n=0;n<mesh->Nfp[N];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->NpMax);
        iint qid = mesh->Nfields*id;

        un[n] = s_q[qidM+0];
        vn[n] = s_q[qidM+1];
        pn[n] = s_q[qidM+2];

        mesh->fQM[qid+0] = un[n];
        mesh->fQM[qid+1] = vn[n];
        mesh->fQM[qid+2] = pn[n];
      }

      // load element neighbour
      iint eP = mesh->EToE[e*mesh->Nfaces+f];
      if (eP<0) eP = e; //boundary
      iint NP = mesh->N[eP]; 

      if (NP > N) { 
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = 0.0;
          vnp[n] = 0.0;
          pnp[n] = 0.0;
          for (iint m=0;m<2;m++){ //apply raise operator sparsly
            dfloat BBRaiseVal = mesh->BBRaiseVals[N][2*n+m];
            iint BBRaiseid = mesh->BBRaiseids[N][2*n+m];
            unp[n] += BBRaiseVal*un[BBRaiseid];
            vnp[n] += BBRaiseVal*vn[BBRaiseid];
            pnp[n] += BBRaiseVal*pn[BBRaiseid];
          }
        }
      } else if (NP < N) { 
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = 0.0;
          vnp[n] = 0.0;
          pnp[n] = 0.0;
          for (iint m=0;m<mesh->Nfp[N];m++){
            iint id = n*mesh->Nfp[N] + m;
            unp[n] += mesh->BBLower[N][id]*un[m];
            vnp[n] += mesh->BBLower[N][id]*vn[m];
            pnp[n] += mesh->BBLower[N][id]*pn[m];
          }
        }
      } else { //equal order neighbor
        for (iint n=0;n<mesh->Nfp[NP];n++){
          unp[n] = un[n];
          vnp[n] = vn[n];
          pnp[n] = pn[n];
        }
      }

      //write new traces to fQ
      for (iint n=0;n<mesh->Nfp[NP];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qid = mesh->Nfields*id;

        mesh->fQP[qid+0] = unp[n];
        mesh->fQP[qid+1] = vnp[n];
        mesh->fQP[qid+2] = pnp[n];
      }
    }
  }
  free(un); free(vn); free(pn);
  free(unp); free(vnp); free(pnp);
  free(s_q);
}


void acousticsMRABpmlUpdate2D_wadg(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){

  dfloat *un = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *unp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *p = (dfloat*) calloc(mesh->NpMax,sizeof(dfloat));
  dfloat *cubp = (dfloat*) calloc(mesh->cubNpMax,sizeof(dfloat));

  for(iint et=0;et<mesh->MRABpmlNelements[lev];et++){
    iint e = mesh->MRABpmlElementIds[lev][et];
    int N = mesh->N[e];
    iint pmlId = mesh->MRABpmlIds[lev][et];

    for(iint n=0;n<mesh->Np;++n){
      iint id = mesh->Nfields*(e*mesh->Np + n);
      iint pid = mesh->pmlNfields*(pmlId*mesh->Np + n);

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      iint pmlrhsId1 = 3*pid + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->pmlNfields;
      iint pmlrhsId2 = 3*pid + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->pmlNfields;
      iint pmlrhsId3 = 3*pid + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->pmlNfields;
      
      mesh->pmlq[pid+0] += dt*(a1*mesh->pmlrhsq[pmlrhsId1+0] + a2*mesh->pmlrhsq[pmlrhsId2+0] + a3*mesh->pmlrhsq[pmlrhsId3+0]);
      mesh->pmlq[pid+1] += dt*(a1*mesh->pmlrhsq[pmlrhsId1+1] + a2*mesh->pmlrhsq[pmlrhsId2+1] + a3*mesh->pmlrhsq[pmlrhsId3+1]);
      mesh->q[id+0] += dt*(a1*mesh->rhsq[rhsId1+0] + a2*mesh->rhsq[rhsId2+0] + a3*mesh->rhsq[rhsId3+0]);
      mesh->q[id+1] += dt*(a1*mesh->rhsq[rhsId1+1] + a2*mesh->rhsq[rhsId2+1] + a3*mesh->rhsq[rhsId3+1]);
      p[n] = mesh->pmlq[pid+0]+mesh->pmlq[pid+1];
    }

    // Interpolate rhs to cubature nodes
    for(iint n=0;n<mesh->cubNp;++n){
      cubp[n] = 0.f;
      for (iint i=0;i<mesh->Np;++i){
        cubp[n] += mesh->cubInterp[n*mesh->Np + i] * p[i];
      }
      // Multiply result by wavespeed c2 at cubature node
      cubp[n] *= mesh->c2[n + e*mesh->cubNp];
    }

    // Increment solution, project result back down
    for(iint n=0;n<mesh->Np;++n){
      // Project scaled rhs down
      dfloat c2p = 0.f;
      for (iint i=0;i<mesh->cubNp;++i){
        c2p += mesh->cubProject[n*mesh->cubNp + i] * cubp[i];
      }
      iint id = mesh->Nfields*(e*mesh->Np + n);
      mesh->q[id+2] = c2p;
    }

    //write new traces to fQ
    for (iint f =0;f<mesh->Nfaces;f++) {
      for (iint n=0;n<mesh->Nfp;n++) {
        iint id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
        iint qidM = mesh->Nfields*mesh->vmapM[id];

        iint qid = mesh->Nfields*id;
        // save trace node values of q
        for (iint fld=0; fld<mesh->Nfields;fld++) {
          mesh->fQP[qid+fld] = mesh->q[qidM+fld];
          mesh->fQM[qid+fld] = mesh->q[qidM+fld];
        }
      }
    }
  }
  free(un); free(vn); free(pn);
  free(unp); free(vnp); free(pnp);
  free(p); free(cubp);
}

void acousticsMRABpmlUpdateTrace2D_wadg(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){

  dfloat *un = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *unp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *s_q = (dfloat*) calloc(mesh->NpMax*mesh->Nfields,sizeof(dfloat));
  dfloat *cubp = (dfloat*) calloc(mesh->cubNpMax,sizeof(dfloat));

  for(iint et=0;et<mesh->MRABpmlNhaloElements[lev];et++){
    iint e = mesh->MRABpmlHaloElementIds[lev][et];
    iint pmlId = mesh->MRABpmlHaloIds[lev][et];

    for(iint n=0;n<mesh->Np;++n){
      iint id = mesh->Nfields*(e*mesh->Np + n);
      iint pid = mesh->pmlNfields*(pmlId*mesh->Np + n);

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      iint pmlrhsId1 = 3*pid + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->pmlNfields;
      iint pmlrhsId2 = 3*pid + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->pmlNfields;
      iint pmlrhsId3 = 3*pid + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->pmlNfields;
      
      // Increment solutions
      dfloat px = mesh->pmlq[pid+0] + dt*(a1*mesh->pmlrhsq[pmlrhsId1+0] + a2*mesh->pmlrhsq[pmlrhsId2+0] + a3*mesh->pmlrhsq[pmlrhsId3+0]);
      dfloat py = mesh->pmlq[pid+1] + dt*(a1*mesh->pmlrhsq[pmlrhsId1+1] + a2*mesh->pmlrhsq[pmlrhsId2+1] + a3*mesh->pmlrhsq[pmlrhsId3+1]);
      s_q[n*mesh->Nfields+0] = mesh->q[id+0] + dt*(a1*mesh->rhsq[rhsId1+0] + a2*mesh->rhsq[rhsId2+0] + a3*mesh->rhsq[rhsId3+0]);
      s_q[n*mesh->Nfields+1] = mesh->q[id+1] + dt*(a1*mesh->rhsq[rhsId1+1] + a2*mesh->rhsq[rhsId2+1] + a3*mesh->rhsq[rhsId3+1]);
      s_q[n*mesh->Nfields+2] = px+py;
    }

    // Interpolate rhs to cubature nodes
    for(iint n=0;n<mesh->cubNp;++n){
      cubp[n] = 0.f;
      for (iint i=0;i<mesh->Np;++i){
        cubp[n] += mesh->cubInterp[n*mesh->Np + i] * s_q[i*mesh->Nfields+2];
      }
      // Multiply result by wavespeed c2 at cubature node
      cubp[n] *= mesh->c2[n + e*mesh->cubNp];
    }

    // Increment solution, project result back down
    for(iint n=0;n<mesh->Np;++n){
      // Project scaled rhs down
      s_q[n*mesh->Nfields+2] = 0.f;
      for (iint i=0;i<mesh->cubNp;++i){
        s_q[n*mesh->Nfields+2] += mesh->cubProject[n*mesh->cubNp + i] * cubp[i];
      }
    }

    //write new traces to fQ
    for (iint f =0;f<mesh->Nfaces;f++) {
      for (iint n=0;n<mesh->Nfp;n++) {
        iint id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
        iint qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->Np);

        iint qid = n*mesh->Nfields + f*mesh->Nfp*mesh->Nfields + e*mesh->Nfaces*mesh->Nfp*mesh->Nfields;
        // save trace node values of q
        for (iint fld=0; fld<mesh->Nfields;fld++) {
          mesh->fQP[qid+fld] = s_q[qidM+fld];
          mesh->fQM[qid+fld] = s_q[qidM+fld];
        }
      }
    }
  }
  free(un); free(vn); free(pn);
  free(unp); free(vnp); free(pnp);
  free(s_q); free(cubp);
}

