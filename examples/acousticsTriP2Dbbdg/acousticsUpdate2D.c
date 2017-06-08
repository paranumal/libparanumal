#include "acoustics2D.h"

void acousticsMRABUpdate2D(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){

  dfloat *un = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *unp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  for(iint et=0;et<mesh->MRABNelements[lev];et++){
    iint e = mesh->MRABelementIds[lev][et];
    iint N = mesh->N[e];

    for(iint n=0;n<mesh->Np[N];++n){
      iint id = mesh->Nfields*(e*mesh->NpMax + n);

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      for (iint fld=0;fld<mesh->Nfields;fld++)
        mesh->q[id+fld] += dt*(a1*mesh->rhsq[rhsId1+fld] + a2*mesh->rhsq[rhsId2+fld] + a3*mesh->rhsq[rhsId3+fld]);
      
    }

    //project traces to proper order for neighbour
    for (iint f =0;f<mesh->Nfaces;f++) {
      //load local traces
      for (iint n=0;n<mesh->Nfp[N];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qidM = mesh->Nfields*mesh->vmapM[id];

        un[n] = mesh->q[qidM+0];
        vn[n] = mesh->q[qidM+1];
        pn[n] = mesh->q[qidM+2];
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
      for (iint n=0;n<mesh->Nfp[N];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qid = mesh->Nfields*id;

        mesh->fQM[qid+0] = un[n];
        mesh->fQM[qid+1] = vn[n];
        mesh->fQM[qid+2] = pn[n];
      }
    }
  }

  //rotate index
  mesh->MRABshiftIndex[lev] = (mesh->MRABshiftIndex[lev]+1)%3;

  free(un); free(vn); free(pn);
  free(unp); free(vnp); free(pnp);
}

void acousticsMRABUpdateTrace2D(mesh2D *mesh,  
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

  for(iint et=0;et<mesh->MRABNelements[lev];et++){
    iint e = mesh->MRABelementIds[lev][et];
    iint N = mesh->N[e];

    for(iint n=0;n<mesh->Np[N];++n){
      iint id = mesh->Nfields*(e*mesh->NpMax + n);

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;

      // Increment solutions
      for (iint fld=0; fld < mesh->Nfields; ++fld){ 
        s_q[n*mesh->Nfields+fld] = mesh->q[id+fld] + dt*(a1*mesh->rhsq[rhsId1+fld] + a2*mesh->rhsq[rhsId2+fld] + a3*mesh->rhsq[rhsId3+fld]);
      }      
    }

    //project traces to proper order for neighbour
    for (iint f =0;f<mesh->Nfaces;f++) {
      //load local traces
      for (iint n=0;n<mesh->Nfp[N];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->NpMax);

        un[n] = s_q[qidM+0];
        vn[n] = s_q[qidM+1];
        pn[n] = s_q[qidM+2];
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
      for (iint n=0;n<mesh->Nfp[N];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qid = mesh->Nfields*id;

        mesh->fQM[qid+0] = un[n];
        mesh->fQM[qid+1] = vn[n];
        mesh->fQM[qid+2] = pn[n];
      }
    }
  }

  free(un); free(vn); free(pn);
  free(unp); free(vnp); free(pnp);
  free(s_q);
}


void acousticsMRABUpdate2D_wadg(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){
  
  dfloat *un = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pn = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *unp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *vnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));
  dfloat *pnp = (dfloat*) calloc(mesh->NfpMax,sizeof(dfloat));

  dfloat *p = (dfloat*) calloc(mesh->cubNpMax,sizeof(dfloat));

  for(iint et=0;et<mesh->MRABNelements[lev];et++){
    iint e = mesh->MRABelementIds[lev][et];
    iint N = mesh->N[e];

    // Interpolate rhs to cubature nodes
    for(iint n=0;n<mesh->cubNp[N];++n){
      p[n] = 0.f;
      iint id = mesh->Nfields*(e*mesh->NpMax);
      iint rhsId = 3*id + mesh->MRABshiftIndex[lev]*mesh->Nfields;

      for (iint i=0;i<mesh->Np[N];++i){
        p[n] += mesh->cubInterp[N][n*mesh->Np[N] + i] * mesh->rhsq[rhsId + 3*mesh->Nfields*i + 2];
      }

      // Multiply result by wavespeed c2 at cubature node
      p[n] *= mesh->c2[n + e*mesh->cubNpMax];
    }

    // Increment solution, project result back down
    for(iint n=0;n<mesh->Np[N];++n){
      // Extract velocity rhs
      iint id = mesh->Nfields*(e*mesh->NpMax + n);
      iint rhsId = 3*id + mesh->MRABshiftIndex[lev]*mesh->Nfields;
      
      // Project scaled rhs down
      dfloat rhsp = 0.f;
      for (iint i=0;i<mesh->cubNp[N];++i){
        rhsp += mesh->cubProject[N][n*mesh->cubNp[N] + i] * p[i];
      }
      mesh->rhsq[rhsId+2] = rhsp;

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      for (iint fld=0;fld<mesh->Nfields;fld++)
        mesh->q[id+fld] += dt*(a1*mesh->rhsq[rhsId1+fld] + a2*mesh->rhsq[rhsId2+fld] + a3*mesh->rhsq[rhsId3+fld]);
    }


    //project traces to proper order for neighbour
    for (iint f=0;f<mesh->Nfaces;f++) {
      //load local traces
      for (iint n=0;n<mesh->Nfp[N];n++) {
        iint id = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qidM = mesh->Nfields*mesh->vmapM[id];

        un[n] = mesh->q[qidM+0];
        vn[n] = mesh->q[qidM+1];
        pn[n] = mesh->q[qidM+2];
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
      for (iint n=0;n<mesh->Nfp[N];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qid = mesh->Nfields*id;

        mesh->fQM[qid+0] = un[n];
        mesh->fQM[qid+1] = vn[n];
        mesh->fQM[qid+2] = pn[n];
      }
    }
  }

  //rotate index
  mesh->MRABshiftIndex[lev] = (mesh->MRABshiftIndex[lev]+1)%3;

  free(un); free(vn); free(pn);
  free(unp); free(vnp); free(pnp);
  free(p);
}

void acousticsMRABUpdateTrace2D_wadg(mesh2D *mesh,  
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
  dfloat *p   = (dfloat*) calloc(mesh->cubNpMax,sizeof(dfloat));

  for(iint et=0;et<mesh->MRABNelements[lev];et++){
    iint e = mesh->MRABelementIds[lev][et];
    iint N = mesh->N[e];

    dfloat rhsqn[mesh->Nfields];

    // Interpolate rhs to cubature nodes
    for(iint n=0;n<mesh->cubNp[N];++n){
      p[n] = 0.f;
      iint id = mesh->Nfields*(e*mesh->NpMax);
      iint rhsId = 3*id + mesh->MRABshiftIndex[lev]*mesh->Nfields;

      for (iint i=0;i<mesh->Np[N];++i){
        p[n] += mesh->cubInterp[N][n*mesh->Np[N] + i] * mesh->rhsq[rhsId + 3*mesh->Nfields*i + 2];
      }

      // Multiply result by wavespeed c2 at cubature node
      p[n] *= mesh->c2[n + e*mesh->cubNpMax];
    }

    // Increment solution, project result back down
    for(iint n=0;n<mesh->Np[N];++n){
      // Extract velocity rhs
      iint id = mesh->Nfields*(e*mesh->NpMax + n);
      iint rhsId = 3*id + mesh->MRABshiftIndex[lev]*mesh->Nfields;
      
      // Project scaled rhs down
      dfloat rhsp = 0.f;
      for (iint i=0;i<mesh->cubNp[N];++i){
        rhsp += mesh->cubProject[N][n*mesh->cubNp[N] + i] * p[i];
      }
      rhsqn[0] = mesh->rhsq[rhsId + 0];
      rhsqn[1] = mesh->rhsq[rhsId + 1];  
      rhsqn[2] = rhsp;

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      for (iint fld=0;fld<mesh->Nfields;fld++)
        s_q[n*mesh->Nfields+fld] = mesh->q[id+fld] + dt*(a1*rhsqn[fld] + a2*mesh->rhsq[rhsId2+fld] + a3*mesh->rhsq[rhsId3+fld]);
    }


    //project traces to proper order for neighbour
    for (iint f =0;f<mesh->Nfaces;f++) {
      //load local traces
      for (iint n=0;n<mesh->Nfp[N];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->NpMax);

        un[n] = s_q[qidM+0];
        vn[n] = s_q[qidM+1];
        pn[n] = s_q[qidM+2];
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
      for (iint n=0;n<mesh->Nfp[N];n++) {
        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
        iint qid = mesh->Nfields*id;

        mesh->fQM[qid+0] = un[n];
        mesh->fQM[qid+1] = vn[n];
        mesh->fQM[qid+2] = pn[n];
      }
    }
  }

  free(un); free(vn); free(pn);
  free(unp); free(vnp); free(pnp);
  free(s_q); free(p);
}