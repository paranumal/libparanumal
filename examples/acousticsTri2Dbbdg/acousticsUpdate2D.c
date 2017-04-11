#include "mesh2D.h"

void acousticsMRABUpdate2D(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){

  for(iint et=0;et<mesh->MRABNelements[lev];et++){
    iint e = mesh->MRABelementIds[lev][et];

    for(iint n=0;n<mesh->Np;++n){
      iint id = mesh->Nfields*(e*mesh->Np + n);

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      for (iint fld=0;fld<mesh->Nfields;fld++)
        mesh->q[id+fld] += dt*(a1*mesh->rhsq[rhsId1+fld] + a2*mesh->rhsq[rhsId2+fld] + a3*mesh->rhsq[rhsId3+fld]);
    }

    //write new traces to fQ
    for (iint f =0;f<mesh->Nfaces;f++) {
      for (iint n=0;n<mesh->Nfp;n++) {
        iint id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
        iint qidM = mesh->Nfields*mesh->vmapM[id];

        iint qid = mesh->Nfields*id;
        // save trace node values of q
        for (iint fld=0; fld<mesh->Nfields;fld++) {
          mesh->fQ[qid+fld] = mesh->q[qidM+fld];
        }
      }
    }
  }

  //rotate index
  mesh->MRABshiftIndex[lev] = (mesh->MRABshiftIndex[lev]+1)%3;
}

void acousticsMRABUpdateTrace2D(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){

  for(iint et=0;et<mesh->MRABNelements[lev];et++){
    iint e = mesh->MRABelementIds[lev][et];

    dfloat s_q[mesh->Np*mesh->Nfields];

    for(iint n=0;n<mesh->Np;++n){
      iint id = mesh->Nfields*(e*mesh->Np + n);

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;

      // Increment solutions
      for (iint fld=0; fld < mesh->Nfields; ++fld){ 
        s_q[n*mesh->Nfields+fld] = mesh->q[id+fld] + dt*(a1*mesh->rhsq[rhsId1+fld] + a2*mesh->rhsq[rhsId2+fld] + a3*mesh->rhsq[rhsId3+fld]);
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
          mesh->fQ[qid+fld] = s_q[qidM+fld];
        }
      }
    }
  }
}


void acousticsMRABUpdate2D_wadg(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){
  
  for(iint et=0;et<mesh->MRABNelements[lev];et++){
    iint e = mesh->MRABelementIds[lev][et];

    dfloat p[mesh->cubNp];

    // Interpolate rhs to cubature nodes
    for(iint n=0;n<mesh->cubNp;++n){
      p[n] = 0.f;
      iint id = mesh->Nfields*(e*mesh->Np);
      iint rhsId = 3*id + mesh->MRABshiftIndex[lev]*mesh->Nfields;

      for (iint i=0;i<mesh->Np;++i){
        p[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->rhsq[rhsId + 3*mesh->Nfields*i + 2];
      }

      // Multiply result by wavespeed c2 at cubature node
      p[n] *= mesh->c2[n + e*mesh->cubNp];
    }

    // Increment solution, project result back down
    for(iint n=0;n<mesh->Np;++n){
      // Extract velocity rhs
      iint id = mesh->Nfields*(e*mesh->Np + n);
      iint rhsId = 3*id + mesh->MRABshiftIndex[lev]*mesh->Nfields;
      
      // Project scaled rhs down
      dfloat rhsp = 0.f;
      for (iint i=0;i<mesh->cubNp;++i){
        rhsp += mesh->cubProject[n*mesh->cubNp + i] * p[i];
      }
      mesh->rhsq[rhsId+2] = rhsp;

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      for (iint fld=0;fld<mesh->Nfields;fld++)
        mesh->q[id+fld] += dt*(a1*mesh->rhsq[rhsId1+fld] + a2*mesh->rhsq[rhsId2+fld] + a3*mesh->rhsq[rhsId3+fld]);
    }

    //write new traces to fQ
    for (iint f =0;f<mesh->Nfaces;f++) {
      for (iint n=0;n<mesh->Nfp;n++) {
        iint id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
        iint qidM = mesh->Nfields*mesh->vmapM[id];

        iint qid = n*mesh->Nfields + f*mesh->Nfp*mesh->Nfields + e*mesh->Nfaces*mesh->Nfp*mesh->Nfields;
        // save trace node values of q
        for (iint fld=0; fld<mesh->Nfields;fld++) {
          mesh->fQ[qid+fld] = mesh->q[qidM+fld];
        }
      }
    }
  }

  //rotate index
  mesh->MRABshiftIndex[lev] = (mesh->MRABshiftIndex[lev]+1)%3;
}

void acousticsMRABUpdateTrace2D_wadg(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){
  
  for(iint et=0;et<mesh->MRABNelements[lev];et++){
    iint e = mesh->MRABelementIds[lev][et];

    dfloat p[mesh->cubNp];
    dfloat s_q[mesh->Np*mesh->Nfields];
    dfloat rhsqn[mesh->Nfields];

    // Interpolate rhs to cubature nodes
    for(iint n=0;n<mesh->cubNp;++n){
      p[n] = 0.f;
      iint id = mesh->Nfields*(e*mesh->Np);
      iint rhsId = 3*id + mesh->MRABshiftIndex[lev]*mesh->Nfields;

      for (iint i=0;i<mesh->Np;++i){
        p[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->rhsq[rhsId + 3*mesh->Nfields*i + 2];
      }

      // Multiply result by wavespeed c2 at cubature node
      p[n] *= mesh->c2[n + e*mesh->cubNp];
    }

    // Increment solution, project result back down
    for(iint n=0;n<mesh->Np;++n){
      // Extract velocity rhs
      iint id = mesh->Nfields*(e*mesh->Np + n);
      iint rhsId = 3*id + mesh->MRABshiftIndex[lev]*mesh->Nfields;
      
      // Project scaled rhs down
      dfloat rhsp = 0.f;
      for (iint i=0;i<mesh->cubNp;++i){
        rhsp += mesh->cubProject[n*mesh->cubNp + i] * p[i];
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

    //write new traces to fQ
    for (iint f =0;f<mesh->Nfaces;f++) {
      for (iint n=0;n<mesh->Nfp;n++) {
        iint id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
        iint qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->Np);

        iint qid = n*mesh->Nfields + f*mesh->Nfp*mesh->Nfields + e*mesh->Nfaces*mesh->Nfp*mesh->Nfields;
        // save trace node values of q
        for (iint fld=0; fld<mesh->Nfields;fld++) {
          mesh->fQ[qid+fld] = s_q[qidM+fld];
        }
      }
    }
  }
}