#include "acoustics3D.h"

void acousticsRickerPulse3D(dfloat x, dfloat y, dfloat z, dfloat t, dfloat f, dfloat c,
                           dfloat *u, dfloat *v, dfloat *w, dfloat *p);

void acousticsMRABUpdate3D(mesh3D *mesh,
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat t, dfloat dt){

  dfloat *qtmp = (dfloat *) calloc(mesh->Nfields*mesh->Nfp,sizeof(dfloat));

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
      //check if this face is an interface of the source injection region
      iint bc = mesh->EToB[e*mesh->Nfaces+f];
      if ((bc==-10)||(bc==-11)) {
        for (iint n=0;n<mesh->Nfp;n++) {
          iint id = n + f*mesh->Nfp + e*mesh->Nfaces*mesh->Nfp;
          iint idM = mesh->vmapM[id];

          //get the nodal values of the incident field along the trace
          dfloat x = mesh->x[idM];
          dfloat y = mesh->y[idM];
          dfloat z = mesh->z[idM];

          dfloat x0 = mesh->sourceX0;
          dfloat y0 = mesh->sourceY0;
          dfloat z0 = mesh->sourceZ0;
          dfloat t0 = mesh->sourceT0;
          dfloat freq = mesh->sourceFreq;

          dfloat c = sqrt(mesh->sourceC2);

          int qid = mesh->Nfields*n;
          acousticsRickerPulse3D(x-x0, y-y0, z-z0, t+t0, freq,c, qtmp+qid+0, qtmp+qid+1, qtmp+qid+2, qtmp+qid+3);
        }

        dfloat s = 0.f;
        if (bc==-10) s= 1.f;
        if (bc==-11) s=-1.f;

        //Transform incident field trace to BB modal space and add into positive trace
        for (iint n=0; n<mesh->Nfp; n++){
          dfloat sourceu = 0.0;
          dfloat sourcev = 0.0;
          dfloat sourcew = 0.0;
          dfloat sourcep = 0.0;
          for (iint m=0; m<mesh->Nfp; m++){
            sourceu += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+0];
            sourcev += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+1];
            sourcew += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+2];
            sourcep += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+3];
          }
          //adjust positive trace values with the incident field
          iint id  = e*mesh->Nfaces*mesh->Nfp + f*mesh->Nfp + n;
          iint qidM = mesh->Nfields*mesh->vmapM[id];

          iint qid = mesh->Nfields*id;
          mesh->fQP[qid+0] = mesh->q[qidM+0] + s*sourceu;
          mesh->fQP[qid+1] = mesh->q[qidM+1] + s*sourcev;
          mesh->fQP[qid+2] = mesh->q[qidM+2] + s*sourcew;
          mesh->fQP[qid+3] = mesh->q[qidM+3] + s*sourcep;
          mesh->fQM[qid+0] = mesh->q[qidM+0];
          mesh->fQM[qid+1] = mesh->q[qidM+1];
          mesh->fQM[qid+2] = mesh->q[qidM+2];
          mesh->fQM[qid+3] = mesh->q[qidM+3];
        }
      } else {

        for (iint n=0;n<mesh->Nfp;n++) {
          iint id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
          iint qidM = mesh->Nfields*mesh->vmapM[id];

          iint qid = mesh->Nfields*id;
          // save trace node values of q
          for (iint fld=0; fld<mesh->Nfields;fld++) {
            mesh->fQM[qid+fld] = mesh->q[qidM+fld];
            mesh->fQP[qid+fld] = mesh->q[qidM+fld];
          }
        }
      }
    }
  }
  free(qtmp);

  //rotate index
  mesh->MRABshiftIndex[lev] = (mesh->MRABshiftIndex[lev]+1)%3;
}

void acousticsMRABUpdateTrace3D(mesh3D *mesh,
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat t, dfloat dt){

  dfloat *qtmp = (dfloat *) calloc(mesh->Nfields*mesh->Nfp,sizeof(dfloat));


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
      //check if this face is an interface of the source injection region
      iint bc = mesh->EToB[e*mesh->Nfaces+f];
      if ((bc==-10)||(bc==-11)) {
        for (iint n=0;n<mesh->Nfp;n++) {
          iint id = n + f*mesh->Nfp + e*mesh->Nfaces*mesh->Nfp;
          iint idM = mesh->vmapM[id];

          //get the nodal values of the incident field along the trace
          dfloat x = mesh->x[idM];
          dfloat y = mesh->y[idM];
          dfloat z = mesh->z[idM];

          dfloat x0 = mesh->sourceX0;
          dfloat y0 = mesh->sourceY0;
          dfloat z0 = mesh->sourceZ0;
          dfloat t0 = mesh->sourceT0;
          dfloat freq = mesh->sourceFreq;

          dfloat c = sqrt(mesh->sourceC2);

          int qid = mesh->Nfields*n;
          acousticsRickerPulse3D(x-x0, y-y0, z-z0, t+t0, freq,c, qtmp+qid+0, qtmp+qid+1, qtmp+qid+2, qtmp+qid+3);
        }

        dfloat s = 0.f;
        if (bc==-10) s= 1.f;
        if (bc==-11) s=-1.f;

        //Transform incident field trace to BB modal space and add into positive trace
        for (iint n=0; n<mesh->Nfp; n++){
          dfloat sourceu = 0.0;
          dfloat sourcev = 0.0;
          dfloat sourcew = 0.0;
          dfloat sourcep = 0.0;
          for (iint m=0; m<mesh->Nfp; m++){
            sourceu += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+0];
            sourcev += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+1];
            sourcew += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+2];
            sourcep += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+3];
          }
          //adjust positive trace values with the incident field
          iint id  = e*mesh->Nfaces*mesh->Nfp + f*mesh->Nfp + n;
          iint qidM = mesh->Nfields*mesh->vmapM[id];

          iint qid = mesh->Nfields*id;
          mesh->fQP[qid+0] = mesh->q[qidM+0] + s*sourceu;
          mesh->fQP[qid+1] = mesh->q[qidM+1] + s*sourcev;
          mesh->fQP[qid+2] = mesh->q[qidM+2] + s*sourcew;
          mesh->fQP[qid+3] = mesh->q[qidM+3] + s*sourcep;
          mesh->fQM[qid+0] = mesh->q[qidM+0];
          mesh->fQM[qid+1] = mesh->q[qidM+1];
          mesh->fQM[qid+2] = mesh->q[qidM+2];
          mesh->fQM[qid+3] = mesh->q[qidM+3];
        }
      } else {

        for (iint n=0;n<mesh->Nfp;n++) {
          iint id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
          iint qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->Np);

          iint qid = n*mesh->Nfields + f*mesh->Nfp*mesh->Nfields + e*mesh->Nfaces*mesh->Nfp*mesh->Nfields;
          // save trace node values of q
          for (iint fld=0; fld<mesh->Nfields;fld++) {
            mesh->fQM[qid+fld] = s_q[qidM+fld];
            mesh->fQP[qid+fld] = s_q[qidM+fld];
          }
        }
      }
    }
  }
  free(qtmp);
}

void acousticsMRABUpdate3D_wadg(mesh3D *mesh,
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat t, dfloat dt){

  dfloat *qtmp = (dfloat *) calloc(mesh->Nfields*mesh->Nfp,sizeof(dfloat));


  for(iint et=0;et<mesh->MRABNelements[lev];et++){
    iint e = mesh->MRABelementIds[lev][et];

    dfloat p[mesh->cubNp];

    // Interpolate rhs to cubature nodes
    for(iint n=0;n<mesh->cubNp;++n){
      p[n] = 0.f;
      iint id = mesh->Nfields*(e*mesh->Np);
      iint rhsId = 3*id + mesh->MRABshiftIndex[lev]*mesh->Nfields;

      for (iint i=0;i<mesh->Np;++i){
        p[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->rhsq[rhsId + 3*mesh->Nfields*i + 3];
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
      mesh->rhsq[rhsId+3] = rhsp;

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      for (iint fld=0;fld<mesh->Nfields;fld++)
        mesh->q[id+fld] += dt*(a1*mesh->rhsq[rhsId1+fld] + a2*mesh->rhsq[rhsId2+fld] + a3*mesh->rhsq[rhsId3+fld]);
    }

    //write new traces to fQ
    for (iint f =0;f<mesh->Nfaces;f++) {
      //check if this face is an interface of the source injection region
      iint bc = mesh->EToB[e*mesh->Nfaces+f];
      if ((bc==-10)||(bc==-11)) {
        for (iint n=0;n<mesh->Nfp;n++) {
          iint id = n + f*mesh->Nfp + e*mesh->Nfaces*mesh->Nfp;
          iint idM = mesh->vmapM[id];

          //get the nodal values of the incident field along the trace
          dfloat x = mesh->x[idM];
          dfloat y = mesh->y[idM];
          dfloat z = mesh->z[idM];

          dfloat x0 = mesh->sourceX0;
          dfloat y0 = mesh->sourceY0;
          dfloat z0 = mesh->sourceZ0;
          dfloat t0 = mesh->sourceT0;
          dfloat freq = mesh->sourceFreq;

          dfloat c = sqrt(mesh->sourceC2);

          int qid = mesh->Nfields*n;
          acousticsRickerPulse3D(x-x0, y-y0, z-z0, t+t0, freq,c, qtmp+qid+0, qtmp+qid+1, qtmp+qid+2, qtmp+qid+3);
        }

        dfloat s = 0.f;
        if (bc==-10) s= 1.f;
        if (bc==-11) s=-1.f;

        //Transform incident field trace to BB modal space and add into positive trace
        for (iint n=0; n<mesh->Nfp; n++){
          dfloat sourceu = 0.0;
          dfloat sourcev = 0.0;
          dfloat sourcew = 0.0;
          dfloat sourcep = 0.0;
          for (iint m=0; m<mesh->Nfp; m++){
            sourceu += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+0];
            sourcev += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+1];
            sourcew += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+2];
            sourcep += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+3];
          }

          //adjust positive trace values with the incident field
          iint id  = e*mesh->Nfaces*mesh->Nfp + f*mesh->Nfp + n;
          iint qidM = mesh->Nfields*mesh->vmapM[id];

          iint qid = mesh->Nfields*id;
          mesh->fQP[qid+0] = mesh->q[qidM+0] + s*sourceu;
          mesh->fQP[qid+1] = mesh->q[qidM+1] + s*sourcev;
          mesh->fQP[qid+2] = mesh->q[qidM+2] + s*sourcew;
          mesh->fQP[qid+3] = mesh->q[qidM+3] + s*sourcep;
          mesh->fQM[qid+0] = mesh->q[qidM+0];
          mesh->fQM[qid+1] = mesh->q[qidM+1];
          mesh->fQM[qid+2] = mesh->q[qidM+2];
          mesh->fQM[qid+3] = mesh->q[qidM+3];
        }
      } else {
        for (iint n=0;n<mesh->Nfp;n++) {
          iint id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
          iint qidM = mesh->Nfields*mesh->vmapM[id];

          iint qid = n*mesh->Nfields + f*mesh->Nfp*mesh->Nfields + e*mesh->Nfaces*mesh->Nfp*mesh->Nfields;
          // save trace node values of q
          for (iint fld=0; fld<mesh->Nfields;fld++) {
            mesh->fQM[qid+fld] = mesh->q[qidM+fld];
            mesh->fQP[qid+fld] = mesh->q[qidM+fld];
          }
        }
      }
    }
  }
  free(qtmp);

  //rotate index
  mesh->MRABshiftIndex[lev] = (mesh->MRABshiftIndex[lev]+1)%3;
}

void acousticsMRABUpdateTrace3D_wadg(mesh3D *mesh,
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat t, dfloat dt){

  dfloat *qtmp = (dfloat *) calloc(mesh->Nfields*mesh->Nfp,sizeof(dfloat));


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
        p[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->rhsq[rhsId + 3*mesh->Nfields*i + 3];
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
      rhsqn[2] = mesh->rhsq[rhsId + 2];
      rhsqn[3] = rhsp;

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      for (iint fld=0;fld<mesh->Nfields;fld++)
        s_q[n*mesh->Nfields+fld] = mesh->q[id+fld] + dt*(a1*rhsqn[fld] + a2*mesh->rhsq[rhsId2+fld] + a3*mesh->rhsq[rhsId3+fld]);
    }

    //write new traces to fQ
    for (iint f =0;f<mesh->Nfaces;f++) {
      //check if this face is an interface of the source injection region
      iint bc = mesh->EToB[e*mesh->Nfaces+f];
      if ((bc==-10)||(bc==-11)) {
        for (iint n=0;n<mesh->Nfp;n++) {
          iint id = n + f*mesh->Nfp + e*mesh->Nfaces*mesh->Nfp;
          iint idM = mesh->vmapM[id];

          //get the nodal values of the incident field along the trace
          dfloat x = mesh->x[idM];
          dfloat y = mesh->y[idM];
          dfloat z = mesh->z[idM];

          dfloat x0 = mesh->sourceX0;
          dfloat y0 = mesh->sourceY0;
          dfloat z0 = mesh->sourceZ0;
          dfloat t0 = mesh->sourceT0;
          dfloat freq = mesh->sourceFreq;

          dfloat c = sqrt(mesh->sourceC2);

          int qid = mesh->Nfields*n;
          acousticsRickerPulse3D(x-x0, y-y0, z-z0, t+t0, freq,c, qtmp+qid+0, qtmp+qid+1, qtmp+qid+2, qtmp+qid+3);
        }

        dfloat s = 0.f;
        if (bc==-10) s= 1.f;
        if (bc==-11) s=-1.f;

        //Transform incident field trace to BB modal space and add into positive trace
        for (iint n=0; n<mesh->Nfp; n++){
          dfloat sourceu = 0.0;
          dfloat sourcev = 0.0;
          dfloat sourcew = 0.0;
          dfloat sourcep = 0.0;
          for (iint m=0; m<mesh->Nfp; m++){
            sourceu += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+0];
            sourcev += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+1];
            sourcew += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+2];
            sourcep += mesh->invVB2D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+3];
          }

          //adjust positive trace values with the incident field
          iint id  = e*mesh->Nfaces*mesh->Nfp + f*mesh->Nfp + n;
          iint qidM = mesh->Nfields*mesh->vmapM[id];

          iint qid = mesh->Nfields*id;
          mesh->fQP[qid+0] = mesh->q[qidM+0] + s*sourceu;
          mesh->fQP[qid+1] = mesh->q[qidM+1] + s*sourcev;
          mesh->fQP[qid+2] = mesh->q[qidM+2] + s*sourcew;
          mesh->fQP[qid+3] = mesh->q[qidM+3] + s*sourcep;
          mesh->fQM[qid+0] = mesh->q[qidM+0];
          mesh->fQM[qid+1] = mesh->q[qidM+1];
          mesh->fQM[qid+2] = mesh->q[qidM+2];
          mesh->fQM[qid+3] = mesh->q[qidM+3];
        }
      } else {
        for (iint n=0;n<mesh->Nfp;n++) {
          iint id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
          iint qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->Np);

          iint qid = n*mesh->Nfields + f*mesh->Nfp*mesh->Nfields + e*mesh->Nfaces*mesh->Nfp*mesh->Nfields;
          // save trace node values of q
          for (iint fld=0; fld<mesh->Nfields;fld++) {
            mesh->fQM[qid+fld] = s_q[qidM+fld];
            mesh->fQP[qid+fld] = s_q[qidM+fld];
          }
        }
      }
    }
  }
  free(qtmp);
}

//Ricker pulse
dfloat ricker(dfloat t, dfloat f) {
  return (1-2*M_PI*M_PI*f*f*t*t)*exp(-M_PI*M_PI*f*f*t*t);
}

//integrated Ricker pulse
dfloat intRicker(dfloat t, dfloat f) {
  return t*exp(-M_PI*M_PI*f*f*t*t);
}

void acousticsRickerPulse3D(dfloat x, dfloat y, dfloat z, dfloat t, dfloat f, dfloat c,
                           dfloat *u, dfloat *v, dfloat *w, dfloat *p) {

  //radial distance
  dfloat r = mymax(sqrt(x*x+y*y+z*z),1e-9);

  *p = ricker(t - r/c,f)/(4*M_PI*r);
  *u = x*(intRicker(t-r/c,f)/r + ricker(t-r/c,f)/c)/(4*M_PI*r*r);
  *v = y*(intRicker(t-r/c,f)/r + ricker(t-r/c,f)/c)/(4*M_PI*r*r);
  *w = z*(intRicker(t-r/c,f)/r + ricker(t-r/c,f)/c)/(4*M_PI*r*r);
}