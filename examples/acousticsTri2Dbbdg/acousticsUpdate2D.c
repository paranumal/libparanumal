#include "mesh2D.h"

void acousticsRickerPulse2D(dfloat x, dfloat y, dfloat t, dfloat f, dfloat c, 
                           dfloat *u, dfloat *v, dfloat *p);

void acousticsMRABUpdate2D(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, int lev, dfloat t,dfloat dt){

  dfloat *qtmp = (dfloat *) calloc(mesh->Nfields*mesh->Nfp,sizeof(dfloat));

  for(int et=0;et<mesh->MRABNelements[lev];et++){
    int e = mesh->MRABelementIds[lev][et];

    for(int n=0;n<mesh->Np;++n){
      int id = mesh->Nfields*(e*mesh->Np + n);

      int rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      int rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      int rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      for (int fld=0;fld<mesh->Nfields;fld++)
        mesh->q[id+fld] += dt*(a1*mesh->rhsq[rhsId1+fld] + a2*mesh->rhsq[rhsId2+fld] + a3*mesh->rhsq[rhsId3+fld]);
    }

    //write new traces to fQ
    for (int f =0;f<mesh->Nfaces;f++) {
      //check if this face is an interface of the source injection region
      int bc = mesh->EToB[e*mesh->Nfaces+f];
      if ((bc==-10)||(bc==-11)) {
        for (int n=0;n<mesh->Nfp;n++) {
          int id = n + f*mesh->Nfp + e*mesh->Nfaces*mesh->Nfp;
          int idM = mesh->vmapM[id];
          
          //get the nodal values of the incident field along the trace
          dfloat x = mesh->x[idM];
          dfloat y = mesh->y[idM];

          dfloat x0 = mesh->sourceX0;
          dfloat y0 = mesh->sourceY0;
          dfloat t0 = mesh->sourceT0;
          dfloat freq = mesh->sourceFreq;
        
          dfloat c = sqrt(mesh->sourceC2);

          int qid = mesh->Nfields*n;
          acousticsRickerPulse2D(x-x0, y-y0, t+t0, freq,c, qtmp+qid+0, qtmp+qid+1, qtmp+qid+2);
        }

        dfloat s = 0.f;
        if (bc==-10) s= 1.f;
        if (bc==-11) s=-1.f;
        
        //Transform incident field trace to BB modal space and add into positive trace
        for (int n=0; n<mesh->Nfp; n++){
          dfloat sourceu = 0.0;
          dfloat sourcev = 0.0;
          dfloat sourcep = 0.0;
          for (int m=0; m<mesh->Nfp; m++){
            sourceu += mesh->invVB1D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+0];
            sourcev += mesh->invVB1D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+1];
            sourcep += mesh->invVB1D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+2];
          }
          //adjust positive trace values with the incident field
          int id  = e*mesh->Nfaces*mesh->Nfp + f*mesh->Nfp + n;
          int qidM = mesh->Nfields*mesh->vmapM[id];

          int qid = mesh->Nfields*id;
          mesh->fQP[qid+0] = mesh->q[qidM+0] + s*sourceu;
          mesh->fQP[qid+1] = mesh->q[qidM+1] + s*sourcev;
          mesh->fQP[qid+2] = mesh->q[qidM+2] + s*sourcep;
          mesh->fQM[qid+0] = mesh->q[qidM+0];
          mesh->fQM[qid+1] = mesh->q[qidM+1];
          mesh->fQM[qid+2] = mesh->q[qidM+2];
        }
      } else {

        for (int n=0;n<mesh->Nfp;n++) {
          int id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
          int qidM = mesh->Nfields*mesh->vmapM[id];

          int qid = mesh->Nfields*id;
          // save trace node values of q
          for (int fld=0; fld<mesh->Nfields;fld++) {
            mesh->fQP[qid+fld] = mesh->q[qidM+fld];
            mesh->fQM[qid+fld] = mesh->q[qidM+fld];
          }
        }
      }
    }
  }
  free(qtmp);

  //rotate index
  mesh->MRABshiftIndex[lev] = (mesh->MRABshiftIndex[lev]+1)%3;
}

void acousticsMRABUpdateTrace2D(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, int lev, dfloat t,dfloat dt){

  dfloat *qtmp = (dfloat *) calloc(mesh->Nfields*mesh->Nfp,sizeof(dfloat));

  for(int et=0;et<mesh->MRABNhaloElements[lev];et++){
    int e = mesh->MRABhaloIds[lev][et];

    dfloat s_q[mesh->Np*mesh->Nfields];

    for(int n=0;n<mesh->Np;++n){
      int id = mesh->Nfields*(e*mesh->Np + n);

      int rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      int rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      int rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;

      // Increment solutions
      for (int fld=0; fld < mesh->Nfields; ++fld){ 
        s_q[n*mesh->Nfields+fld] = mesh->q[id+fld] + dt*(a1*mesh->rhsq[rhsId1+fld] + a2*mesh->rhsq[rhsId2+fld] + a3*mesh->rhsq[rhsId3+fld]);
      }      
    }

    //write new traces to fQ
    for (int f =0;f<mesh->Nfaces;f++) {
      //check if this face is an interface of the source injection region
      int bc = mesh->EToB[e*mesh->Nfaces+f];
      if ((bc==-10)||(bc==-11)) {
        for (int n=0;n<mesh->Nfp;n++) {
          int id = n + f*mesh->Nfp + e*mesh->Nfaces*mesh->Nfp;
          int idM = mesh->vmapM[id];
          
          //get the nodal values of the incident field along the trace
          dfloat x = mesh->x[idM];
          dfloat y = mesh->y[idM];

          dfloat x0 = mesh->sourceX0;
          dfloat y0 = mesh->sourceY0;
          dfloat t0 = mesh->sourceT0;
          dfloat freq = mesh->sourceFreq;
        
          dfloat c = sqrt(mesh->sourceC2);

          int qid = mesh->Nfields*n;
          acousticsRickerPulse2D(x-x0, y-y0, t+t0, freq,c, qtmp+qid+0, qtmp+qid+1, qtmp+qid+2);
        }

        dfloat s = 0.f;
        if (bc==-10) s= 1.f;
        if (bc==-11) s=-1.f;

        //Transform incident field trace to BB modal space and add into positive trace
        for (int n=0; n<mesh->Nfp; n++){
          dfloat sourceu = 0.0;
          dfloat sourcev = 0.0;
          dfloat sourcep = 0.0;
          for (int m=0; m<mesh->Nfp; m++){
            sourceu += mesh->invVB1D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+0];
            sourcev += mesh->invVB1D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+1];
            sourcep += mesh->invVB1D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+2];
          }

          //adjust positive trace values with the incident field
          int id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
          int qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->Np);

          int qid = n*mesh->Nfields + f*mesh->Nfp*mesh->Nfields + e*mesh->Nfaces*mesh->Nfp*mesh->Nfields;
          mesh->fQP[qid+0] = s_q[qidM+0] + s*sourceu;
          mesh->fQP[qid+1] = s_q[qidM+1] + s*sourcev;
          mesh->fQP[qid+2] = s_q[qidM+2] + s*sourcep;
          mesh->fQM[qid+0] = s_q[qidM+0];
          mesh->fQM[qid+1] = s_q[qidM+1];
          mesh->fQM[qid+2] = s_q[qidM+2];
        }
      } else {
        for (int n=0;n<mesh->Nfp;n++) {
          int id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
          int qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->Np);

          int qid = n*mesh->Nfields + f*mesh->Nfp*mesh->Nfields + e*mesh->Nfaces*mesh->Nfp*mesh->Nfields;
          // save trace node values of q
          for (int fld=0; fld<mesh->Nfields;fld++) {
            mesh->fQM[qid+fld] = s_q[qidM+fld];
            mesh->fQP[qid+fld] = s_q[qidM+fld];
          }
        }
      }
    }
  }
  free(qtmp);
}


void acousticsMRABUpdate2D_wadg(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, int lev,dfloat t, dfloat dt){
  
  dfloat *qtmp = (dfloat *) calloc(mesh->Nfields*mesh->Nfp,sizeof(dfloat));

  for(int et=0;et<mesh->MRABNelements[lev];et++){
    int e = mesh->MRABelementIds[lev][et];

    dfloat p[mesh->cubNp];

    // Interpolate rhs to cubature nodes
    for(int n=0;n<mesh->cubNp;++n){
      p[n] = 0.f;
      int id = mesh->Nfields*(e*mesh->Np);
      int rhsId = 3*id + mesh->MRABshiftIndex[lev]*mesh->Nfields;

      for (int i=0;i<mesh->Np;++i){
        p[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->rhsq[rhsId + 3*mesh->Nfields*i + 2];
      }

      // Multiply result by wavespeed c2 at cubature node
      p[n] *= mesh->c2[n + e*mesh->cubNp];
    }

    // Increment solution, project result back down
    for(int n=0;n<mesh->Np;++n){
      // Extract velocity rhs
      int id = mesh->Nfields*(e*mesh->Np + n);
      int rhsId = 3*id + mesh->MRABshiftIndex[lev]*mesh->Nfields;
      
      // Project scaled rhs down
      dfloat rhsp = 0.f;
      for (int i=0;i<mesh->cubNp;++i){
        rhsp += mesh->cubProject[n*mesh->cubNp + i] * p[i];
      }
      mesh->rhsq[rhsId+2] = rhsp;

      int rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      int rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      int rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      for (int fld=0;fld<mesh->Nfields;fld++)
        mesh->q[id+fld] += dt*(a1*mesh->rhsq[rhsId1+fld] + a2*mesh->rhsq[rhsId2+fld] + a3*mesh->rhsq[rhsId3+fld]);
    }

    //write new traces to fQ
    for (int f =0;f<mesh->Nfaces;f++) {
      //check if this face is an interface of the source injection region
      int bc = mesh->EToB[e*mesh->Nfaces+f];
      if ((bc==-10)||(bc==-11)) {
        for (int n=0;n<mesh->Nfp;n++) {
          int id = n + f*mesh->Nfp + e*mesh->Nfaces*mesh->Nfp;
          int idM = mesh->vmapM[id];
          
          //get the nodal values of the incident field along the trace
          dfloat x = mesh->x[idM];
          dfloat y = mesh->y[idM];

          dfloat x0 = mesh->sourceX0;
          dfloat y0 = mesh->sourceY0;
          dfloat t0 = mesh->sourceT0;
          dfloat freq = mesh->sourceFreq;
        
          dfloat c = sqrt(mesh->sourceC2);

          int qid = mesh->Nfields*n;
          acousticsRickerPulse2D(x-x0, y-y0, t+t0, freq,c, qtmp+qid+0, qtmp+qid+1, qtmp+qid+2);
        }

        dfloat s = 0.f;
        if (bc==-10) s= 1.f;
        if (bc==-11) s=-1.f;

        //Transform incident field trace to BB modal space and add into positive trace
        for (int n=0; n<mesh->Nfp; n++){
          dfloat sourceu = 0.0;
          dfloat sourcev = 0.0;
          dfloat sourcep = 0.0;
          for (int m=0; m<mesh->Nfp; m++){
            sourceu += mesh->invVB1D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+0];
            sourcev += mesh->invVB1D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+1];
            sourcep += mesh->invVB1D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+2];
          }

          //adjust positive trace values with the incident field
          int id  = e*mesh->Nfaces*mesh->Nfp + f*mesh->Nfp + n;
          int qidM = mesh->Nfields*mesh->vmapM[id];

          int qid = mesh->Nfields*id;
          mesh->fQP[qid+0] = mesh->q[qidM+0] + s*sourceu;
          mesh->fQP[qid+1] = mesh->q[qidM+1] + s*sourcev;
          mesh->fQP[qid+2] = mesh->q[qidM+2] + s*sourcep;
          mesh->fQM[qid+0] = mesh->q[qidM+0];
          mesh->fQM[qid+1] = mesh->q[qidM+1];
          mesh->fQM[qid+2] = mesh->q[qidM+2];
        }
      } else {

        for (int n=0;n<mesh->Nfp;n++) {
          int id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
          int qidM = mesh->Nfields*mesh->vmapM[id];

          int qid = n*mesh->Nfields + f*mesh->Nfp*mesh->Nfields + e*mesh->Nfaces*mesh->Nfp*mesh->Nfields;
          // save trace node values of q
          for (int fld=0; fld<mesh->Nfields;fld++) {
            mesh->fQP[qid+fld] = mesh->q[qidM+fld];
            mesh->fQM[qid+fld] = mesh->q[qidM+fld];
          }
        }
      }
    }
  }
  free(qtmp);

  //rotate index
  mesh->MRABshiftIndex[lev] = (mesh->MRABshiftIndex[lev]+1)%3;
}

void acousticsMRABUpdateTrace2D_wadg(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, int lev, dfloat t, dfloat dt){
  
  dfloat *qtmp = (dfloat *) calloc(mesh->Nfields*mesh->Nfp,sizeof(dfloat));

  for(int et=0;et<mesh->MRABNhaloElements[lev];et++){
    int e = mesh->MRABhaloIds[lev][et];

    dfloat p[mesh->cubNp];
    dfloat s_q[mesh->Np*mesh->Nfields];
    dfloat rhsqn[mesh->Nfields];

    // Interpolate rhs to cubature nodes
    for(int n=0;n<mesh->cubNp;++n){
      p[n] = 0.f;
      int id = mesh->Nfields*(e*mesh->Np);
      int rhsId = 3*id + mesh->MRABshiftIndex[lev]*mesh->Nfields;

      for (int i=0;i<mesh->Np;++i){
        p[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->rhsq[rhsId + 3*mesh->Nfields*i + 2];
      }

      // Multiply result by wavespeed c2 at cubature node
      p[n] *= mesh->c2[n + e*mesh->cubNp];
    }

    // Increment solution, project result back down
    for(int n=0;n<mesh->Np;++n){
      // Extract velocity rhs
      int id = mesh->Nfields*(e*mesh->Np + n);
      int rhsId = 3*id + mesh->MRABshiftIndex[lev]*mesh->Nfields;
      
      // Project scaled rhs down
      dfloat rhsp = 0.f;
      for (int i=0;i<mesh->cubNp;++i){
        rhsp += mesh->cubProject[n*mesh->cubNp + i] * p[i];
      }
      rhsqn[0] = mesh->rhsq[rhsId + 0];
      rhsqn[1] = mesh->rhsq[rhsId + 1];  
      rhsqn[2] = rhsp;

      int rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      int rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      int rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      for (int fld=0;fld<mesh->Nfields;fld++)
        s_q[n*mesh->Nfields+fld] = mesh->q[id+fld] + dt*(a1*rhsqn[fld] + a2*mesh->rhsq[rhsId2+fld] + a3*mesh->rhsq[rhsId3+fld]);
    }

    //write new traces to fQ
    for (int f =0;f<mesh->Nfaces;f++) {
      //check if this face is an interface of the source injection region
      int bc = mesh->EToB[e*mesh->Nfaces+f];
      if ((bc==-10)||(bc==-11)) {
        for (int n=0;n<mesh->Nfp;n++) {
          int id = n + f*mesh->Nfp + e*mesh->Nfaces*mesh->Nfp;
          int idM = mesh->vmapM[id];
          
          //get the nodal values of the incident field along the trace
          dfloat x = mesh->x[idM];
          dfloat y = mesh->y[idM];

          dfloat x0 = mesh->sourceX0;
          dfloat y0 = mesh->sourceY0;
          dfloat t0 = mesh->sourceT0;
          dfloat freq = mesh->sourceFreq;
        
          dfloat c = sqrt(mesh->sourceC2);

          int qid = mesh->Nfields*n;
          acousticsRickerPulse2D(x-x0, y-y0, t+t0, freq,c, qtmp+qid+0, qtmp+qid+1, qtmp+qid+2);
        }

        dfloat s = 0.f;
        if (bc==-10) s= 1.f;
        if (bc==-11) s=-1.f;

        //Transform incident field trace to BB modal space and add into positive trace
        for (int n=0; n<mesh->Nfp; n++){
          dfloat sourceu = 0.0;
          dfloat sourcev = 0.0;
          dfloat sourcep = 0.0;
          for (int m=0; m<mesh->Nfp; m++){
            sourceu += mesh->invVB1D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+0];
            sourcev += mesh->invVB1D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+1];
            sourcep += mesh->invVB1D[n*mesh->Nfp+m]*qtmp[m*mesh->Nfields+2];
          }

          //adjust positive trace values with the incident field
          int id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
          int qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->Np);

          int qid = n*mesh->Nfields + f*mesh->Nfp*mesh->Nfields + e*mesh->Nfaces*mesh->Nfp*mesh->Nfields;
          mesh->fQP[qid+0] = s_q[qidM+0] + s*sourceu;
          mesh->fQP[qid+1] = s_q[qidM+1] + s*sourcev;
          mesh->fQP[qid+2] = s_q[qidM+2] + s*sourcep;
          mesh->fQM[qid+0] = s_q[qidM+0];
          mesh->fQM[qid+1] = s_q[qidM+1];
          mesh->fQM[qid+2] = s_q[qidM+2];
        }
      } else {
        for (int n=0;n<mesh->Nfp;n++) {
          int id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
          int qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->Np);

          int qid = n*mesh->Nfields + f*mesh->Nfp*mesh->Nfields + e*mesh->Nfaces*mesh->Nfp*mesh->Nfields;
          // save trace node values of q
          for (int fld=0; fld<mesh->Nfields;fld++) {
            mesh->fQP[qid+fld] = s_q[qidM+fld];
            mesh->fQM[qid+fld] = s_q[qidM+fld];
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

void acousticsRickerPulse2D(dfloat x, dfloat y, dfloat t, dfloat f, dfloat c, 
                           dfloat *u, dfloat *v, dfloat *p) {

  //radial distance
  dfloat r = mymax(sqrt(x*x+y*y),1e-9);

  *p = ricker(t - r/c,f)/(4*M_PI*c*c*r);
  *u = x*(intRicker(t-r/c,f)/r + ricker(t-r/c,f)/c)/(4*M_PI*c*c*r*r);
  *v = y*(intRicker(t-r/c,f)/r + ricker(t-r/c,f)/c)/(4*M_PI*c*c*r*r);
}