#include "acoustics2D.h"

void acousticsRikerPulse2D(mesh2D *mesh, dfloat x,dfloat y,dfloat t,dfloat *u, dfloat *v, dfloat *p);


void acousticsSourceUpdate2D(mesh2D *mesh, iint lev, dfloat t){

  dfloat *qtmp = (dfloat *) calloc(mesh->Nfields*mesh->Np,sizeof(dfloat));

  for(iint m=0;m<mesh->MRABsourceNelements[lev];m++){
    iint e = mesh->MRABsourceElementIds[lev][m];
    iint sourceId = mesh->MRABsourceIds[lev][m];

    for(iint n=0;n<mesh->Np;++n){
      iint id = mesh->Nfields*(e*mesh->Np + n);
      iint sid = mesh->Nfields*(sourceId*mesh->Np + n);

      //subrtract the old incident field
      for (iint fld=0;fld<mesh->Nfields;fld++)
        mesh->q[id+fld] -= mesh->sourceq[sid+fld];

      //compute new the incident field (nodally)
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];
      acousticsRikerPulse2D(mesh, x, y, t, qtmp+0, qtmp+1, qtmp+2);
    }

    //Transform to BB modal space
    for (iint n=0; n<mesh->Np; n++){
      iint sid = mesh->Nfields*(sourceId*mesh->Np + n);
      mesh->sourceq[sid+0] = 0.0;
      mesh->sourceq[sid+1] = 0.0;
      mesh->sourceq[sid+2] = 0.0;
    }
    for (iint n=0;n<mesh->Np;n++){
      iint sid = mesh->Nfields*(sourceId*mesh->Np + n);
      for (iint m=0; m<mesh->Np; m++){
        mesh->sourceq[sid+0] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+0];
        mesh->sourceq[sid+1] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+1];
        mesh->sourceq[sid+2] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+2];
      }

      //add new incident field
      iint id = mesh->Nfields*(e*mesh->Np + n);
      for (iint fld=0;fld<mesh->Nfields;fld++)
        mesh->q[id+fld] += mesh->sourceq[sid+fld];
    }
  }

  free(qtmp);
}

//Riker pulse
dfloat riker(dfloat t, dfloat f) {
  return (1-2*M_PI*M_PI*f*f*t*t)*exp(-M_PI*M_PI*f*f*t*t);
}

//integrated Riker pulse
dfloat intRiker(dfloat t, dfloat f) {
  return t*exp(-M_PI*M_PI*f*f*t*t);
}

void acousticsRikerPulse2D(dfloat x, dfloat y, dfloat t, dfloat f, dfloat c, 
                           dfloat *u, dfloat *v, dfloat *p) {

  

  //radial distance
  dfloat r = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));

  *p = riker(t+t0 - r/c,f)/(4*M_PI*c*c*r);
  *u = (x-x0)*(intRiker(t+t0-r/c,f)/r + riker(t+t0-r/c,f)/c)/(4*M_PI*c*c*r*r);
  *v = (y-y0)*(intRiker(t+t0-r/c,f)/r + riker(t+t0-r/c,f)/c)/(4*M_PI*c*c*r*r);
}