#include "acoustics2D.h"

void acousticsSourceSetup2D(mesh2D *mesh) {

  // location of source
  dfloat x0 = 0.f; dfloat y0 = 0.2;

  iint sourceId = -1;
  mesh->sourceNelements = 0;

  //find the element which contains this point
  for (iint e=0;e<mesh->Nelements;e++) {
    iint id = e*mesh->Nverts;

    dfloat x1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat x2 = mesh->EX[id+1];
    dfloat x3 = mesh->EX[id+2];

    dfloat y1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat y2 = mesh->EY[id+1];
    dfloat y3 = mesh->EY[id+2];

    //find the local coordinates of (x0,y0) in this element's coordinate frame.
    // if 0<=r<=1 && 0<=s<=1 && 0<=r+s<=1 the point is in this element
    dfloat J = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);

    dfloat r = ( (y3-y1)*(x0-x1)-(x3-x1)*(y0-y1))/J;
    dfloat s = (-(y2-y1)*(x0-x1)+(x2-x1)*(y0-y1))/J;

    if ((r<0)||(r>1)) continue;
    if ((s<0)||(s>1)) continue;
    if (((r+s)<0)||((r+s)>1)) continue;

    //this element contains the source point
    sourceId = e;

    //find the node which is closest to the source point and use the c2 from that node
    int minId = 0
    dfloat dist = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));    
    for(iint n=0;n<mesh->cubNp;++n){
      // cubature node coordinates
      dfloat rn = mesh->cubr[n];
      dfloat sn = mesh->cubs[n];

      /* physical coordinate of interpolation node */
      dfloat x = -0.5*(rn+sn)*x1 + 0.5*(1+rn)*x2 + 0.5*(1+sn)*x3;
      dfloat y = -0.5*(rn+sn)*y1 + 0.5*(1+rn)*y2 + 0.5*(1+sn)*y3;

      dfloat dist2 = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));

      if (dist2 < dist) {
        dist = dist2;
        minId = n;
      }
    }

    sourceC2 = mesh->c2[n+ e*mesh->cubNp];

    break;
  }

  //take the patch of elements sharing the vertices of the source element as our scatter field patch
  iint sourceV1, sourceV2, sourceV3;
  mesh->MRABsourceNelements = (int *) calloc(mesh->MRABNlevels,sizeof(int));
  mesh->MRABsourceElementIds = (iint **) calloc(mesh->MRABNlevels,sizeof(iint*));
  mesh->MRABsourceIds = (iint **) calloc(mesh->MRABNlevels,sizeof(iint*));
  if (sourceId > -1) {
    sourceV1 = mesh->EToV[sourceId*mesh->Nverts+0];
    sourceV2 = mesh->EToV[sourceId*mesh->Nverts+1];
    sourceV3 = mesh->EToV[sourceId*mesh->Nverts+2];

    for (iint e=0;e<mesh->Nelements;e++) {
      for (int n=0;n<mesh->Nverts) {
        V = mesh->EToV[e*mesh->Nverts+n];
        if ((V==sourceV1)||(V==sourceV2)||(V==sourceV3))
          mesh->MRABsourceNelements[mesh->MRABlevel[e]]++;
      }
    }

    iint cnt =0;
    for (iint lev=0;lev<mesh->MRABNlevels;lev++) {
      if (mesh->sourceNelements[lev]) {
        mesh->MRABsourceIds[lev] = (iint *) calloc(mesh->MRABsourceNelements[lev],sizeof(init));
        mesh->MRABsourceElementIds[lev] = (iint *) calloc(mesh->MRABsourceNelements[lev],sizeof(init));
        mesh->MRABsourceNelements[lev]=0;
        for (iint e=0;e<mesh->Nelements;e++) {
          for (int n=0;n<mesh->Nverts) {
            V = mesh->EToV[e*mesh->Nverts+n];
            if ((V==sourceV1)||(V==sourceV2)||(V==sourceV3)) {
              mesh->MRABsourceIds[lev][mesh->MRABsourceNelements[lev]] = cnt++;
              mesh->MRABsourceElementIds[lev][mesh->MRABsourceNelements[lev]++] = e;
            }
          }
        }
      }
    }
  }
}