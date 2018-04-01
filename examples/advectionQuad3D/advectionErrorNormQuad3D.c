
#include "advectionQuad3D.h"

void advectionErrorNormQuad3D(mesh_t *mesh, dfloat t, char *fileBase, int slice){

  dfloat l2 = 0;
  for (iint e = 0; e < mesh->Nelements; ++e) {
    for (iint n = 0; n < mesh->Np; ++n) {
      dfloat x = mesh->x[e*mesh->Np + n];
      dfloat y = mesh->y[e*mesh->Np + n];
      dfloat z = mesh->z[e*mesh->Np + n];

      //rotate reference frame back to original
      dfloat xrotP = x*cos(t) + y*sin(t);
      dfloat yrotP = -1*x*sin(t) + y*cos(t);
      dfloat zrotP = z;

      dfloat xrotM = x*cos(-t) + y*sin(-t);
      dfloat yrotM = -1*x*sin(-t) + y*cos(-t);
      dfloat zrotM = z;
      
      //current q0 is a gaussian pulse
      dfloat qrefM = 1 + .1*exp(-20*((xrotM-1)*(xrotM-1)+yrotM*yrotM+zrotM*zrotM));
      dfloat qrefP = 1 + .1*exp(-20*((xrotP-1)*(xrotP-1)+yrotP*yrotP+zrotP*zrotP));

      dfloat qref = 0.5*(qrefM+qrefP);
      
      dfloat JW = mesh->vgeo[mesh->Nvgeo*mesh->Np*e + JWID*mesh->Np + n];

      dfloat err = qref - mesh->q[e*mesh->Np*mesh->Nfields + n];
      
      l2 += JW*err*err;

      if(fileBase!=NULL)
	mesh->q[e*mesh->Np*mesh->Nfields+n] = err;
    }
  }

  printf("%7.5lg %.2e (t,L2 norm err)\n", t, sqrt(l2));

  if(fileBase!=NULL)
    advectionPlotVTUQuad3DV2(mesh, fileBase, slice);  
}
