
#include "advectionQuad3D.h"

void advectionErrorNormQuad3D(solver_t *solver, dfloat t, char *fileBase, int slice){

  mesh_t *mesh = solver->mesh;
  
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

      dfloat err = qref - solver->q[e*mesh->Np*solver->Nfields + n];
      
      l2 += JW*err*err;

      if(fileBase!=NULL)
	solver->q[e*mesh->Np*solver->Nfields+n] = l2;
    }
  }
  
  printf("%7.5lg %.2e %d (t,L2 norm err,Nlevels)\n", t, sqrt(l2),mesh->MRABNlevels);

  if(fileBase!=NULL)
    advectionPlotVTUQuad3DV2(solver, fileBase, slice);  
}
