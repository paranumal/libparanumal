
#include "advectionQuad3D.h"

void advectionErrorNormQuad3D(mesh_t *mesh, dfloat t){

  dfloat l2 = 0;
  for (iint e = 0; e < mesh->Nelements; ++e) {
    for (iint n = 0; n < mesh->Np; ++n) {
      dfloat x = mesh->x[e*mesh->Np + n];
      dfloat y = mesh->y[e*mesh->Np + n];
      dfloat z = mesh->z[e*mesh->Np + n];

      //rotate reference frame back to original
      dfloat xrot = x*cos(t) + y*sin(t);
      dfloat yrot = -1*x*sin(t) + y*cos(t);
      dfloat zrot = z;
      
      //current q0 is a gaussian pulse
      dfloat qref = 1 + .1*exp(-20*((xrot-1)*(xrot-1)+yrot*yrot+zrot*zrot));
      
      dfloat JW = mesh->vgeo[mesh->Nvgeo*mesh->Np*e + JWID*mesh->Np + n];
      
      l2 += JW*(qref - mesh->q[e*mesh->Np*mesh->Nfields + n])*(qref - mesh->q[e*mesh->Np*mesh->Nfields + n]);
    //else printf("success %.15lf %.15lf\n", qref, mesh->q[e*mesh->Np*mesh->Nfields + n]);
    }
  }
  printf("%7.5lg %7.5lg (t,L2 norm err)\n", t, sqrt(l2));

}
