#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void meshSurfaceGeometricFactorsQuad2D(mesh2D *mesh){

  /* unified storage array for geometric factors */
  mesh->Nsgeo = 4;
  mesh->sgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nsgeo*mesh->Nfp*mesh->Nfaces, 
				sizeof(dfloat));
  
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    iint id = e*mesh->Nverts;

    dfloat *xe = mesh->EX + id;
    dfloat *ye = mesh->EY + id;
    
    for(iint f=0;f<mesh->Nfaces;++f){ // for each face
      
      for(iint i=0;i<mesh->Nfp;++i){  // for each node on face

	/* volume index of face node */
	iint n = mesh->faceNodes[f*mesh->Nfp+i];

	/* local node coordinates */
	dfloat rn = mesh->r[n]; 
	dfloat sn = mesh->s[n];
    
	/* Jacobian matrix */
	dfloat xr = 0.25*( (1-sn)*(xe[1]-xe[0]) + (1+sn)*(xe[2]-xe[3]) );
	dfloat xs = 0.25*( (1-rn)*(xe[3]-xe[0]) + (1+rn)*(xe[2]-xe[1]) );
	dfloat yr = 0.25*( (1-sn)*(ye[1]-ye[0]) + (1+sn)*(ye[2]-ye[3]) );
	dfloat ys = 0.25*( (1-rn)*(ye[3]-ye[0]) + (1+rn)*(ye[2]-ye[1]) );
	
	/* compute geometric factors for affine coordinate transform*/
	dfloat J = xr*ys - xs*yr;
	
	/* face f normal and length */
	dfloat nx =   ye[(f+1)%mesh->Nverts]-ye[f];
	dfloat ny = -(xe[(f+1)%mesh->Nverts]-xe[f]);
	dfloat  d = norm(nx,ny);

	/* output index */
	iint base = mesh->Nsgeo*(mesh->Nfaces*mesh->Nfp*e + mesh->Nfp*f + i);

	/* store normal, surface Jacobian, and reciprocal of volume Jacobian */
	mesh->sgeo[base+NXID] = nx/d;
	mesh->sgeo[base+NYID] = ny/d;
	mesh->sgeo[base+SJID] = d/2.;
	mesh->sgeo[base+IJID] = 1./J;
      }
    }
  }
}
