#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void meshSurfaceGeometricFactorsQuad2D(mesh2D *mesh){

  /* unified storage array for geometric factors */
  mesh->Nsgeo = 6;
  mesh->sgeo = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
				mesh->Nsgeo*mesh->Nfp*mesh->Nfaces, 
				sizeof(dfloat));
  
  for(int e=0;e<mesh->Nelements+mesh->totalHaloPairs;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    int id = e*mesh->Nverts;

    dfloat *xe = mesh->EX + id;
    dfloat *ye = mesh->EY + id;
    
    for(int f=0;f<mesh->Nfaces;++f){ // for each face
      
      for(int i=0;i<mesh->Nfp;++i){  // for each node on face

	/* volume index of face node */
	int n = mesh->faceNodes[f*mesh->Nfp+i];

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
	int base = mesh->Nsgeo*(mesh->Nfaces*mesh->Nfp*e + mesh->Nfp*f + i);

	/* store normal, surface Jacobian, and reciprocal of volume Jacobian */
	mesh->sgeo[base+NXID] = nx/d;
	mesh->sgeo[base+NYID] = ny/d;
	mesh->sgeo[base+SJID] = d/2.;
	mesh->sgeo[base+IJID] = 1./J;

	mesh->sgeo[base+WSJID] = (d/2.)*mesh->gllw[i];
      }
    }
  }

  for(int e=0;e<mesh->Nelements;++e){ /* for each non-halo element */
    for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      int baseM = e*mesh->Nfp*mesh->Nfaces + n;
      int baseP = mesh->mapP[baseM];
      if(baseP<0) baseP = baseM;
      
      // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)
      dfloat hinvM = mesh->sgeo[baseM*mesh->Nsgeo + SJID]*mesh->sgeo[baseM*mesh->Nsgeo + IJID];
      dfloat hinvP = mesh->sgeo[baseP*mesh->Nsgeo + SJID]*mesh->sgeo[baseP*mesh->Nsgeo + IJID];
      
      mesh->sgeo[baseM*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
      mesh->sgeo[baseP*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
    }
  }

  
}
