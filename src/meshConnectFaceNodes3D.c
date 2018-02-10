#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mesh3D.h"

iint findBestMatch(dfloat x1, dfloat y1, dfloat z1,
		   iint Np2, iint *nodeList, dfloat *x2, dfloat *y2, dfloat *z2, int *nP,int face, int faceP){
  
  iint matchIndex;
  dfloat mindist2=1e9;

  for(int n=0;n<Np2;++n){
    
    /* next node */
    const iint i2 = nodeList[n];

    /* distance between target and next node */
    const dfloat dist2 = pow(x1-x2[i2],2) + pow(y1-y2[i2],2) + pow(z1-z2[i2],2);
    
    /* if next node is closer to target update match */
    if(n==0 || dist2<mindist2){
      mindist2 = dist2;
      matchIndex = i2;
      *nP = n;
    }
  }
  if(mindist2>1e-3) printf("arggh - bad match: x,y,z=%g,%g,%g,   %d,%d\n", x1,y1,z1,face,faceP);

  return matchIndex;
}
		

// serial face-node to face-node connection
void meshConnectFaceNodes3D(mesh3D *mesh){
  
  /* volume indices of the interior and exterior face nodes for each element */
  mesh->vmapM = (iint*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(iint));
  mesh->vmapP = (iint*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(iint));
  mesh->mapP  = (iint*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(iint));
  
  /* assume elements already connected */
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint f=0;f<mesh->Nfaces;++f){
      iint eP = mesh->EToE[e*mesh->Nfaces+f];
      iint fP = mesh->EToF[e*mesh->Nfaces+f];
      if(eP<0 || fP<0){ // fake connections for unconnected faces
	eP = e;
	fP = f;
      }
      /* for each node on this face find the neighbor node */
      for(iint n=0;n<mesh->Nfp;++n){
	iint  idM = mesh->faceNodes[f*mesh->Nfp+n] + e*mesh->Np;
	dfloat xM = mesh->x[idM];
	dfloat yM = mesh->y[idM];
	dfloat zM = mesh->z[idM];
	iint nP;
	
	iint  idP = findBestMatch(xM, yM, zM,
				  mesh->Nfp, 
				  mesh->faceNodes+fP*mesh->Nfp,
				  mesh->x+eP*mesh->Np,
				  mesh->y+eP*mesh->Np,
				  mesh->z+eP*mesh->Np, &nP,mesh->cubeFaceNumber[e],mesh->cubeFaceNumber[eP]);

	iint   id = mesh->Nfaces*mesh->Nfp*e + f*mesh->Nfp + n;
	mesh->vmapM[id] = idM;
	mesh->vmapP[id] = idP + eP*mesh->Np;
	mesh->mapP[id] = eP*mesh->Nfaces*mesh->Nfp + fP*mesh->Nfp + nP;
      }
    }
  }
}

