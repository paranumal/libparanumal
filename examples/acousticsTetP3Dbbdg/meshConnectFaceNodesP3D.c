#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mesh3D.h"

iint findBestMatch(dfloat x1, dfloat y1, dfloat z1,
		   iint Np2, iint *nodeList, dfloat *x2, dfloat *y2, dfloat *z2, int *nP){
  
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
  if(mindist2>1e-3) printf("arggh - bad match: x,y,z=%g,%g,%g\n", x1,y1,z1);

  return matchIndex;
}
		

// serial face-node to face-node connection
void meshConnectFaceNodesP3D(mesh3D *mesh){
  
  /* volume indices of the interior and exterior face nodes for each element */
  mesh->vmapM = (iint*) calloc(mesh->NfpMax*mesh->Nfaces*mesh->Nelements, sizeof(iint));
  mesh->vmapP = (iint*) calloc(mesh->NfpMax*mesh->Nfaces*mesh->Nelements, sizeof(iint));
  mesh->mapP  = (iint*) calloc(mesh->NfpMax*mesh->Nfaces*mesh->Nelements, sizeof(iint));
  
  dfloat xConnect[mesh->NpMax];
  dfloat yConnect[mesh->NpMax];
  dfloat zConnect[mesh->NpMax];

  /* assume elements already connected */
  for(iint e=0;e<mesh->Nelements;++e){
    iint N = mesh->N[e];

    /* for each node on this face find the neighbor node */
    for (iint f=0;f<mesh->Nfaces;f++) {      
      iint eP = mesh->EToE[e*mesh->Nfaces+f];
      iint fP = mesh->EToF[e*mesh->Nfaces+f];
      
      if(eP<0 || fP<0){ // fake connections for unconnected faces
      	eP = e;
      	fP = f;
      }
      iint NP = mesh->N[eP];

      for(iint n=0;n<mesh->Nfp[N];++n){
        iint  idM = mesh->faceNodes[N][f*mesh->Nfp[N]+n] + e*mesh->NpMax;
        iint   id = mesh->Nfaces*mesh->NfpMax*e + f*mesh->NfpMax + n;
        mesh->vmapM[id] = idM;
      }


      //Construct node ordering which agrees with neighbour's degree
      iint id = e*mesh->Nverts;

      dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
      dfloat xe2 = mesh->EX[id+1];
      dfloat xe3 = mesh->EX[id+2];
      dfloat xe4 = mesh->EX[id+3];

      dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
      dfloat ye2 = mesh->EY[id+1];
      dfloat ye3 = mesh->EY[id+2];
      dfloat ye4 = mesh->EY[id+3];
      
      dfloat ze1 = mesh->EZ[id+0]; /* z-coordinates of vertices */
      dfloat ze2 = mesh->EZ[id+1];
      dfloat ze3 = mesh->EZ[id+2];
      dfloat ze4 = mesh->EZ[id+3];
      
      for(iint n=0;n<mesh->Np[NP];++n){ /* for each node */
        
        /* (r,s,t) coordinates of interpolation nodes*/
        dfloat rn = mesh->r[NP][n]; 
        dfloat sn = mesh->s[NP][n];
        dfloat tn = mesh->t[NP][n];

        /* physical coordinate of interpolation node */
        xConnect[n] = -0.5*(1+rn+sn+tn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3 + 0.5*(1+tn)*xe4;
        yConnect[n] = -0.5*(1+rn+sn+tn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3 + 0.5*(1+tn)*ye4;
        zConnect[n] = -0.5*(1+rn+sn+tn)*ze1 + 0.5*(1+rn)*ze2 + 0.5*(1+sn)*ze3 + 0.5*(1+tn)*ze4;
      }

      for(iint n=0;n<mesh->Nfp[NP];++n){
        iint  idM = mesh->faceNodes[NP][f*mesh->Nfp[NP]+n]; 
        dfloat xM = xConnect[idM];
        dfloat yM = yConnect[idM];
        dfloat zM = zConnect[idM];
        iint nP;
        
        iint  idP = findBestMatch(xM, yM, zM,
                mesh->Nfp[NP], 
                mesh->faceNodes[NP]+fP*mesh->Nfp[NP],
                mesh->x+eP*mesh->NpMax,
                mesh->y+eP*mesh->NpMax,
                mesh->z+eP*mesh->NpMax, &nP);

        iint   id = mesh->Nfaces*mesh->NfpMax*e + f*mesh->NfpMax + n;
        mesh->vmapP[id] = idP + eP*mesh->NpMax;
        mesh->mapP[id] = eP*mesh->Nfaces*mesh->NfpMax + fP*mesh->NfpMax + nP;
      }
    }
  }
}

