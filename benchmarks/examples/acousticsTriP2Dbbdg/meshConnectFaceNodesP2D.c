#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mesh2D.h"

iint findBestMatch(dfloat x1, dfloat y1, 
		   iint Np2, iint *nodeList, dfloat *x2, dfloat *y2, int *nP){
  
  iint matchIndex = nodeList[0];
  dfloat mindist2 = pow(x1-x2[nodeList[0]],2) + pow(y1-y2[nodeList[0]],2);

  for(int n=1;n<Np2;++n){
    
    /* next node */
    const iint i2 = nodeList[n];

    /* distance between target and next node */
    const dfloat dist2 = pow(x1-x2[i2],2) + pow(y1-y2[i2],2);
    
    /* if next node is closer to target update match */
    if(dist2<mindist2){
      mindist2 = dist2;
      matchIndex = i2;
      *nP = n;
    }
  }
  if(mindist2>1e-3) printf("arggh - bad match: x,y=%g,%g\n", x1,y1);
  return matchIndex;
}
		

// serial face-node to face-node connection
void meshConnectFaceNodesP2D(mesh2D *mesh){
  
  /* volume indices of the interior and exterior face nodes for each element */
  mesh->vmapM = (iint*) calloc(mesh->NfpMax*mesh->Nfaces*mesh->Nelements, sizeof(iint));
  mesh->vmapP = (iint*) calloc(mesh->NfpMax*mesh->Nfaces*mesh->Nelements, sizeof(iint));
  mesh->mapP  = (iint*) calloc(mesh->NfpMax*mesh->Nfaces*mesh->Nelements, sizeof(iint));
  
  dfloat xConnect[mesh->NpMax];
  dfloat yConnect[mesh->NpMax];

  /* assume elements already connected */
  for(iint e=0;e<mesh->Nelements;++e){
    iint N = mesh->N[e];
    for(iint f=0;f<mesh->Nfaces;++f){
      iint eP = mesh->EToE[e*mesh->Nfaces+f];
      iint fP = mesh->EToF[e*mesh->Nfaces+f];
 
      if(eP<0 || fP<0){ // fake connections for unconnected faces
         eP = e;
         fP = f;
      }
      iint NP = mesh->N[eP];

      /* for each node on this face find the neighbor node */
      for(iint n=0;n<mesh->Nfp[N];++n){
        iint   id = e*mesh->Nfaces*mesh->NfpMax + f*mesh->NfpMax + n;     //face node index
      	iint  idM = e*mesh->NpMax + mesh->faceNodes[N][f*mesh->Nfp[N]+n]; //global index
      	mesh->vmapM[id] = idM;
      }
    
      //Construct node ordering of neighbour which agrees with local element's degree
      iint id = eP*mesh->Nverts;

      dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
      dfloat xe2 = mesh->EX[id+1];
      dfloat xe3 = mesh->EX[id+2];

      dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
      dfloat ye2 = mesh->EY[id+1];
      dfloat ye3 = mesh->EY[id+2];
      
      for(iint n=0;n<mesh->Np[N];++n){ /* for each node */
        
        /* (r,s,t) coordinates of interpolation nodes*/
        dfloat rn = mesh->r[N][n]; 
        dfloat sn = mesh->s[N][n];

        /* physical coordinate of interpolation node */
        xConnect[n] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
        yConnect[n] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
      }

      for(iint n=0;n<mesh->Nfp[N];++n){
        iint  idM = mesh->faceNodes[N][f*mesh->Nfp[N]+n]; 
        dfloat xM = mesh->x[e*mesh->NpMax+idM];
        dfloat yM = mesh->y[e*mesh->NpMax+idM];
        iint nP = 0;
        
        iint  idP = findBestMatch(xM, yM,
                mesh->Nfp[N], 
                mesh->faceNodes[N]+fP*mesh->Nfp[N],
                xConnect,
                yConnect, &nP);

        iint   id = mesh->Nfaces*mesh->NfpMax*e + f*mesh->NfpMax + n;
        mesh->vmapP[id] = idP + eP*mesh->NpMax;
        mesh->mapP[id] = eP*mesh->Nfaces*mesh->NfpMax + fP*mesh->NfpMax + nP;
      }
    }
  }
}
