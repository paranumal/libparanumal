#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mesh3D.h"

int findBestMatch(dfloat x1, dfloat y1, dfloat z1,
		   int Np2, int *nodeList, dfloat *x2, dfloat *y2, dfloat *z2, int *nP){
  
  int matchIndex;
  dfloat mindist2=1e9;

  for(int n=0;n<Np2;++n){
    
    /* next node */
    const int i2 = nodeList[n];

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
  mesh->vmapM = (int*) calloc(mesh->NfpMax*mesh->Nfaces*mesh->Nelements, sizeof(int));
  mesh->vmapP = (int*) calloc(mesh->NfpMax*mesh->Nfaces*mesh->Nelements, sizeof(int));
  mesh->mapP  = (int*) calloc(mesh->NfpMax*mesh->Nfaces*mesh->Nelements, sizeof(int));
  
  dfloat xConnect[mesh->NpMax];
  dfloat yConnect[mesh->NpMax];
  dfloat zConnect[mesh->NpMax];

  /* assume elements already connected */
  for(int e=0;e<mesh->Nelements;++e){
    int N = mesh->N[e];

    /* for each node on this face find the neighbor node */
    for (int f=0;f<mesh->Nfaces;f++) {      
      int eP = mesh->EToE[e*mesh->Nfaces+f];
      int fP = mesh->EToF[e*mesh->Nfaces+f];
      
      if(eP<0 || fP<0){ // fake connections for unconnected faces
      	eP = e;
      	fP = f;
      }
      int NP = mesh->N[eP];

      for(int n=0;n<mesh->Nfp[N];++n){
        int   id = mesh->Nfaces*mesh->NfpMax*e + f*mesh->NfpMax + n;
        int  idM = mesh->faceNodes[N][f*mesh->Nfp[N]+n] + e*mesh->NpMax;
        mesh->vmapM[id] = idM;
      }


      //Construct node ordering which agrees with neighbour's degree
      int id = eP*mesh->Nverts;

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
      
      for(int n=0;n<mesh->Np[N];++n){ /* for each node */
        
        /* (r,s,t) coordinates of interpolation nodes*/
        dfloat rn = mesh->r[N][n]; 
        dfloat sn = mesh->s[N][n];
        dfloat tn = mesh->t[N][n];

        /* physical coordinate of interpolation node */
        xConnect[n] = -0.5*(1+rn+sn+tn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3 + 0.5*(1+tn)*xe4;
        yConnect[n] = -0.5*(1+rn+sn+tn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3 + 0.5*(1+tn)*ye4;
        zConnect[n] = -0.5*(1+rn+sn+tn)*ze1 + 0.5*(1+rn)*ze2 + 0.5*(1+sn)*ze3 + 0.5*(1+tn)*ze4;
      }

      for(int n=0;n<mesh->Nfp[N];++n){
        int  idM = mesh->faceNodes[N][f*mesh->Nfp[N]+n]; 
        dfloat xM = mesh->x[e*mesh->NpMax+idM];
        dfloat yM = mesh->y[e*mesh->NpMax+idM];
        dfloat zM = mesh->z[e*mesh->NpMax+idM];
        int nP=0;
        
        int  idP = findBestMatch(xM, yM, zM,
                mesh->Nfp[N], 
                mesh->faceNodes[N]+fP*mesh->Nfp[N],
                xConnect,
                yConnect,
                zConnect, &nP);

        int   id = mesh->Nfaces*mesh->NfpMax*e + f*mesh->NfpMax + n;
        mesh->vmapP[id] = idP + eP*mesh->NpMax;
        mesh->mapP[id] = eP*mesh->Nfaces*mesh->NfpMax + fP*mesh->NfpMax + nP;
      }
    }
  }
}

