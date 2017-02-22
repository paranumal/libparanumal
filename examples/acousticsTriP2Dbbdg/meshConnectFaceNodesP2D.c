#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mesh2D.h"

iint findBestMatch(dfloat x1, dfloat y1, 
		   iint Np2, iint *nodeList, dfloat *x2, dfloat *y2){
  
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
  
  /* assume elements already connected */
  for(iint e=0;e<mesh->Nelements;++e){
    iint N = mesh->N[e];
    for(iint f=0;f<mesh->Nfaces;++f){
      iint eP = mesh->EToE[e*mesh->Nfaces+f];
      iint fP = mesh->EToF[e*mesh->Nfaces+f];
      iint NP = mesh->N[eP];
 
      if(eP<0 || fP<0){ // fake connections for unconnected faces
         eP = e;
         fP = f;
         NP = N;
      }

      /* for each node on this face find the neighbor node */
      for(iint n=0;n<mesh->Nfp[N];++n){
        iint   id = e*mesh->Nfaces*mesh->NfpMax + f*mesh->NfpMax + n;     //face node index
      	iint  idM = e*mesh->NpMax + mesh->faceNodes[N][f*mesh->Nfp[N]+n]; //global index
      	mesh->vmapM[id] = idM;
      }
    
      if (eP<0 || fP <0){
        // make vmapP map to the same point as vmapM 
        for(iint n=0;n<mesh->Nfp[N];++n){
          iint   id = e*mesh->Nfaces*mesh->NfpMax + f*mesh->NfpMax + n;     //face node index
          iint  idM = e*mesh->NpMax + mesh->faceNodes[N][f*mesh->Nfp[N]+n]; //global index
          mesh->vmapP[id] = idM;
        }
      } else {
        for (iint n=0;n<mesh->Nfp[NP];++n){
          iint   id = e*mesh->Nfaces*mesh->NfpMax + f*mesh->NfpMax + n;     //face node index
          iint idP;
          // The ordering of nodes along side f=2 is reversed, which makes this part awkward 
          if (f==0 || f==1) {
            if (fP==0 || fP==1) {
              idP = eP*mesh->NpMax + mesh->faceNodes[NP][fP*mesh->Nfp[NP]+mesh->Nfp[NP]-1-n]; //global index    
            } else {
              idP = eP*mesh->NpMax + mesh->faceNodes[NP][fP*mesh->Nfp[NP]+n]; 
            }
          } else {
            if (fP==2) {
              idP = eP*mesh->NpMax + mesh->faceNodes[NP][fP*mesh->Nfp[NP]+mesh->Nfp[NP]-1-n]; 
            } else {
              idP = eP*mesh->NpMax + mesh->faceNodes[NP][fP*mesh->Nfp[NP]+n]; 
            }
          }
          mesh->vmapP[id] = idP;
        }
      }
    }
  }
}

//	printf("connecting (%d,%d) to (%d,%d) [ vmapM %d to vmapP %d ]\n", 
//	       e,f,eP,fP, mesh->vmapM[id], mesh->vmapP[id]);
