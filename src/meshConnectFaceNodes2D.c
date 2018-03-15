#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mesh2D.h"

int findBestMatch(dfloat x1, dfloat y1,
                   int Np2, int *nodeList, dfloat *x2, dfloat *y2, int *nP){

  int matchIndex = nodeList[0];
  dfloat mindist2 = pow(x1-x2[nodeList[0]],2) + pow(y1-y2[nodeList[0]],2);

  *nP = 0;
  for(int n=1;n<Np2;++n){

    /* next node */
    const int i2 = nodeList[n];

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
void meshConnectFaceNodes2D(mesh2D *mesh){

  /* volume indices of the interior and exterior face nodes for each element */
  mesh->vmapM = (int*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(int));
  mesh->vmapP = (int*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(int));
  mesh->mapP  = (int*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(int));

  /* assume elements already connected */
  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      int eP = mesh->EToE[e*mesh->Nfaces+f];
      int fP = mesh->EToF[e*mesh->Nfaces+f];
      if(eP<0 || fP<0){ // fake connections for unconnected faces
        eP = e;
        fP = f;
      }
      /* for each node on this face find the neighbor node */
      for(int n=0;n<mesh->Nfp;++n){
        int  idM = mesh->faceNodes[f*mesh->Nfp+n] + e*mesh->Np;
        dfloat xM = mesh->x[idM];
        dfloat yM = mesh->y[idM];
        int   id = mesh->Nfaces*mesh->Nfp*e + f*mesh->Nfp + n;
        int nP;

        int  idP = findBestMatch(xM, yM,
                                  mesh->Nfp,
                                  mesh->faceNodes+fP*mesh->Nfp,
                                  mesh->x+eP*mesh->Np,
                                  mesh->y+eP*mesh->Np, &nP);

        mesh->vmapM[id] = idM;
        mesh->vmapP[id] = idP + eP*mesh->Np;
        mesh->mapP[id] = eP*mesh->Nfaces*mesh->Nfp + fP*mesh->Nfp + nP;

      }
    }
  }
}

//      printf("connecting (%d,%d) to (%d,%d) [ vmapM %d to vmapP %d ]\n",
//             e,f,eP,fP, mesh->vmapM[id], mesh->vmapP[id]);
