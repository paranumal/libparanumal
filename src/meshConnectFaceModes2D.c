#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mesh2D.h"


// serial face-mode to face-mode connection
void meshConnectFaceModes2D(mesh2D *mesh, int *faceModes, dfloat *V){

  /* volume indices of the interior and exterior face modes for each element */
  mesh->mmapM = (iint*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(iint));
  mesh->mmapP = (iint*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(iint));
  mesh->mmapS = (iint*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(iint));

  dfloat *VM = (dfloat *) calloc(mesh->Np,sizeof(dfloat));
  dfloat *VP = (dfloat *) calloc(mesh->Np,sizeof(dfloat));

  /* assume elements already connected */
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint f=0;f<mesh->Nfaces;++f){
      iint eP = mesh->EToE[e*mesh->Nfaces+f];
      iint fP = mesh->EToF[e*mesh->Nfaces+f];
      if(eP<0 || fP<0){ // fake connections for unconnected faces
        eP = e;
        fP = f;
      }
      
      /* for each mode on this face find the neighbor mode */
      for(int n=0;n<mesh->Nfp;++n){
        int m = faceModes[n+f*mesh->Nfp]; //get face mode number

        for (int i=0;i<mesh->Nfp;i++) {
          int k = mesh->faceNodes[i+f*mesh->Nfp];
          VM[i] = V[m+k*mesh->Np]; //evaluate mode at WB nodes on face  
        }
        
        dfloat mindist = 1E9;
        int s; 
        int mMatch;
        for (int nP=0;nP<mesh->Nfp;nP++) {
          //test the modes on face fP
          int mP = faceModes[nP+fP*mesh->Nfp];

          for (int i=0;i<mesh->Nfp;i++) {
            //get neighbouring node
            iint id = i+f*mesh->Nfp+e*mesh->Nfp*mesh->Nfaces;
            int k = mesh->vmapP[id]%mesh->Np;

            VP[i] = V[mP+k*mesh->Np]; //evaluate mode at WB nodes on face     
          }

          dfloat dist1=0, dist2=0;
          for (int i=0;i<mesh->Nfp;i++){
            dist1 += pow(VM[i]-VP[i],2);
            dist2 += pow(VM[i]+VP[i],2);
          }
          dist1 = sqrt(dist1); 
          dist2 = sqrt(dist2);

          /* if next node is closer to target update match */
          if(dist1<mindist){
            mindist = dist1;
            mMatch = mP;
            s = 1;
          }
          if(dist2<mindist){
            mindist = dist2;
            mMatch = mP;
            s = -1;
          }
        }
        if(mindist2>1e-3) printf("arggh - bad match: e=%d,f=%d, mode=%d\n", e,f, m);

        iint id  = mesh->Nfaces*mesh->Nfp*e + f*mesh->Nfp + n;
        iint idM = faceModes[f*mesh->Nfp+n] + e*mesh->Np;
        iint idP = mMatch + eP*mesh->Np;

        mesh->mmapM[id] = idM;
        mesh->mmapP[id] = idP;
        mesh->mmapS[id] = s;
      }
    }
  }
}