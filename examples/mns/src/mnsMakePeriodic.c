#include "mns.h"


dfloat findFaceMatch(mesh_t *mesh, int k1, int f1, int k2, int f2, dfloat xper, dfloat yper){
  dfloat xc1 = 0.0, yc1 = 0.0, xc2=0.0 , yc2=0.0;
  for(int n=0;n<mesh->Nfp;++n){
    int  id1 = mesh->faceNodes[f1*mesh->Nfp+n] + k1*mesh->Np;
    dfloat x1 = mesh->x[id1];
    dfloat y1 = mesh->y[id1];

    int  id2 = mesh->faceNodes[f2*mesh->Nfp+n] + k2*mesh->Np;
    dfloat x2 = mesh->x[id2];
    dfloat y2 = mesh->y[id2];

    xc1 += x1;
    yc1 += y1;
    xc2 += x2;
    yc2 += y2;
  }

  xc1 /= mesh->Nfp;
  yc1 /= mesh->Nfp;
  xc2 /= mesh->Nfp;
  yc2 /= mesh->Nfp;

  // Find distance
  const dfloat dx = pow(fabs(xc1-xc2)-xper,2) + pow(yc1-yc2,2);
  const dfloat dy = pow(xc1-xc2,2) + pow(fabs(yc1-yc2)-yper,2);
  //
  return mymin(dx,dy);
}

int findBestPeriodicMatch(dfloat xper, dfloat yper, dfloat x1, dfloat y1, 
       int Np2, int *nodeList, dfloat *x2, dfloat *y2, int *nP){
  
  int matchIndex = nodeList[0];
  dfloat mindist2 = pow(fabs(x1-x2[0])-xper,2) + pow(fabs(y1-y2[0])-yper,2);

  *nP = 0;
  for(int n=1;n<Np2;++n){
    
    /* next node */
    const int i2 = nodeList[n];

    /* distance between target and next node */
    const dfloat dist2 = pow(fabs(x1-x2[i2])-xper,2) + pow(fabs(y1-y2[i2])-yper,2);
    //const dfloat disty = pow(x1-x2[i2],2) + pow(fabs(y1-y2[i2])-yper,2);
    
    //const dfloat dist2 = mymin(distx, disty);
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
  void mnsMakePeriodic(mesh_t *mesh, dfloat xper, dfloat yper){ 

int Nconnected = 0;    
/* assume elements already connected */
for(int e=0;e<mesh->Nelements;++e){
  for(int f=0;f<mesh->Nfaces;++f){
    int e1 = mesh->EToE[e*mesh->Nfaces+f];
    int f1 = mesh->EToF[e*mesh->Nfaces+f];
    int bc1 = mesh->EToB[e*mesh->Nfaces+f];

    if( e1<0 && f1<0 && bc1<0){ 
      // printf("Catch a possible periodic connection on element: %d and face: %d",e,f);
      bool connected = false;
      /* assume elements already connected */
      for(int k=0;k<mesh->Nelements;++k){
        for(int l=0;l<mesh->Nfaces;++l){
          int e2  = mesh->EToE[k*mesh->Nfaces+l];
          int f2  = mesh->EToF[k*mesh->Nfaces+l];
          int bc2 = mesh->EToB[k*mesh->Nfaces+l];
          if( e2<0 && f2<0 && bc2<0 && k!=e){   
            //printf("With: %d and face: %d \n",k,l);
            dfloat dist = findFaceMatch(mesh,e,f,k,l,xper,yper);
            //printf("With: %d and face: %d and distance :%.4e \n",k,l,dist);
            if(dist<1e-10){
              // printf(" Connecting With: %d and face: %d \n",k,l);
              connected = true; Nconnected++; 
              mesh->EToE[e*mesh->Nfaces+f] = k;
              mesh->EToE[k*mesh->Nfaces+l] = e;
              //
              mesh->EToF[e*mesh->Nfaces+f] = l;
              mesh->EToF[k*mesh->Nfaces+l] = f;
              for(int n=0;n<mesh->Nfp;++n){
                int  id1M = mesh->faceNodes[f*mesh->Nfp+n] + e*mesh->Np;
                int  id2M = mesh->faceNodes[l*mesh->Nfp+n] + k*mesh->Np;

                dfloat x1M = mesh->x[id1M];
                dfloat y1M = mesh->y[id1M];
                //
                dfloat x2M = mesh->x[id2M];
                dfloat y2M = mesh->y[id2M];

                int   id1 = mesh->Nfaces*mesh->Nfp*e + f*mesh->Nfp + n;
                int   id2 = mesh->Nfaces*mesh->Nfp*k + l*mesh->Nfp + n;
                int n1P;
                int n2P;

                int e1P = k; 
                int f1P = l; 
                //
                int e2P = e; 
                int f2P = f; 

                int  id1P = findBestPeriodicMatch(xper,yper,x1M, y1M, 
                mesh->Nfp, 
                mesh->faceNodes+f1P*mesh->Nfp,
                mesh->x+e1P*mesh->Np,
                mesh->y+e1P*mesh->Np, &n1P);


                //No need to do this but just for double check, TODO: Change it later!!!!!
                int  id2P = findBestPeriodicMatch(xper,yper,x2M, y2M, 
                mesh->Nfp, 
                mesh->faceNodes+f2P*mesh->Nfp,
                mesh->x+e2P*mesh->Np,
                mesh->y+e2P*mesh->Np, &n2P);
                //Correct First Connection
                mesh->vmapM[id1] = id1M;
                mesh->vmapP[id1] = id1P + e1P*mesh->Np;

                //
                mesh->vmapM[id2] = id2M;
                mesh->vmapP[id2] = id2P + e2P*mesh->Np;

                //
                mesh->mapP[id1] = e1P*mesh->Nfaces*mesh->Nfp + f1P*mesh->Nfp + n1P;  

                //                              
                mesh->mapP[id2] = e2P*mesh->Nfaces*mesh->Nfp + f2P*mesh->Nfp + n2P;
              }
            } 
          }
        }
      }

      // if(!connected) 
      //   printf("missed \n");
    }
  }
} 

printf("Connected %d periodic pairs\n", Nconnected); 

}



      

  



