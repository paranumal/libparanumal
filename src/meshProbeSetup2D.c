#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"

void meshProbeSetup2D(mesh2D *mesh, dfloat *pX, dfloat *pY){
 

 //dfloat tol     =-1.0e-8;        // set tolerance
 dfloat mindist = -1.e-9; // Set minum distance
 //

 mesh->probeN = 0; 
 
 dfloat *probeR           = (dfloat *) calloc(mesh->probeNTotal,sizeof(dfloat));
 dfloat *probeS           = (dfloat *) calloc(mesh->probeNTotal,sizeof(dfloat));
 mesh->probeElementIds  = (iint *)   calloc(mesh->probeNTotal,sizeof(iint));

// mesh->Nprobes =3; 

 dfloat *A        = (dfloat *) calloc((mesh->dim+1)*mesh->Nverts,sizeof(dfloat)); 
 iint   *IPIV     = (iint *)   calloc(mesh->Nverts,sizeof(iint)); 
 iint   *IPIV2    = (iint *)   calloc((mesh->Np+1),sizeof(iint)); 
 
 dfloat *b       = (dfloat *) calloc((mesh->dim+1)*mesh->probeNTotal,sizeof(dfloat));
 dfloat *q       = (dfloat *) calloc(mesh->Nverts*mesh->probeNTotal,sizeof(dfloat));

  for(iint n=0; n<mesh->Nverts;n++)
    IPIV[n] = 1; 

  for(iint n=0; n<(mesh->Np+1);n++)
    IPIV2[n] = 1; 


 iint N    = (mesh->dim+1); // A->Nrwos
 iint NRHS =  mesh->probeNTotal; // B->Ncolumns
 iint LDA  = N; 
 iint LDB  = (mesh->dim + 1); // B->Nrows
 iint INFO;

  // for(iint n=0; n<mesh->mesh->probeNTotal; n++){
  // // Coordinates of probe
  // printf("Probe %d  pX: %g pY:%g \n",n, pX[n], pY[n]);
  // }

  // Assumes each probe is in one element, may change later 
  for(iint n=0; n<mesh->probeNTotal; n++){
    // Coordinates of probe
    b[n*mesh->probeN + 0] = 1.0; 
    b[n*mesh->probeN + 1] = pX[n]; 
    b[n*mesh->probeN + 2] = pY[n]; 
  }
  
  
  //
  for (iint e=0;e<mesh->Nelements;e++) {
  // Create A[1 vx vy]
    for (iint n=0;n<mesh->Nverts;n++) {
      dfloat vx = mesh->EX[e*mesh->Nverts+n];
      dfloat vy = mesh->EY[e*mesh->Nverts+n];
      //
      A[n*mesh->Nverts + 0] = 1.0;
      A[n*mesh->Nverts + 1] = vx;
      A[n*mesh->Nverts + 2] = vy;
    } 


    for(iint n=0; n<mesh->probeNTotal*(mesh->dim+1); n++)
      q[n] = b[n]; 
    
  
  // Find barycentric coordinates
  // q = A^-1*b

   if(sizeof(dfloat) == sizeof(double))
    dgesv_(&N,&NRHS,A,&LDA,IPIV,q,&LDB,&INFO);

   if(INFO)
    printf("DGSEV error: %d \n", INFO);

    // if(sizeof(dfloat) == sizeof(float))
    // sgesv_(&N,&NRHS,A,&LDA,IPIV,q,&LDB,&INFO);

    // Check all non-negative barycentric coordinates 
    // Assumes a probe can be represented by single element!!!!
    for(iint n=0; n<mesh->probeNTotal; n++){

       dfloat qmin = q[n*mesh->probeN + 0]; 
       for(iint i =1; i<(mesh->dim+1); i++)
          qmin = mymin(qmin, q[n*mesh->probeNTotal + i]);
        
        // Catch the element
        if(qmin>mindist){
         mesh->probeN++;
         // hold element ids
         mesh->probeElementIds[n] = e; 
         // hold local r,s coordinates
         dfloat l1 =  q[n*mesh->probeNTotal + 2]; 
         dfloat l2 =  q[n*mesh->probeNTotal + 0]; 
         dfloat l3 =  q[n*mesh->probeNTotal + 1]; 

         probeR[n] = 2.*l3 -1.; 
         probeS[n] = 2.*l1 -1.;

         printf("element: %d probe %d qmin:%.5e R: %.5e S:%.5e\n", e, n, qmin, probeR[n],probeS[n]); 

        }
    }
  }
 
  if(mesh->probeN){
  // Compute Vandermonde Matrix and invert  it
  dfloat *V       = (dfloat *) calloc(mesh->Np* (mesh->N+1)*(mesh->N+2)/2, sizeof(dfloat));
  meshVandermonde2D(mesh->N, mesh->Np, mesh->r, mesh->s, V);

  //
  N    = mesh->Np; 
  iint LWORK = mesh->Np*mesh->Np;
  dfloat *WORK = (dfloat *) calloc(LWORK, sizeof(dfloat));
  dgetrf_(&N,&N,V,&N,IPIV2,&INFO);
  dgetri_(&N,V,&N,IPIV2,WORK,&LWORK,&INFO);
  if(INFO)
    printf("DGE_TRI/TRF error: %d \n", INFO);

  
  // Compute Vandermonde matrix of probes
  dfloat *Vprobe = (dfloat *) calloc(mesh->probeN*mesh->Np,sizeof(dfloat));
  meshVandermonde2D(mesh->N, mesh->probeN, probeR, probeS, Vprobe);
  

  mesh->probeI = (dfloat *) calloc(mesh->probeN*mesh->Np, sizeof(dfloat));

  for(iint r=0; r<mesh->probeN; r++){
    for(iint c=0; c<mesh->Np; c++){
      dfloat s = Vprobe[r*mesh->Np+0]*V[0*mesh->Np + c];
      for(iint i=1; i<mesh->Np; i++){
        s += Vprobe[r*mesh->Np+i]*V[i*mesh->Np + c];
      }
      mesh->probeI[r*mesh->Np + c] = s;
    }

  }

  free(V);
  free(WORK);

}

free(IPIV);
free(IPIV2);
free(probeR);
free(probeS);
free(A);
free(b);
free(q);

}


void meshVandermonde2D(iint N, iint Npoints, dfloat *r, dfloat *s, dfloat *V){

// First convert to ab coordinates
  dfloat *a = (dfloat *) calloc(Npoints, sizeof(dfloat));
  dfloat *b = (dfloat *) calloc(Npoints, sizeof(dfloat));
  for(iint n=0; n<Npoints; n++){

    if(fabs(s[n]-1.0)>1e-8)
      a[n] = 2.0*(1.+r[n])/(1.0-s[n])-1.0;
    else
      a[n] = -1.0; 

    b[n] = s[n];

  }
  
  iint sk=0;

  iint Np = (N+1)*(N+2)/2; 

  for(iint i=0; i<=N; i++){
    for(iint j=0; j<=N-i; j++){
      for(iint n=0; n<Npoints; n++){
        V[n*Np + sk] = meshSimplex2D(a[n], b[n], i, j);
      }
    sk++;
    }
  }


  free(a);
  free(b); 

}



dfloat meshSimplex2D(dfloat a, dfloat b, iint i, iint j){
// 
 dfloat p1 = meshJacobiP(a,0,0,i);
 dfloat p2 = meshJacobiP(b,2*i+1,0,j);
 dfloat P = sqrt(2.0)*p1*p2*pow(1-b,i);

 return P; 
}


dfloat meshJacobiP(dfloat a, dfloat alpha, dfloat beta, iint N){

dfloat ax = a; 

dfloat *P = (dfloat *) calloc((N+1), sizeof(dfloat));

// Zero order
dfloat gamma0 = pow(2,(alpha+beta+1))/(alpha+beta+1)*meshFactorial(alpha)*meshFactorial(beta)/meshFactorial(alpha+beta);
dfloat p0     = 1.0/sqrt(gamma0);

  if (N==0){ free(P); return p0;}
P[0] = p0; 

// first order
dfloat gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
dfloat p1     = ((alpha+beta+2)*ax/2 + (alpha-beta)/2)/sqrt(gamma1);
  if (N==1){free(P); return p1;} 

P[1] = p1;

 /// Repeat value in recurrence.
  dfloat aold = 2/(2+alpha+beta)*sqrt((alpha+1.)*(beta+1.)/(alpha+beta+3.));
  /// Forward recurrence using the symmetry of the recurrence.
  for(iint i=1;i<=N-1;++i){
    dfloat h1 = 2.*i+alpha+beta;
    dfloat anew = 2./(h1+2.)*sqrt( (i+1.)*(i+1.+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3));
    dfloat bnew = -(alpha*alpha-beta*beta)/h1/(h1+2);
    P[i+1] = 1./anew*( -aold*P[i-1] + (ax-bnew)*P[i]);
    aold =anew;
  }
  
 dfloat pN = P[N]; 
 free(P);
   return pN;

}


dfloat meshFactorial(int n){

  if(n==0)
    return 1;
  else
    return n*meshFactorial(n-1);
}
