
#include "ellipticTri2D.h"

void ellipticBuildJacobiIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda,
                                   iint *BCType, dfloat **invDiagA,
                                   const char *options){
  
  iint size, rankM;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);
  

  if(!basis) { // default to degree N Lagrange basis
    basisNp = mesh->Np;
    basis = (dfloat*) calloc(basisNp*basisNp, sizeof(dfloat));
    for(iint n=0;n<basisNp;++n){
      basis[n+n*basisNp] = 1;
    }
  }

  iint diagNnum = basisNp*mesh->Nelements;

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(mesh->Nfaces*mesh->Nfp*mesh->Nfp,sizeof(dfloat));
  for (iint f=0;f<mesh->Nfaces;f++) {
    for (iint n=0;n<mesh->Nfp;n++) {
      iint fn = mesh->faceNodes[f*mesh->Nfp+n];
      
      for (iint m=0;m<mesh->Nfp;m++) {
        dfloat MSnm = 0;
        
        for (iint i=0;i<mesh->Np;i++){
          MSnm += mesh->MM[fn+i*mesh->Np]*mesh->LIFT[i*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+m];
        }
        
        MS[m+n*mesh->Nfp + f*mesh->Nfp*mesh->Nfp]  = MSnm;
      }
    }
  }
  
  // reset non-zero counter
  *invDiagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));

  dfloat *SM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

  // loop over all elements
  for(iint eM=0;eM<mesh->Nelements;++eM){
    
    iint vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat J = mesh->vgeo[vbase+JID];

    /* start with stiffness matrix  */
    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Np;++m){
        SM[n*mesh->Np+m]  = J*lambda*mesh->MM[n*mesh->Np+m];
        SM[n*mesh->Np+m] += J*drdx*drdx*mesh->Srr[n*mesh->Np+m];
        SM[n*mesh->Np+m] += J*drdx*dsdx*mesh->Srs[n*mesh->Np+m];
        SM[n*mesh->Np+m] += J*dsdx*drdx*mesh->Ssr[n*mesh->Np+m];
        SM[n*mesh->Np+m] += J*dsdx*dsdx*mesh->Sss[n*mesh->Np+m];
                                      
        SM[n*mesh->Np+m] += J*drdy*drdy*mesh->Srr[n*mesh->Np+m];
        SM[n*mesh->Np+m] += J*drdy*dsdy*mesh->Srs[n*mesh->Np+m];
        SM[n*mesh->Np+m] += J*dsdy*drdy*mesh->Ssr[n*mesh->Np+m];
        SM[n*mesh->Np+m] += J*dsdy*dsdy*mesh->Sss[n*mesh->Np+m];
      }
    }

    for (iint fM=0;fM<mesh->Nfaces;fM++) {
      
      // load surface geofactors for this face
      iint sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat hinv = mesh->sgeo[sid+IHID];
      dfloat penalty = tau*hinv; 
      
      iint eP = mesh->EToE[eM*mesh->Nfaces+fM];
      if (eP < 0) eP = eM;
      
      iint vbaseP = eP*mesh->Nvgeo;
      dfloat drdxP = mesh->vgeo[vbaseP+RXID];
      dfloat drdyP = mesh->vgeo[vbaseP+RYID];
      dfloat dsdxP = mesh->vgeo[vbaseP+SXID];
      dfloat dsdyP = mesh->vgeo[vbaseP+SYID];
      
      int bcD = 0, bcN =0;
      int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag
      iint bcType = 0;

      if(bc>0) bcType = BCType[bc];          //find its type (Dirichlet/Neumann)

      // this needs to be double checked (and the code where these are used)
      if(bcType==1){ // Dirichlet
        bcD = 1;
        bcN = 0;
      } else if(bcType==2){ // Neumann
        bcD = 0;
        bcN = 1;
      }
      
      // reset eP
      eP = mesh->EToE[eM*mesh->Nfaces+fM];

      // mass matrix for this face
      dfloat *MSf = MS+fM*mesh->Nfp*mesh->Nfp;

      // penalty term just involves face nodes
      for(iint n=0;n<mesh->Nfp;++n){
        for(iint m=0;m<mesh->Nfp;++m){
          iint nM = mesh->faceNodes[fM*mesh->Nfp+n];
          iint mM = mesh->faceNodes[fM*mesh->Nfp+m];
          
          // OP11 = OP11 + 0.5*( gtau*mmE )
          dfloat MSfnm = sJ*MSf[n*mesh->Nfp+m];
          SM[nM*mesh->Np+mM] += 0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
        }
      }

      // now add differential surface terms
      for(iint n=0;n<mesh->Nfp;++n){
        for(iint m=0;m<mesh->Np;++m){
          iint nM = mesh->faceNodes[fM*mesh->Nfp+n];
          
          for(iint i=0;i<mesh->Nfp;++i){
            iint iM = mesh->faceNodes[fM*mesh->Nfp+i];
            iint iP = mesh->vmapP[i + fM*mesh->Nfp+eM*mesh->Nfp*mesh->Nfaces]%mesh->Np;
              
            dfloat MSfni = sJ*MSf[n*mesh->Nfp+i]; // surface Jacobian built in
            
            dfloat DxMim = drdx*mesh->Dr[iM*mesh->Np+m] + dsdx*mesh->Ds[iM*mesh->Np+m];
            dfloat DyMim = drdy*mesh->Dr[iM*mesh->Np+m] + dsdy*mesh->Ds[iM*mesh->Np+m];

            // OP11 = OP11 + 0.5*( - mmE*Dn1)       
            SM[nM*mesh->Np+m] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
            SM[nM*mesh->Np+m] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;
          }
        }
      }

      for(iint n=0;n<mesh->Np;++n){
        for(iint m=0;m<mesh->Nfp;++m){
          iint mM = mesh->faceNodes[fM*mesh->Nfp+m];
          iint mP = mesh->vmapP[m + fM*mesh->Nfp+eM*mesh->Nfp*mesh->Nfaces]%mesh->Np;
          
          for(iint i=0;i<mesh->Nfp;++i){
            iint iM = mesh->faceNodes[fM*mesh->Nfp+i];  

            dfloat MSfim = sJ*MSf[i*mesh->Nfp+m];
            
            dfloat DxMin = drdx*mesh->Dr[iM*mesh->Np+n] + dsdx*mesh->Ds[iM*mesh->Np+n];
            dfloat DyMin = drdy*mesh->Dr[iM*mesh->Np+n] + dsdy*mesh->Ds[iM*mesh->Np+n];
          
            // OP11 = OP11 + (- Dn1'*mmE );
            SM[n*mesh->Np+mM] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
            SM[n*mesh->Np+mM] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
          }
        }
      }
    }
    // compute the diagonal entries
    for(iint i=0;i<basisNp;++i){
      dfloat val = 0;
      for(iint n=0;n<mesh->Np;++n){
        for(iint m=0;m<mesh->Np;++m){
          val += basis[n*mesh->Np+i]*SM[n*mesh->Np+m]*basis[m*mesh->Np+i];
        }
      }
  
      (*invDiagA)[eM*mesh->Np + i] = 1/val; //store the inverse diagonal entry
    }
  }

  free(SM);
}
