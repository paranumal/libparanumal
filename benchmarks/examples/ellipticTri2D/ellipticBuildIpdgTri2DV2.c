
#include "ellipticTri2D.h"

typedef struct{

  iint row;
  iint col;
  iint ownerRank;
  dfloat val;

} nonZero_t;

void ellipticBuildIpdgTri2D(mesh2D *mesh, dfloat lambda, nonZero_t **A, iint *nnzA, const char *options){

  iint size, rankM;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);

  /* find number of elements on each rank */
  iint *rankNelements = (iint*) calloc(size, sizeof(iint));
  iint *rankStarts = (iint*) calloc(size+1, sizeof(iint));
  MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT,
		rankNelements, 1, MPI_IINT, MPI_COMM_WORLD);
  for(iint r=0;r<size;++r){
    rankStarts[r+1] = rankStarts[r]+(rankNelements[r])*mesh->Np;
  }
  
  /* do a halo exchange of local node numbers */

  iint nnzLocalBound = mesh->Np*mesh->Np*(1+mesh->Nfaces)*mesh->Nelements;

  *A = (nonZero_t*) calloc(nnzLocalBound, sizeof(nonZero_t));

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;
  dfloat tau = 2; // hackery

  dfloat *BM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

  dfloat *qmP = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *qmM = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *ndotgradqmM = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *ndotgradqmP = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  for (iint f=0;f<mesh->Nfaces;f++) {
    for (iint n=0;n<mesh->Np;n++) {
      for (iint m=0;m<mesh->Nfp;m++) {
	dfloat MSnm = 0;
		  
	for (iint i=0;i<mesh->Np;i++) 
	  MSnm += mesh->MM[n+i*mesh->Np]*mesh->LIFT[i*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+m];
	  
	MS[m+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]  = MSnm;
      }
    }
  }

  // DrT*MS, DsT*MS
  dfloat *DrTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  dfloat *DsTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  for (iint f=0;f<mesh->Nfaces;f++) {
    for (iint n=0;n<mesh->Np;n++) {
      for (iint i=0;i<mesh->Nfp;i++) {
        DrTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        DsTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        for (iint m=0;m<mesh->Np;m++) {
          DrTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] 
            += mesh->Dr[n+m*mesh->Np]*MS[i+m*mesh->Nfp+f*mesh->Nfp*mesh->Np];
          DsTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] 
            += mesh->Ds[n+m*mesh->Np]*MS[i+m*mesh->Nfp+f*mesh->Nfp*mesh->Np];
        }
      }
    }
  }

#if 0
  for (iint n=0;n<mesh->Np;n++) {
    iint m=0;
    for(iint j=0;j<mesh->N+1;++j){
      for(iint i=0;i<mesh->N+1-j;++i){
	// try rescaling
	mesh->V[n*mesh->Np+m] *= (i+j+1)*(i+j+1);
	++m;
      }
    }
  }
#endif
  
  // reset non-zero counter
  int nnz = 0;

  // loop over all elements
  for(iint eM=0;eM<mesh->Nelements;++eM){

    dfloat *BP = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Np,sizeof(dfloat));
    
    iint gbase = eM*mesh->Nggeo;
    dfloat Grr = mesh->ggeo[gbase+G00ID];
    dfloat Grs = mesh->ggeo[gbase+G01ID];
    dfloat Gss = mesh->ggeo[gbase+G11ID];
    dfloat J   = mesh->ggeo[gbase+GWJID];

    /* start with stiffness matrix  */
    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Np;++m){
	iint id = m + n*mesh->Np;
        BM[id]  = J*lambda*mesh->MM[id];
        BM[id] += Grr*mesh->Srr[id];
        BM[id] += Grs*mesh->Srs[id];
        BM[id] += Grs*mesh->Ssr[id];
        BM[id] += Gss*mesh->Sss[id];
      }
    }

    iint vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];

    for (iint m=0;m<mesh->Np;m++) {
      for (iint fM=0;fM<mesh->Nfaces;fM++) {
        // load surface geofactors for this face
        iint sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
        dfloat nx = mesh->sgeo[sid+NXID];
        dfloat ny = mesh->sgeo[sid+NYID];
        dfloat sJ = mesh->sgeo[sid+SJID];
        dfloat hinv = mesh->sgeo[sid+IHID];

        iint eP = mesh->EToE[eM*mesh->Nfaces+fM];
        if (eP < 0) eP = eM;
        iint vbaseP = eP*mesh->Nvgeo;
        dfloat drdxP = mesh->vgeo[vbaseP+RXID];
        dfloat drdyP = mesh->vgeo[vbaseP+RYID];
        dfloat dsdxP = mesh->vgeo[vbaseP+SXID];
        dfloat dsdyP = mesh->vgeo[vbaseP+SYID];
      
        // extract trace nodes
        for (iint i=0;i<mesh->Nfp;i++) {
          // double check vol geometric factors are in halo storage of vgeo
          iint idM    = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+i;
          iint vidM   = mesh->faceNodes[i+fM*mesh->Nfp];
          iint vidP   = mesh->vmapP[idM]%mesh->Np; // only use this to identify location of positive trace vgeo

          qmM[i] =0;
          if (vidM == m) qmM[i] =1;
          qmP[i] =0;
          if (vidP == m) qmP[i] =1;
          
          ndotgradqmM[i] = (nx*drdx+ny*drdy)*mesh->Dr[m+vidM*mesh->Np]   +(nx*dsdx+ny*dsdy)*mesh->Ds[m+vidM*mesh->Np];
          ndotgradqmP[i] = (nx*drdxP+ny*drdyP)*mesh->Dr[m+vidP*mesh->Np] +(nx*dsdxP+ny*dsdyP)*mesh->Ds[m+vidP*mesh->Np];
        }

        dfloat penalty = tau*(mesh->N+1)*(mesh->N+1)*hinv; 
        eP = mesh->EToE[eM*mesh->Nfaces+fM];
        for (iint n=0;n<mesh->Np;n++) {
	  iint id = m+n*mesh->Np;
          for (iint i=0;i<mesh->Nfp;i++) {
            BM[id] += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmM[i];
            BM[id] +=
	      -0.5*sJ*(nx*drdx+ny*drdy)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
    	      -0.5*sJ*(nx*dsdx+ny*dsdy)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]; 
            BM[id] += +0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*penalty*qmM[i];
          }

          dfloat AnmP = 0;
          if (eP < 0) {
            for (iint i=0;i<mesh->Nfp;i++) {
              BM[id] += +0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmM[i];
              BM[id] +=
		+0.5*sJ*(nx*drdx+ny*drdy)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
		+0.5*sJ*(nx*dsdx+ny*dsdy)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]; 
              BM[id] += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*penalty*qmM[i];
            }
          } else {
            for (iint i=0;i<mesh->Nfp;i++) {
              BP[fM*mesh->Np*mesh->Np+id] += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmP[i];
              BP[fM*mesh->Np*mesh->Np+id] +=
		+0.5*sJ*(nx*drdx+ny*drdy)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmP[i]
		+0.5*sJ*(nx*dsdx+ny*dsdy)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmP[i]; 
              BP[fM*mesh->Np*mesh->Np+id] += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*penalty*qmP[i];
            }
          }
	}
      }
    }

    for (iint n=0;n<mesh->Np;n++) {
      for (iint m=0;m<mesh->Np;m++) {
#if 1
	dfloat Anm = BM[m+n*mesh->Np];
#else
	// use this to transform the diagonal blocks
	dfloat Anm = 0;
	for(iint j=0;j<mesh->Np;++j){
	  for(iint i=0;i<mesh->Np;++i){
	    // V_{jn} B-_{ji} V_{im}
	    Anm += mesh->V[j*mesh->Np+n]*BM[j*mesh->Np+i]*mesh->V[i*mesh->Np+m];
	  }
	}
#endif
        if(fabs(Anm)>tol){
          (*A)[nnz].row = n + eM*mesh->Np + rankStarts[rankM];
          (*A)[nnz].col = m + eM*mesh->Np + rankStarts[rankM]; // this is wrong
          (*A)[nnz].val = Anm;
          (*A)[nnz].ownerRank = rankM;
          ++nnz;
        }
      } 
    }

    for (iint fM=0;fM<mesh->Nfaces;fM++) {
      iint eP = mesh->EToE[eM*mesh->Nfaces+fM];
      if(eP>=0){
	for (iint n=0;n<mesh->Np;n++) {
	  for (iint m=0;m<mesh->Np;m++) {
#if 1
	    dfloat AnmP = BP[fM*mesh->Np*mesh->Np+n*mesh->Np+m];
#else
	    // use this to transform the off diagonal blocks
	    dfloat AnmP = 0;
	    for(iint j=0;j<mesh->Np;++j){
	      for(iint i=0;i<mesh->Np;++i){
		// V_{jn} B+_{ji} V_{im}
		AnmP +=
		  mesh->V[j*mesh->Np+n]*BP[fM*mesh->Np*mesh->Np+j*mesh->Np+i]*mesh->V[i*mesh->Np+m];
	      }
	    }
#endif
	    if(fabs(AnmP)>tol){
	      // remote info
	      iint rankP = mesh->EToP[eM*mesh->Nfaces+fM]; 
	      if (rankP<0) rankP = rankM;
	      (*A)[nnz].row = n + eM*mesh->Np + rankStarts[rankM];
	      (*A)[nnz].col = m + eP*mesh->Np + rankStarts[rankP]; // this is wrong
	      (*A)[nnz].val = AnmP;
	      (*A)[nnz].ownerRank = rankM;
	      ++nnz;
	    }
	  }
	}
      }
    }
    
    free(BP);
  }
  
  // free up unused storage
  *A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  
  *nnzA = nnz;

  FILE *fp = fopen("sparseA.m", "w");
  fprintf(fp, "spA = [ \n");
  for(iint n=0;n<nnz;++n){
    fprintf(fp, "%d %d %17.15g\n", (*A)[n].row+1, (*A)[n].col+1, (*A)[n].val);
  }
  fprintf(fp, "];\n");
  fprintf(fp, "A = spconvert(spA);\n");
  fprintf(fp, "B = (A-transpose(A));\n");
  fclose(fp);
  free(BM);  free(MS);
  free(DrTMS); free(DsTMS);  
  
  free(qmM); free(qmP);
  free(ndotgradqmM); free(ndotgradqmP);
}
