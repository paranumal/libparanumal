
#include "ellipticTri2D.h"

int parallelCompareRowColumn(const void *a, const void *b);

void ellipticBuildBRdgTri2D(mesh2D *mesh, dfloat tau, dfloat lambda, iint *BCType, nonZero_t **A,
                            iint *nnzA, iint *globalStarts, const char *options){

  iint size, rankM;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);

  int Np = mesh->Np;
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;
  iint Nelements = mesh->Nelements;

  iint Nnum = Np*Nelements;

  // create a global numbering system
  iint *globalIds = (iint *) calloc((Nelements+mesh->totalHaloPairs)*Np,sizeof(iint));

  // every degree of freedom has its own global id
  MPI_Allgather(&(Nelements), 1, MPI_IINT, globalStarts+1, 1, MPI_IINT, MPI_COMM_WORLD);
    for(iint r=0;r<size;++r)
      globalStarts[r+1] = globalStarts[r]+globalStarts[r+1]*Np;

  /* so find number of elements on each rank */
  iint *rankNelements = (iint*) calloc(size, sizeof(iint));
  iint *rankStarts = (iint*) calloc(size+1, sizeof(iint));
  MPI_Allgather(&(Nelements), 1, MPI_IINT,
    rankNelements, 1, MPI_IINT, MPI_COMM_WORLD);
  //find offsets
  for(iint r=0;r<size;++r){
    rankStarts[r+1] = rankStarts[r]+rankNelements[r];
  }
  //use the offsets to set a global id
  for (iint e =0;e<Nelements;e++) {
    for (int n=0;n<Np;n++) {
      globalIds[e*Np +n] = n + (e + rankStarts[rankM])*Np;
    }
  }

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    iint *idSendBuffer = (iint *) calloc(Np*mesh->totalHaloPairs,sizeof(iint));
    meshHaloExchange(mesh, Np*sizeof(iint), globalIds, idSendBuffer, globalIds + Nelements*Np);
    free(idSendBuffer);
  }

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(Nfaces*Nfp*Nfp,sizeof(dfloat));
  for (iint f=0;f<Nfaces;f++) {
    for (iint n=0;n<Nfp;n++) {
      iint fn = mesh->faceNodes[f*Nfp+n];

      for (iint m=0;m<Nfp;m++) {
        dfloat MSnm = 0;

        for (iint i=0;i<Np;i++){
          MSnm += mesh->MM[fn+i*Np]*mesh->LIFT[i*Nfp*Nfaces+f*Nfp+m];
        }
        MS[m+n*Nfp + f*Nfp*Nfp]  = MSnm;
      }
    }
  }

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-10;

  /* Construct gradient as block matrix and load it to the halo */
  int GblockSize = Np*Np*(Nfaces+1);
  iint gradNNZ = GblockSize*(Nelements+mesh->totalHaloPairs);
  iint *Gids = (iint *) calloc(Np*(Nfaces+1)*(Nelements+mesh->totalHaloPairs),sizeof(iint));
  dfloat  *Gx = (dfloat *) calloc(gradNNZ,sizeof(dfloat));
  dfloat  *Gy = (dfloat *) calloc(gradNNZ,sizeof(dfloat));

  for (iint eM=0;eM<Nelements;eM++) {

    iint blockId = eM*GblockSize;

    iint vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat J = mesh->vgeo[vbase+JID];

    for(iint n=0;n<Np;++n){
      for(iint m=0;m<Np;++m){
        Gx[m+n*Np+blockId] = drdx*mesh->Dr[m+n*Np]+dsdx*mesh->Ds[m+n*Np];
        Gy[m+n*Np+blockId] = drdy*mesh->Dr[m+n*Np]+dsdy*mesh->Ds[m+n*Np];
      }
      Gids[n+eM*(Nfaces+1)*Np] = globalIds[eM*Np+n];
    }

    for (iint fM=0;fM<Nfaces;fM++) {
      // load surface geofactors for this face
      iint sid = mesh->Nsgeo*(eM*Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat invJ = mesh->sgeo[sid+IJID];

      iint eP = mesh->EToE[eM*Nfaces+fM];
      if (eP < 0) eP = eM;

      int bcD = 0, bcN =0;
      int bc = mesh->EToB[fM+Nfaces*eM]; //raw boundary flag
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

      // lift term
      for(iint n=0;n<Np;++n){
        for(iint m=0;m<Nfp;++m){
          iint idM = eM*Nfp*Nfaces+fM*Nfp+m;
          //iint nM = mesh->faceNodes[fM*Nfp+n];
          iint mM = mesh->faceNodes[fM*Nfp+m];
          iint mP = mesh->vmapP[idM]%Np;

          dfloat LIFTfnm = sJ*invJ*mesh->LIFT[m + fM*Nfp + n*Nfp*Nfaces];

          // G = sJ/J*LIFT*n*[[ uP-uM ]]
          Gx[mM+n*Np+             blockId] += -0.5*(1-bcN)*(1+bcD)*nx*LIFTfnm;
          Gy[mM+n*Np+             blockId] += -0.5*(1-bcN)*(1+bcD)*ny*LIFTfnm;

          Gx[mP+n*Np+(fM+1)*Np*Np+blockId] +=  0.5*(1-bcN)*(1-bcD)*nx*LIFTfnm;
          Gy[mP+n*Np+(fM+1)*Np*Np+blockId] +=  0.5*(1-bcN)*(1-bcD)*ny*LIFTfnm;
        }
        Gids[n+(fM+1)*Np+eM*(Nfaces+1)*Np] = globalIds[eP*Np+n];
      }
    }
  }

  //halo exchange the grad operators
  dfloat *GSendBuffer = (dfloat *) calloc(GblockSize*mesh->totalHaloPairs,sizeof(dfloat));
  iint   *GidsSendBuffer  = (iint *)   calloc(Np*Nfaces*mesh->totalHaloPairs,sizeof(iint));
  meshHaloExchange(mesh, GblockSize*sizeof(dfloat), Gx, GSendBuffer, Gx + Nelements*GblockSize);
  meshHaloExchange(mesh, GblockSize*sizeof(dfloat), Gy, GSendBuffer, Gy + Nelements*GblockSize);
  meshHaloExchange(mesh, Np*(Nfaces+1)*sizeof(iint), Gids, GidsSendBuffer, Gids + Np*(Nfaces+1)*Nelements);


  int AblockSize = Np*Np*(Nfaces+1)*(Nfaces+1);
  iint nnzLocalBound = AblockSize*Nelements;
  *A = (nonZero_t*) calloc(nnzLocalBound,sizeof(nonZero_t));

  dfloat *Ae = (dfloat*) calloc(AblockSize,sizeof(dfloat));

  // reset non-zero counter
  iint nnz = 0;

  // loop over all elements
  for(iint eM=0;eM<Nelements;++eM){

    iint vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat J = mesh->vgeo[vbase+JID];

    //zero out Ae
    for (int n=0;n<AblockSize;n++) Ae[n] = 0;

    /* start with stiffness matrix  */
    for(int n=0;n<Np;++n){
      for(int m=0;m<Np;++m){
        Ae[n*Np+m]  = J*lambda*mesh->MM[n*Np+m];
        Ae[n*Np+m] += J*drdx*drdx*mesh->Srr[n*Np+m];
        Ae[n*Np+m] += J*drdx*dsdx*mesh->Srs[n*Np+m];
        Ae[n*Np+m] += J*dsdx*drdx*mesh->Ssr[n*Np+m];
        Ae[n*Np+m] += J*dsdx*dsdx*mesh->Sss[n*Np+m];

        Ae[n*Np+m] += J*drdy*drdy*mesh->Srr[n*Np+m];
        Ae[n*Np+m] += J*drdy*dsdy*mesh->Srs[n*Np+m];
        Ae[n*Np+m] += J*dsdy*drdy*mesh->Ssr[n*Np+m];
        Ae[n*Np+m] += J*dsdy*dsdy*mesh->Sss[n*Np+m];
      }
    }

    for (iint fM=0;fM<Nfaces;fM++) {

      dfloat *AeP = Ae + (fM+1)*Np*Np*(Nfaces+1);

      // load surface geofactors for this face
      iint sid = mesh->Nsgeo*(eM*Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat hinv = mesh->sgeo[sid+IHID];

      iint eP = mesh->EToE[eM*Nfaces+fM];
      if (eP < 0) eP = eM;

      int bcD = 0, bcN =0;
      int bc = mesh->EToB[fM+Nfaces*eM]; //raw boundary flag
      int bcType = 0;

      if(bc>0) bcType = BCType[bc];          //find its type (Dirichlet/Neumann)

      // this needs to be double checked (and the code where these are used)
      if(bcType==1){ // Dirichlet
        bcD = 1;
        bcN = 0;
      } else if(bcType==2){ // Neumann
        bcD = 0;
        bcN = 1;
      }

      // mass matrix for this face
      dfloat *MSf = MS + fM*Nfp*Nfp;

      // penalty term just involves face nodes
      for(iint n=0;n<Nfp;++n){
        for(iint m=0;m<Nfp;++m){
          iint idM = eM*Nfp*Nfaces+fM*Nfp+m;
          int nM = mesh->faceNodes[fM*Nfp+n];
          int mM = mesh->faceNodes[fM*Nfp+m];
          int mP = mesh->vmapP[idM]%Np;

          dfloat MSfnm = sJ*MSf[n*Nfp+m];
          Ae [nM*Np+mM] +=  0.5*(1.-bcN)*(1.+bcD)*tau*MSfnm;
          AeP[nM*Np+mP] += -0.5*(1.-bcN)*(1.-bcD)*tau*MSfnm;
        }
      }

      // now add differential surface terms
      for(iint n=0;n<Nfp;++n){
        for(iint m=0;m<Np;++m){
          int nM = mesh->faceNodes[fM*Nfp+n];

          for(iint i=0;i<Nfp;++i){
            int iM = mesh->faceNodes[fM*Nfp+i];
            int iP = mesh->vmapP[i+fM*Nfp+eM*Nfp*Nfaces]%Np;

            dfloat MSfni = sJ*MSf[n*Nfp+i]; // surface Jacobian built in

            for (int fP=0;fP<Nfaces+1;fP++) {
              dfloat DxMim = Gx[m+iM*Np+fP*Np*Np+eM*GblockSize];
              dfloat DyMim = Gy[m+iM*Np+fP*Np*Np+eM*GblockSize];

              dfloat DxPim = Gx[m+iP*Np+fP*Np*Np+eP*GblockSize];
              dfloat DyPim = Gy[m+iP*Np+fP*Np*Np+eP*GblockSize];

              Ae [m+nM*Np+fP*Np*Np] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
              Ae [m+nM*Np+fP*Np*Np] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;

              AeP[m+nM*Np+fP*Np*Np] += -0.5*nx*(1-bcD)*(1-bcN)*MSfni*DxPim;
              AeP[m+nM*Np+fP*Np*Np] += -0.5*ny*(1-bcD)*(1-bcN)*MSfni*DyPim;
            }
          }
        }
      }

      for(iint n=0;n<Np;++n){
        for(iint m=0;m<Nfp;++m){
          int mM = mesh->faceNodes[fM*Nfp+m];
          int mP = mesh->vmapP[m + fM*Nfp+eM*Nfp*Nfaces]%Np;

          for(iint i=0;i<Nfp;++i){
            int iM = mesh->faceNodes[fM*Nfp+i];

            dfloat MSfim = sJ*MSf[i*Nfp+m];

            dfloat DxMin = drdx*mesh->Dr[iM*Np+n] + dsdx*mesh->Ds[iM*Np+n];
            dfloat DyMin = drdy*mesh->Dr[iM*Np+n] + dsdy*mesh->Ds[iM*Np+n];

            Ae [mM+n*Np] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
            Ae [mM+n*Np] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;

            AeP[mP+n*Np] +=  +0.5*nx*(1-bcD)*(1-bcN)*DxMin*MSfim;
            AeP[mP+n*Np] +=  +0.5*ny*(1-bcD)*(1-bcN)*DyMin*MSfim;
          }
        }
      }
    }

    for (int fM=0;fM<Nfaces+1;fM++) {
      iint eP = (fM==0) ? eM : mesh->EToE[eM*Nfaces+fM-1];
      if (eP<0) continue;

      for (int fP=0;fP<Nfaces+1;fP++) {
        for (int n=0;n<Np;n++) {
          for (int m=0;m<Np;m++) {
            dfloat Anm = Ae[m+n*Np+fP*Np*Np+fM*(Nfaces+1)*Np*Np];

            if(fabs(Anm)>tol){
              (*A)[nnz].row = globalIds[eM*Np+n];
              (*A)[nnz].col = Gids[m+fP*Np+eP*(Nfaces+1)*Np]; //use the column ids from the gradient lift

              (*A)[nnz].val = Anm;
              (*A)[nnz].ownerRank = rankM;
              ++nnz;
            }
          }
        }
      }
    }
  }

  // sort non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), nnz, sizeof(nonZero_t), parallelCompareRowColumn);

  // compress duplicates
  iint cnt = 0;
  for(iint n=1;n<nnz;++n){
    if((*A)[n].row == (*A)[cnt].row &&
       (*A)[n].col == (*A)[cnt].col){
      (*A)[cnt].val += (*A)[n].val;
    }
    else{
      ++cnt;
      (*A)[cnt] = (*A)[n];
    }
  }
  *nnzA = cnt+1;

  printf("nnz = %d\n", *nnzA);

  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));

  free(globalIds);

  free(Gx); free(Gy);
  free(Gids);


}
