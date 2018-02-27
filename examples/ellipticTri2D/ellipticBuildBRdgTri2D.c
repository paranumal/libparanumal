
#include "ellipticTri2D.h"

int parallelCompareRowColumn(const void *a, const void *b);

void ellipticBuildBRdgTri2D(mesh2D *mesh, int basisNp, dfloat *basis,
                            dfloat tau, dfloat lambda, int *BCType, nonZero_t **A,
                            long long int *nnzA, hlong *globalStarts, const char *options){

  int size, rankM;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);

  int Np = mesh->Np;
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;
  dlong Nelements = mesh->Nelements;

  if(!basis) { // default to degree N Lagrange basis
    basisNp = Np;
    basis = (dfloat*) calloc(basisNp*basisNp, sizeof(dfloat));
    for(int n=0;n<basisNp;++n){
      basis[n+n*basisNp] = 1;
    }
  }

  hlong Nnum = Np*Nelements;

  // create a global element numbering system
  hlong *globalIds = (hlong *) calloc(Nelements+mesh->totalHaloPairs,sizeof(hlong));

  // every degree of freedom has its own global id
  MPI_Allgather(&Nnum, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, MPI_COMM_WORLD);
  for(int r=0;r<size;++r)
    globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

  /* so find number of elements on each rank */
  dlong *rankNelements = (dlong*) calloc(size, sizeof(dlong));
  hlong *rankStarts = (hlong*) calloc(size+1, sizeof(hlong));
  MPI_Allgather(&Nelements, 1, MPI_DLONG, rankNelements, 1, MPI_DLONG, MPI_COMM_WORLD);
  //find offsets
  for(int r=0;r<size;++r)
    rankStarts[r+1] = rankStarts[r]+rankNelements[r];

  //use the offsets to set a global id
  for (dlong e =0;e<Nelements;e++)
    globalIds[e] = e + rankStarts[rankM];

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    hlong *idSendBuffer = (hlong *) calloc(mesh->totalHaloPairs,sizeof(hlong));
    meshHaloExchange(mesh, sizeof(hlong), globalIds, idSendBuffer, globalIds + Nelements);
    free(idSendBuffer);
  }

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(Nfaces*Nfp*Nfp,sizeof(dfloat));
  for (int f=0;f<Nfaces;f++) {
    for (int n=0;n<Nfp;n++) {
      int fn = mesh->faceNodes[f*Nfp+n];

      for (int m=0;m<Nfp;m++) {
        dfloat MSnm = 0;

        for (int i=0;i<Np;i++){
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
  long long int gradNNZ = GblockSize*(Nelements+mesh->totalHaloPairs);
  hlong *Gids = (hlong *) calloc((Nfaces+1)*(Nelements+mesh->totalHaloPairs),sizeof(hlong));
  dfloat  *Gx = (dfloat *) calloc(gradNNZ,sizeof(dfloat));
  dfloat  *Gy = (dfloat *) calloc(gradNNZ,sizeof(dfloat));

  for (dlong eM=0;eM<Nelements;eM++) {

    long long int blockId = eM*GblockSize;

    dlong vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat J = mesh->vgeo[vbase+JID];

    Gids[eM*(Nfaces+1)] = globalIds[eM]; //record global element number
    for(int n=0;n<Np;++n){
      for(int m=0;m<Np;++m){
        Gx[m+n*Np+blockId] = drdx*mesh->Dr[m+n*Np]+dsdx*mesh->Ds[m+n*Np];
        Gy[m+n*Np+blockId] = drdy*mesh->Dr[m+n*Np]+dsdy*mesh->Ds[m+n*Np];
      }
    }

    for (int fM=0;fM<Nfaces;fM++) {
      // load surface geofactors for this face
      dlong sid = mesh->Nsgeo*(eM*Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat invJ = mesh->sgeo[sid+IJID];

      dlong eP = mesh->EToE[eM*Nfaces+fM];
      if (eP < 0) eP = eM;

      Gids[fM+1 + eM*(Nfaces+1)] = globalIds[eP]; //record global element number

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

      // lift term
      for(int n=0;n<Np;++n){
        for(int m=0;m<Nfp;++m){
          dlong idM = eM*Nfp*Nfaces+fM*Nfp+m;
          //int nM = mesh->faceNodes[fM*Nfp+n];
          int mM = mesh->faceNodes[fM*Nfp+m];
          int mP = (int) (mesh->vmapP[idM]%Np);

          dfloat LIFTfnm = sJ*invJ*mesh->LIFT[m + fM*Nfp + n*Nfp*Nfaces];

          // G = sJ/J*LIFT*n*[[ uP-uM ]]
          Gx[mM+n*Np+             blockId] += -0.5*(1-bcN)*(1+bcD)*nx*LIFTfnm;
          Gy[mM+n*Np+             blockId] += -0.5*(1-bcN)*(1+bcD)*ny*LIFTfnm;

          Gx[mP+n*Np+(fM+1)*Np*Np+blockId] +=  0.5*(1-bcN)*(1-bcD)*nx*LIFTfnm;
          Gy[mP+n*Np+(fM+1)*Np*Np+blockId] +=  0.5*(1-bcN)*(1-bcD)*ny*LIFTfnm;
        }
      }
    }
  }

  //halo exchange the grad operators
  dfloat *GSendBuffer = (dfloat *) calloc(GblockSize*mesh->totalHaloPairs,sizeof(dfloat));
  hlong  *GidsSendBuffer  = (hlong *) calloc((Nfaces+1)*mesh->totalHaloPairs,sizeof(hlong));
  meshHaloExchange(mesh, GblockSize*sizeof(dfloat), Gx, GSendBuffer, Gx + Nelements*GblockSize);
  meshHaloExchange(mesh, GblockSize*sizeof(dfloat), Gy, GSendBuffer, Gy + Nelements*GblockSize);
  meshHaloExchange(mesh, (Nfaces+1)*sizeof(hlong), Gids, GidsSendBuffer, Gids + (Nfaces+1)*Nelements);

  hlong *globalPatchIds = (hlong*) calloc((Nfaces+1)*(Nfaces+1),sizeof(hlong));
  int *patchIds = (int*) calloc((Nfaces+1)*(Nfaces+1),sizeof(int));

  int AblockSize = Np*Np*(Nfaces*Nfaces+1);
  long long int nnzLocalBound = AblockSize*Nelements;
  *A = (nonZero_t*) calloc(nnzLocalBound,sizeof(nonZero_t));

  dfloat *Ae = (dfloat*) calloc(AblockSize,sizeof(dfloat));

  // reset non-zero counter
  long long int nnz = 0;

  // loop over all elements
  for(dlong eM=0;eM<Nelements;++eM){

    dlong vbase = eM*mesh->Nvgeo;
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

    //check and record all the global element numbers in this patch
    // and give them a local ordering so each element knows what patch to contribute to
    int NpatchElements=1;
    globalPatchIds[0] = globalIds[eM];
    for (int fM=0;fM<Nfaces+1;fM++) {
      dlong eP = (fM==0) ? eM : mesh->EToE[eM*Nfaces+fM-1];
      if (eP < 0) eP = eM;

      for (int fP=0;fP<Nfaces+1;fP++) {
        hlong gid = Gids[fP+eP*(Nfaces+1)];
        int i=0;
        for (;i<NpatchElements;i++) { //loop through the list of ids in the patch
          //check if we have this elements id already
          if (gid==globalPatchIds[i]) {
            patchIds[fP + fM*(Nfaces+1)] = i;
            break;
          }
        }
        if (i==NpatchElements) {
          globalPatchIds[NpatchElements] = gid;
          patchIds[fP + fM*(Nfaces+1)] = NpatchElements;
          NpatchElements++;
        }
      }
    }

    for (int fM=0;fM<Nfaces;fM++) {
      // load surface geofactors for this face
      dlong sid = mesh->Nsgeo*(eM*Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat hinv = mesh->sgeo[sid+IHID];

      dlong eP = mesh->EToE[eM*Nfaces+fM];
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

      //local patch id for this neighbor
      int pid = patchIds[fM+1];

      // penalty term just involves face nodes
      for(int n=0;n<Nfp;++n){
        for(int m=0;m<Nfp;++m){
          dlong idM = eM*Nfp*Nfaces+fM*Nfp+m;
          int nM = mesh->faceNodes[fM*Nfp+n];
          int mM = mesh->faceNodes[fM*Nfp+m];
          int mP = (int) (mesh->vmapP[idM]%Np);

          dfloat MSfnm = sJ*MSf[n*Nfp+m];
          Ae[nM*Np+mM          ] +=  0.5*(1.-bcN)*(1.+bcD)*tau*MSfnm;
          Ae[nM*Np+mP+pid*Np*Np] += -0.5*(1.-bcN)*(1.-bcD)*tau*MSfnm;
        }
      }

      // now add differential surface terms
      for(int n=0;n<Nfp;++n){
        for(int m=0;m<Np;++m){
          int nM = mesh->faceNodes[fM*Nfp+n];

          for(int i=0;i<Nfp;++i){
            int iM = mesh->faceNodes[fM*Nfp+i];
            int iP = (int) (mesh->vmapP[i+fM*Nfp+eM*Nfp*Nfaces]%Np);

            dfloat MSfni = sJ*MSf[n*Nfp+i]; // surface Jacobian built in

            for (int fP=0;fP<Nfaces+1;fP++) {
              //local patch id for this neighbor
              int pidM = patchIds[fP];
              int pidP = patchIds[fP+(fM+1)*(Nfaces+1)];

              dfloat DxMim = Gx[m+iM*Np+fP*Np*Np+eM*GblockSize];
              dfloat DyMim = Gy[m+iM*Np+fP*Np*Np+eM*GblockSize];

              dfloat DxPim = Gx[m+iP*Np+fP*Np*Np+eP*GblockSize];
              dfloat DyPim = Gy[m+iP*Np+fP*Np*Np+eP*GblockSize];

              Ae[m+nM*Np+pidM*Np*Np] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
              Ae[m+nM*Np+pidM*Np*Np] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;

              Ae[m+nM*Np+pidP*Np*Np] += -0.5*nx*(1-bcD)*(1-bcN)*MSfni*DxPim;
              Ae[m+nM*Np+pidP*Np*Np] += -0.5*ny*(1-bcD)*(1-bcN)*MSfni*DyPim;
            }
          }
        }
      }

      for(int n=0;n<Np;++n){
        for(int m=0;m<Nfp;++m){
          int mM = mesh->faceNodes[fM*Nfp+m];
          int mP = (int) (mesh->vmapP[m + fM*Nfp+eM*Nfp*Nfaces]%Np);

          for(int i=0;i<Nfp;++i){
            int iM = mesh->faceNodes[fM*Nfp+i];

            dfloat MSfim = sJ*MSf[i*Nfp+m];

            dfloat DxMin = drdx*mesh->Dr[iM*Np+n] + dsdx*mesh->Ds[iM*Np+n];
            dfloat DyMin = drdy*mesh->Dr[iM*Np+n] + dsdy*mesh->Ds[iM*Np+n];

            Ae[mM+n*Np          ] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
            Ae[mM+n*Np          ] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;

            Ae[mP+n*Np+pid*Np*Np] +=  +0.5*nx*(1-bcD)*(1-bcN)*DxMin*MSfim;
            Ae[mP+n*Np+pid*Np*Np] +=  +0.5*ny*(1-bcD)*(1-bcN)*DyMin*MSfim;
          }
        }
      }
    }

    for (int e=0;e<NpatchElements;e++) {
      dlong eP = globalPatchIds[e];

      for(int j=0;j<basisNp;++j){
        for(int i=0;i<basisNp;++i){
          dfloat val = 0;
          for (int n=0;n<Np;n++) {
            for (int m=0;m<Np;m++) {
              val += basis[n*Np+j]*Ae[n*Np+m+e*Np*Np]*basis[m*Np+i];
            }
          }

          if(fabs(val)>tol){
            (*A)[nnz].row = eM*Np+j;
            (*A)[nnz].col = eP*Np+i;

            (*A)[nnz].val = val;
            (*A)[nnz].ownerRank = rankM;
            ++nnz;
          }
        }
      }
    }
  }

  // sort non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), nnz, sizeof(nonZero_t), parallelCompareRowColumn);

  *nnzA = nnz;

#if 0
  dfloat* Ap = (dfloat *) calloc(Np*Np*Nelements*Nelements,sizeof(dfloat));
  for (int n=0;n<nnz;n++) {
    int row = (*A)[n].row;
    int col = (*A)[n].col;

    Ap[col+row*Np*Nelements] = (*A)[n].val;
  }

  for (int i=0;i<Np*Nelements;i++) {
    for (int j =0;j<Nelements*Np;j++) {
      printf("%4.2f \t", Ap[j+i*Np*Nelements]);
    }
    printf("\n");
  }
#endif

  printf("nnz = %lld\n", *nnzA);

  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));

  free(globalIds);
  free(globalPatchIds); free(patchIds);

  free(Gx); free(Gy);
  free(Gids);


}
