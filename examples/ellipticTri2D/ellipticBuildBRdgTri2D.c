
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


  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-10;

  dfloat *qmP = (dfloat *) calloc(Nfp,sizeof(dfloat));
  dfloat *qmM = (dfloat *) calloc(Nfp,sizeof(dfloat));
  dfloat *QmP = (dfloat *) calloc(Nfp*(Nfaces+1)*Np,sizeof(dfloat));
  dfloat *QmM = (dfloat *) calloc(Nfp*(Nfaces+1)*Np,sizeof(dfloat));
  
  /* Construct gradient as block matrix and load it to the halo */
  //storage for gradient operator
  iint gradNNZ = Np*Np*(1+Nfaces)*(Nelements+mesh->totalHaloPairs);
  nonZero_t *Gx = (nonZero_t *) calloc(gradNNZ,sizeof(nonZero_t));
  nonZero_t *Gy = (nonZero_t *) calloc(gradNNZ,sizeof(nonZero_t));

  for (iint eM=0;eM<Nelements;eM++) {
    iint id = eM*Np*Np*(Nfaces+1);

    iint vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];

    for(iint n=0;n<Np;++n){
      for(iint m=0;m<Np;++m){
        Gx[n+m*Np+id].row = globalIds[eM*Np + m];
        Gx[n+m*Np+id].col = globalIds[eM*Np + n];
        Gx[n+m*Np+id].ownerRank = rankM;
        Gy[n+m*Np+id].row = globalIds[eM*Np + m];
        Gy[n+m*Np+id].col = globalIds[eM*Np + n];
        Gy[n+m*Np+id].ownerRank = rankM;

        Gx[n+m*Np+id].val = drdx*mesh->Dr[n+m*Np]+dsdx*mesh->Ds[n+m*Np];
        Gy[n+m*Np+id].val = drdy*mesh->Dr[n+m*Np]+dsdy*mesh->Ds[n+m*Np];
      }
    }

    for (iint m=0;m<Np;m++) {
      for (iint fM=0;fM<Nfaces;fM++) {

        id = (fM+1)*Np*Np + eM*Np*Np*(Nfaces+1);

        // load surface geofactors for this face
        iint sid = mesh->Nsgeo*(eM*Nfaces+fM);
        dfloat nx = mesh->sgeo[sid+NXID];
        dfloat ny = mesh->sgeo[sid+NYID];
        dfloat sJ = mesh->sgeo[sid+SJID];

        iint eP = mesh->EToE[eM*Nfaces+fM];
      
        // extract trace nodes
        for (iint i=0;i<Nfp;i++) {
          // double check vol geometric factors are in halo storage of vgeo
          iint idM    = eM*Nfp*Nfaces+fM*Nfp+i;
          iint vidM   = mesh->faceNodes[i+fM*Nfp];
          iint vidP   = mesh->vmapP[idM]%Np; // only use this to identify location of positive trace vgeo

          qmM[i] =0;
          if (vidM == m) qmM[i] =1;
          qmP[i] =0;
          if (vidP == m) qmP[i] =1;
        }

        id = Np*Np + eM*Np*Np*(Nfaces+1);
        iint fid = (fM+1)*Np*Np + eM*Np*Np*(Nfaces+1);

        int bcD = 0;
        int bcN = 0;
        if (eP < 0) {
          int bc = mesh->EToB[fM+Nfaces*eM]; //raw boundary flag
          iint bcType = BCType[bc];          //find its type (Dirichlet/Neumann)
          if(bcType==1){ // Dirichlet
            bcD = 1;
            bcN = 0;
          } else if (bcType==2){ // Neumann
            bcD = 0;
            bcN = 1;
          } else { // Neumann for now
            bcD = 0;
            bcN = 1;
          }
          eP=eM;
        }

        for (iint n=0;n<Np;n++) {
          Gx[m+n*Np+fid].row = globalIds[eM*Np + n];
          Gx[m+n*Np+fid].col = globalIds[eP*Np + m];
          
          for (iint i=0;i<Nfp;i++) {
            Gx[m+n*Np+id].val += -0.5*(1-bcN)*(1+bcD)*sJ*nx*mesh->LIFT[i+n*Nfp+fM*Nfp*Np]*qmM[i];
            Gy[m+n*Np+id].val += -0.5*(1-bcN)*(1+bcD)*sJ*ny*mesh->LIFT[i+n*Nfp+fM*Nfp*Np]*qmM[i];
            
            Gx[m+n*Np+fid].val += 0.5*(1-bcN)*(1-bcD)*sJ*nx*mesh->LIFT[i+n*Nfp+fM*Nfp*Np]*qmP[i];
            Gy[m+n*Np+fid].val += 0.5*(1-bcN)*(1-bcD)*sJ*ny*mesh->LIFT[i+n*Nfp+fM*Nfp*Np]*qmP[i];           
          }
        }
      }
    }
  }

  //halo exchange the grad operator
  int GblockSize = Np*Np*(Nfaces+1);
  nonZero_t *GSendBuffer = (nonZero_t *) calloc(GblockSize*mesh->totalHaloPairs,sizeof(nonZero_t));
  meshHaloExchange(mesh, GblockSize*sizeof(nonZero_t), Gx, GSendBuffer, Gx + Nelements*GblockSize);
  meshHaloExchange(mesh, GblockSize*sizeof(nonZero_t), Gy, GSendBuffer, Gy + Nelements*GblockSize);


  int AblockSize = Np*Np*(Nfaces+1)*(Nfaces+1);
  iint nnzLocalBound = AblockSize*Nelements;
  *A = (nonZero_t*) calloc(nnzLocalBound,sizeof(nonZero_t));

  dfloat *Ae = (dfloat*) calloc(AblockSize,sizeof(dfloat));

  // reset non-zero counter
  iint nnz = 0;

  // loop over all elements
  for(iint eM=0;eM<Nelements;++eM){

    //zero out Ae
    for (int n=0;n<AblockSize;n++) Ae[n] = 0.;

    iint vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat J = mesh->vgeo[vbase+JID];

    for(int n=0;n<Np;++n){
      for(int m=0;m<Np;++m){
        for (int fP = 0; fP<Nfaces+1;fP++) {
          iint id = fP*Np*Np + eM*GblockSize;          
          for (int k=0;k<Np;k++) {
            //Ae[n+m*Np+fP*(Nfaces+1)*Np*Np] += (drdx*mesh->Dr[k+m*Np]+dsdx*mesh->Ds[k+m*Np])*Gx[n+k*Np+id].val;
            //Ae[n+m*Np+fP*(Nfaces+1)*Np*Np] += (drdy*mesh->Dr[k+m*Np]+dsdy*mesh->Ds[k+m*Np])*Gy[n+k*Np+id].val;
          }  
        }
      }
      Ae[n+n*Np] -= lambda;  
    }

    for (iint m=0;m<Np;m++) {
      for (iint fM=0;fM<Nfaces;fM++) {
        // load surface geofactors for this face
        iint sid = mesh->Nsgeo*(eM*Nfaces+fM);
        dfloat nx = mesh->sgeo[sid+NXID];
        dfloat ny = mesh->sgeo[sid+NYID];
        dfloat sJ = mesh->sgeo[sid+SJID];

        iint eP = mesh->EToE[eM*Nfaces+fM];
        if (eP<0) eP = eM;
      
        // extract trace matrix from Gx and Gy
        for (iint i=0;i<Nfp;i++) {
          // double check vol geometric factors are in halo storage of vgeo
          iint idM    = eM*Nfp*Nfaces+fM*Nfp+i;
          iint vidM   = mesh->faceNodes[i+fM*Nfp];
          iint vidP   = mesh->vmapP[idM]%Np; // only use this to identify location of positive trace vgeo

          qmM[i] = 0;          
          for (int fP=0;fP<Nfaces+1;fP++) {
            for (int j=0;j<Np;j++) 
              QmM[j+fP*Np+i*(Nfaces+1)*Np] = 0;
          }
          if (vidM == m) {
            qmM[i] = 1;
            for (int fP=0;fP<Nfaces+1;fP++) {
              iint id = fP*Np*Np + eM*Np*Np*(Nfaces+1);
              for (int j=0;j<Np;j++) 
                QmM[j+fP*Np+i*Np] = nx*Gx[j+m*Np+id].val+ny*Gy[j+m*Np+id].val;              
            }
          }
          
          qmP[i] = 0;          
          for (int fP=0;fP<Nfaces+1;fP++) {
            for (int j=0;j< Np;j++) 
              QmP[j+fP*Np+i*(Nfaces+1)*Np] = 0;
          }
          if (vidP == m) {
            qmP[i] = 1;
            for (int fP=0;fP<Nfaces+1;fP++) {
              iint id = fP*Np*Np + eP*Np*Np*(Nfaces+1);
              for (int j=0;j<Np;j++) 
                QmP[j+fP*Np+i*Np] = nx*Gx[j+m*Np+id].val + ny*Gy[j+m*Np+id].val;
            }
          }
        }

        int bcD = 0;
        int bcN = 0;
        eP = mesh->EToE[eM*Nfaces+fM];
        if (eP < 0) {
          int bc = mesh->EToB[fM+Nfaces*eM]; //raw boundary flag
          iint bcType = BCType[bc];          //find its type (Dirichlet/Neumann)
          if(bcType==1){ // Dirichlet
            bcD = 1;
            bcN = 0;
          } else if (bcType==2){ // Neumann
            bcD = 0;
            bcN = 1;
          } else { // Neumann for now
            bcD = 0;
            bcN = 1;
          }
          eP=eM;
        }

        for (int n=0;n<Np;n++) {
          for (int fP=0;fP<Nfaces;fP++) {
            for (iint i=0;i<Nfp;i++) {
              //Ae[m+n*Np+             fP*(Nfaces+1)*Np*Np] += -0.5*(1+bcN)*(1-bcD)*sJ*mesh->LIFT[i+n*Nfp+fM*Nfp*Np]*QmM[n+fP*Np+i*Np];
              //Ae[m+n*Np+(fM+1)*Np*Np+fP*(Nfaces+1)*Np*Np] +=  0.5*(1-bcN)*(1-bcD)*sJ*mesh->LIFT[i+n*Nfp+fM*Nfp*Np]*QmP[n+fP*Np+i*Np];
            }
          }
          for (iint i=0;i<Nfp;i++) {
            Ae[m+n*Np             ] += -0.5*(1-bcN)*(1+bcD)*sJ*mesh->LIFT[i+n*Nfp+fM*Nfp*Np]*qmM[i];
            Ae[m+n*Np+(fM+1)*Np*Np] +=  0.5*(1-bcN)*(1-bcD)*sJ*mesh->LIFT[i+n*Nfp+fM*Nfp*Np]*qmP[i];
          }
        }
      }
    }
    
    //multiply by mass matrix 
    for (int fM=0;fM<Nfaces+1;fM++) {
      iint eP = (fM==0) ? eM : mesh->EToE[eM*Nfaces+fM];
      if (eP<0) eP = eM;

      for (int fP=0;fP<Nfaces+1;fP++) {
        for (int n=0;n<Np;n++) {
          for (int m=0;m<Np;m++) {
            dfloat Anm = 0.;
            for (int k=0;k<Np;k++) {
              Anm += mesh->MM[k+n*Np]*Ae[m+k*Np+fM*Np*Np+fP*(Nfaces+1)*Np*Np]; 
            }
            Anm *= J;
            
            if(fabs(Anm)>tol){
              (*A)[nnz].row = globalIds[eM*Np+n];
              (*A)[nnz].col = Gx[m+fP*Np*Np+eP*GblockSize].col;
              (*A)[nnz].val = -Anm;
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


  printf("BR matrix \n");
  dfloat *Ap = (dfloat*) calloc(Np*Nelements*Np*Nelements,sizeof(dfloat));
  for (int n=0;n<*nnzA;n++){
    int row = (*A)[n].row;
    int col = (*A)[n].col;
    Ap[col+row*mesh->Np*Nelements] = (*A)[n].val;
  }

  for (int n=0;n<Np*Nelements;n++) {
    for (int m=0;m<Np*Nelements;m++) {
      printf("%4.2f \t", Ap[m + n*Np*Nelements]);
    }
    printf("\n");
  }

  printf("nnz = %d\n", *nnzA);
  
  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));

  free(globalIds);

  free(Gx); free(Gy); free(Ae);

  free(qmM); free(qmP);
  free(QmM); free(QmP);
}
