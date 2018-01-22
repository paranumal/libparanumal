#include "ellipticTri2D.h"


// compare on global indices
int parallelCompareRowColumn(const void *a, const void *b){

  nonZero_t *fa = (nonZero_t*) a;
  nonZero_t *fb = (nonZero_t*) b;

  if(fa->row < fb->row) return -1;
  if(fa->row > fb->row) return +1;

  if(fa->col < fb->col) return -1;
  if(fa->col > fb->col) return +1;

  return 0;
}

void ellipticBuildContinuousTri2D(mesh2D *mesh, dfloat lambda, nonZero_t **A, iint *nnz, ogs_t **ogs, iint *globalStarts, const char* options) {

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Build a gather-scatter to assemble the global masked problem */
  iint Ntotal = mesh->Np*mesh->Nelements;

  iint *globalNumbering = (iint *) calloc(Ntotal,sizeof(iint));
  memcpy(globalNumbering,mesh->globalIds,Ntotal*sizeof(iint)); 
  for (iint n=0;n<mesh->Nmasked;n++) 
    globalNumbering[mesh->maskIds[n]] = -1;

  // squeeze node numbering
  meshParallelConsecutiveGlobalNumbering(mesh, Ntotal, globalNumbering, mesh->globalOwners, globalStarts);

  iint *gatherMaskedBaseIds   = (iint *) calloc(Ntotal,sizeof(iint));
  for (iint n=0;n<Ntotal;n++) {
    iint id = mesh->gatherLocalIds[n];
    gatherMaskedBaseIds[n] = globalNumbering[id];
  }

  //build gather scatter with masked nodes
  *ogs = meshParallelGatherScatterSetup(mesh, Ntotal, 
                                        mesh->gatherLocalIds,  gatherMaskedBaseIds, 
                                        mesh->gatherBaseRanks, mesh->gatherHaloFlags);

  // Build non-zeros of stiffness matrix (unassembled)
  iint nnzLocal = mesh->Np*mesh->Np*mesh->Nelements;

  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocal, sizeof(nonZero_t));
  iint *AsendCounts  = (iint*) calloc(size, sizeof(iint));
  iint *ArecvCounts  = (iint*) calloc(size, sizeof(iint));
  iint *AsendOffsets = (iint*) calloc(size+1, sizeof(iint));
  iint *ArecvOffsets = (iint*) calloc(size+1, sizeof(iint));

  dfloat *Srr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *Srs = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *Sss = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *MM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

  if (strstr(options,"SPARSE")) {
    for (int k=0;k<mesh->SparseNnzPerRow;k++) {
      for (int n=0;n<mesh->Np;n++) {
        int m = mesh->sparseStackedNZ[n+k*mesh->Np]-1;
        if (m>-1) {
          Srr[m+n*mesh->Np] = mesh->sparseSrrT[n+k*mesh->Np];
          Srs[m+n*mesh->Np] = mesh->sparseSrsT[n+k*mesh->Np];
          Sss[m+n*mesh->Np] = mesh->sparseSssT[n+k*mesh->Np];  
        }
      }
    }
    for (iint n=0;n<mesh->Np;n++) {
      for (iint m=0;m<mesh->Np;m++) {
        MM[m+n*mesh->Np] = mesh->sparseMM[m+n*mesh->Np];
      }
    }
  } else {
    for (iint n=0;n<mesh->Np;n++) {
      for (iint m=0;m<mesh->Np;m++) {
        Srr[m+n*mesh->Np] = mesh->Srr[m+n*mesh->Np];
        Srs[m+n*mesh->Np] = mesh->Srs[m+n*mesh->Np] + mesh->Ssr[m+n*mesh->Np];
        Sss[m+n*mesh->Np] = mesh->Sss[m+n*mesh->Np];
        MM[m+n*mesh->Np] = mesh->MM[m+n*mesh->Np];
      }
    }
  }

  int *mask = (int *) calloc(mesh->Np*mesh->Nelements,sizeof(int));
  for (iint n=0;n<mesh->Nmasked;n++) mask[mesh->maskIds[n]] = 1;

  if(rank==0) printf("Building full FEM matrix...");fflush(stdout);

  //Build unassembed non-zeros
  iint cnt =0;
  for (iint e=0;e<mesh->Nelements;e++) {
    for (iint n=0;n<mesh->Np;n++) {
      if (mask[e*mesh->Np + n]) continue; //skip masked nodes
      for (iint m=0;m<mesh->Np;m++) {
        if (mask[e*mesh->Np + m]) continue; //skip masked nodes

        dfloat val = 0.;

        dfloat Grr = mesh->ggeo[e*mesh->Nggeo + G00ID];
        dfloat Grs = mesh->ggeo[e*mesh->Nggeo + G01ID];
        dfloat Gss = mesh->ggeo[e*mesh->Nggeo + G11ID];
        dfloat J   = mesh->ggeo[e*mesh->Nggeo + GWJID];

        val += Grr*Srr[m+n*mesh->Np];
        val += Grs*Srs[m+n*mesh->Np];
        val += Gss*Sss[m+n*mesh->Np];
        val += J*lambda*MM[m+n*mesh->Np];

        if (strstr(options,"SPARSE")) 
          val *= mesh->mapSgn[m + e*mesh->Np]*mesh->mapSgn[n + e*mesh->Np];

        dfloat nonZeroThreshold = 1e-7;
        if (fabs(val)>nonZeroThreshold) {
          // pack non-zero
          sendNonZeros[cnt].val = val;
          sendNonZeros[cnt].row = globalNumbering[e*mesh->Np + n];
          sendNonZeros[cnt].col = globalNumbering[e*mesh->Np + m];
          sendNonZeros[cnt].ownerRank = mesh->globalOwners[e*mesh->Np + n];
          cnt++;
        }
      }
    }
  }

  // count how many non-zeros to send to each process
  for(iint n=0;n<cnt;++n)
    AsendCounts[sendNonZeros[n].ownerRank] += sizeof(nonZero_t);

  // sort by row ordering
  qsort(sendNonZeros, cnt, sizeof(nonZero_t), parallelCompareRowColumn);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(AsendCounts, 1, MPI_IINT, ArecvCounts, 1, MPI_IINT, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  *nnz = 0;
  for(iint r=0;r<size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    *nnz += ArecvCounts[r]/sizeof(nonZero_t);
  }

  *A = (nonZero_t*) calloc(*nnz, sizeof(nonZero_t));

  // determine number to receive
  MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_CHAR,
    (*A), ArecvCounts, ArecvOffsets, MPI_CHAR,
    MPI_COMM_WORLD);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), *nnz, sizeof(nonZero_t), parallelCompareRowColumn);

  // compress duplicates
  cnt = 0;
  for(iint n=1;n<*nnz;++n){
    if((*A)[n].row == (*A)[cnt].row &&
       (*A)[n].col == (*A)[cnt].col){
      (*A)[cnt].val += (*A)[n].val;
    }
    else{
      ++cnt;
      (*A)[cnt] = (*A)[n];
    }
  }
  *nnz = cnt+1;

  if(rank==0) printf("done.\n");

  free(sendNonZeros);
  free(globalNumbering);

  free(AsendCounts);
  free(ArecvCounts);
  free(AsendOffsets);
  free(ArecvOffsets);

  free(Srr);
  free(Srs);
  free(Sss);
  free(MM );
  free(mask);
}