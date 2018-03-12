#include "ellipticHex3D.h"

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

void ellipticBuildContinuousHex3D(solver_t *solver, dfloat lambda, nonZero_t **A, dlong *nnz, ogs_t **ogs, hlong *globalStarts, const char* options) {  
  
  mesh3D *mesh = solver->mesh;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Build a gather-scatter to assemble the global masked problem */
  dlong Ntotal = mesh->Np*mesh->Nelements;

  hlong *globalNumbering = (hlong *) calloc(Ntotal,sizeof(hlong));
  memcpy(globalNumbering,mesh->globalIds,Ntotal*sizeof(hlong)); 
  for (dlong n=0;n<solver->Nmasked;n++) 
    globalNumbering[solver->maskIds[n]] = -1;

  // squeeze node numbering
  meshParallelConsecutiveGlobalNumbering(mesh, Ntotal, globalNumbering, mesh->globalOwners, globalStarts);

  hlong *gatherMaskedBaseIds   = (hlong *) calloc(Ntotal,sizeof(hlong));
  for (dlong n=0;n<Ntotal;n++) {
    dlong id = mesh->gatherLocalIds[n];
    gatherMaskedBaseIds[n] = globalNumbering[id];
  }

  //build gather scatter with masked nodes
  int verbose = strstr(options,"VERBOSE") ? 1:0;
  *ogs = meshParallelGatherScatterSetup(mesh, Ntotal, 
                                        mesh->gatherLocalIds,  gatherMaskedBaseIds, 
                                        mesh->gatherBaseRanks, mesh->gatherHaloFlags,verbose);


    // 2. Build non-zeros of stiffness matrix (unassembled)
  dlong nnzLocal = mesh->Np*mesh->Np*mesh->Nelements;
  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocal, sizeof(nonZero_t));
  int *AsendCounts  = (int*) calloc(size, sizeof(int));
  int *ArecvCounts  = (int*) calloc(size, sizeof(int));
  int *AsendOffsets = (int*) calloc(size+1, sizeof(int));
  int *ArecvOffsets = (int*) calloc(size+1, sizeof(int));

  int *mask = (int *) calloc(mesh->Np*mesh->Nelements,sizeof(int));
  for (dlong n=0;n<solver->Nmasked;n++) mask[solver->maskIds[n]] = 1;

  if(rank==0) printf("Building full FEM matrix...");fflush(stdout);

  dlong cnt =0;
  for (dlong e=0;e<mesh->Nelements;e++) {
    for (int nz=0;nz<mesh->Nq;nz++) {
    for (int ny=0;ny<mesh->Nq;ny++) {
    for (int nx=0;nx<mesh->Nq;nx++) {
      int idn = nx+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
      if (mask[e*mesh->Np + idn]) continue; //skip masked nodes

        for (int mz=0;mz<mesh->Nq;mz++) {
        for (int my=0;my<mesh->Nq;my++) {
        for (int mx=0;mx<mesh->Nq;mx++) {
          int idm = mx+my*mesh->Nq+mz*mesh->Nq*mesh->Nq;
          if (mask[e*mesh->Np + idm]) continue; //skip masked nodes
          
            int id;
            dfloat val = 0.;
            
            if ((ny==my)&&(nz==mz)) {
              for (int k=0;k<mesh->Nq;k++) {
                id = k+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
                dfloat Grr = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G00ID*mesh->Np];

                val += Grr*mesh->D[nx+k*mesh->Nq]*mesh->D[mx+k*mesh->Nq];
              }
            }

            if (nz==mz) {
              id = mx+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
              dfloat Grs = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G01ID*mesh->Np];
              val += Grs*mesh->D[nx+mx*mesh->Nq]*mesh->D[my+ny*mesh->Nq];

              id = nx+my*mesh->Nq+nz*mesh->Nq*mesh->Nq;
              dfloat Gsr = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G01ID*mesh->Np];
              val += Gsr*mesh->D[mx+nx*mesh->Nq]*mesh->D[ny+my*mesh->Nq];
            }

            if (ny==my) {
              id = mx+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
              dfloat Grt = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G02ID*mesh->Np];
              val += Grt*mesh->D[nx+mx*mesh->Nq]*mesh->D[mz+nz*mesh->Nq];

              id = nx+ny*mesh->Nq+mz*mesh->Nq*mesh->Nq;
              dfloat Gst = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G02ID*mesh->Np];
              val += Gst*mesh->D[mx+nx*mesh->Nq]*mesh->D[nz+mz*mesh->Nq];
            }

            if ((nx==mx)&&(nz==mz)) {
              for (int k=0;k<mesh->Nq;k++) {
                id = nx+k*mesh->Nq+nz*mesh->Nq*mesh->Nq;
                dfloat Gss = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G11ID*mesh->Np];

                val += Gss*mesh->D[ny+k*mesh->Nq]*mesh->D[my+k*mesh->Nq];
              }
            }
            
            if (nx==mx) {
              id = nx+my*mesh->Nq+nz*mesh->Nq*mesh->Nq;
              dfloat Gst = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G12ID*mesh->Np];
              val += Gst*mesh->D[ny+my*mesh->Nq]*mesh->D[mz+nz*mesh->Nq];

              id = nx+ny*mesh->Nq+mz*mesh->Nq*mesh->Nq;
              dfloat Gts = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G12ID*mesh->Np];
              val += Gts*mesh->D[my+ny*mesh->Nq]*mesh->D[nz+mz*mesh->Nq];
            }

            if ((nx==mx)&&(ny==my)) {
              for (int k=0;k<mesh->Nq;k++) {
                id = nx+ny*mesh->Nq+k*mesh->Nq*mesh->Nq;
                dfloat Gtt = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G22ID*mesh->Np];

                val += Gtt*mesh->D[nz+k*mesh->Nq]*mesh->D[mz+k*mesh->Nq];
              }
            }
            
            if ((nx==mx)&&(ny==my)&&(nz==mz)) {
              id = nx + ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
              dfloat JW = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + GWJID*mesh->Np];
              val += JW*lambda;
            }
            
            // pack non-zero
            dfloat nonZeroThreshold = 1e-7;
            if (fabs(val) >= nonZeroThreshold) {
              sendNonZeros[cnt].val = val;
              sendNonZeros[cnt].row = globalNumbering[e*mesh->Np + idn];
              sendNonZeros[cnt].col = globalNumbering[e*mesh->Np + idm];
              sendNonZeros[cnt].ownerRank = mesh->globalOwners[e*mesh->Np + idn];
              cnt++;
            }
        }
        }
        }
      }
      }
      }
  }

  // Make the MPI_NONZERO_T data type
  MPI_Datatype MPI_NONZERO_T;
  MPI_Datatype dtype[4] = {MPI_HLONG, MPI_HLONG, MPI_INT, MPI_DFLOAT};
  int blength[4] = {1, 1, 1, 1};
  MPI_Aint addr[4], displ[4];
  MPI_Get_address ( &(sendNonZeros[0]          ), addr+0);
  MPI_Get_address ( &(sendNonZeros[0].col      ), addr+1);
  MPI_Get_address ( &(sendNonZeros[0].ownerRank), addr+2);
  MPI_Get_address ( &(sendNonZeros[0].val      ), addr+3);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  displ[3] = addr[3] - addr[0];
  MPI_Type_create_struct (4, blength, displ, dtype, &MPI_NONZERO_T);
  MPI_Type_commit (&MPI_NONZERO_T);

  // count how many non-zeros to send to each process
  for(dlong n=0;n<cnt;++n)
    AsendCounts[sendNonZeros[n].ownerRank]++;

  // sort by row ordering
  qsort(sendNonZeros, cnt, sizeof(nonZero_t), parallelCompareRowColumn);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(AsendCounts, 1, MPI_INT, ArecvCounts, 1, MPI_INT, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  *nnz = 0;
  for(int r=0;r<size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    *nnz += ArecvCounts[r];
  }

  *A = (nonZero_t*) calloc(*nnz, sizeof(nonZero_t));

  // determine number to receive
  MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_NONZERO_T,
                        (*A), ArecvCounts, ArecvOffsets, MPI_NONZERO_T,
                        MPI_COMM_WORLD);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), *nnz, sizeof(nonZero_t), parallelCompareRowColumn);

  // compress duplicates
  cnt = 0;
  for(dlong n=1;n<*nnz;++n){
    if((*A)[n].row == (*A)[cnt].row &&
       (*A)[n].col == (*A)[cnt].col){
      (*A)[cnt].val += (*A)[n].val;
    }
    else{
      ++cnt;
      (*A)[cnt] = (*A)[n];
    }
  }
  if (*nnz) cnt++;
  *nnz = cnt;

  if(rank==0) printf("done.\n");

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Type_free(&MPI_NONZERO_T);

  free(sendNonZeros);
  free(globalNumbering);

  free(AsendCounts);
  free(ArecvCounts);
  free(AsendOffsets);
  free(ArecvOffsets);
}