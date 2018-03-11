#include "ellipticQuad2D.h"

int parallelCompareRowColumn(const void *a, const void *b);

void ellipticSEMFEMSetupQuad2D(solver_t *solver, precon_t* precon,
                              dfloat tau, dfloat lambda, int *BCType,
                              const char *options, const char *parAlmondOptions) {

  if (!strstr(options, "CONTINUOUS")) {
    printf("SEMFEM is supported for CONTINUOUS only\n");
    exit(-1);
  }

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D* mesh = solver->mesh; //original mesh

  precon->femMesh = (mesh2D*) calloc (1,sizeof(mesh2D)); //full fem mesh
  mesh2D *femMesh = precon->femMesh;

  memcpy(femMesh,mesh,sizeof(mesh2D));

  //now build the full degree 1 fem mesh
  int femN = 1; //degree of fem approximation

  /* allocate space for node coordinates */
  femMesh->Nelements = mesh->NelFEM*mesh->Nelements;
  femMesh->EToV = (hlong*) calloc(femMesh->Nelements*femMesh->Nverts, sizeof(hlong));
  femMesh->EX = (dfloat*) calloc(femMesh->Nverts*femMesh->Nelements, sizeof(dfloat));
  femMesh->EY = (dfloat*) calloc(femMesh->Nverts*femMesh->Nelements, sizeof(dfloat));

  dlong *localIds = (dlong *) calloc(femMesh->Nverts*femMesh->Nelements,sizeof(dlong));

  // dlong NFEMverts = mesh->Nelements*mesh->NpFEM;
  for(dlong e=0;e<mesh->Nelements;++e){
    for (int n=0;n<mesh->NelFEM;n++) {
      //local ids in the subelement fem grid
      dlong id1 = e*mesh->Np + mesh->FEMEToV[n*mesh->Nverts+0];
      dlong id2 = e*mesh->Np + mesh->FEMEToV[n*mesh->Nverts+1];
      dlong id3 = e*mesh->Np + mesh->FEMEToV[n*mesh->Nverts+2];
      dlong id4 = e*mesh->Np + mesh->FEMEToV[n*mesh->Nverts+3];

      /* read vertex quadruplet for quad */
      dlong femId = e*mesh->NelFEM*mesh->Nverts + n*mesh->Nverts;
      femMesh->EToV[femId+0] = mesh->globalIds[id1];
      femMesh->EToV[femId+1] = mesh->globalIds[id2];
      femMesh->EToV[femId+2] = mesh->globalIds[id3];
      femMesh->EToV[femId+3] = mesh->globalIds[id4];

      femMesh->EX[femId+0] = mesh->x[id1];
      femMesh->EX[femId+1] = mesh->x[id2];
      femMesh->EX[femId+2] = mesh->x[id3];
      femMesh->EX[femId+3] = mesh->x[id4];

      femMesh->EY[femId+0] = mesh->y[id1];
      femMesh->EY[femId+1] = mesh->y[id2];
      femMesh->EY[femId+2] = mesh->y[id3];
      femMesh->EY[femId+3] = mesh->y[id4];

      localIds[femId+0] = id1;
      localIds[femId+1] = id2;
      localIds[femId+2] = id4;  //need to swap this as the Np nodes are ordered [0,1,4,3] in a degree 1 element
      localIds[femId+3] = id3;
    }
  }

  // connect elements using parallel sort
  meshParallelConnect(femMesh);

  // load reference (r,s) element nodes
  meshLoadReferenceNodesQuad2D(femMesh, femN);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesQuad2D(femMesh);

  // compute geometric factors
  meshGeometricFactorsQuad2D(femMesh);

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetup(femMesh);

  // connect face nodes (find trace indices)
  meshConnectFaceNodes2D(femMesh);

  // compute surface geofacs
  meshSurfaceGeometricFactorsQuad2D(femMesh);

  // global nodes
  meshParallelConnectNodes(femMesh);


  dlong Ntotal = mesh->Np*mesh->Nelements;

  hlong *globalNumbering = (hlong *) calloc(Ntotal,sizeof(hlong));
  hlong *globalStarts = (hlong *) calloc(size+1,sizeof(hlong));
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
  precon->FEMogs = meshParallelGatherScatterSetup(mesh, Ntotal, 
                                        mesh->gatherLocalIds,  gatherMaskedBaseIds, 
                                        mesh->gatherBaseRanks, mesh->gatherHaloFlags,verbose);

  printf("Building full SEMFEM matrix..."); fflush(stdout);

  // Build non-zeros of stiffness matrix (unassembled)
  dlong nnzLocal = femMesh->Np*femMesh->Np*femMesh->Nelements;

  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocal, sizeof(nonZero_t));
  int *AsendCounts  = (int*) calloc(size, sizeof(int));
  int *ArecvCounts  = (int*) calloc(size, sizeof(int));
  int *AsendOffsets = (int*) calloc(size+1, sizeof(int));
  int *ArecvOffsets = (int*) calloc(size+1, sizeof(int));

  //Build unassembed non-zeros
  dlong cnt =0;
  #pragma omp parallel for
  for (dlong e=0;e<femMesh->Nelements;e++) {
    for (int ny=0;ny<femMesh->Nq;ny++) {
      for (int nx=0;nx<femMesh->Nq;nx++) {

        dlong idn = localIds[e*femMesh->Np + nx+ny*femMesh->Nq];
        if (globalNumbering[idn]<0) continue; //skip masked nodes

        for (int my=0;my<femMesh->Nq;my++) {
          for (int mx=0;mx<femMesh->Nq;mx++) {

            dlong idm = localIds[e*femMesh->Np + mx+my*femMesh->Nq];
            if (globalNumbering[idm]<0) continue; //skip masked nodes

            int id;
            dfloat val = 0.;

            if (ny==my) {
              for (int k=0;k<femMesh->Nq;k++) {
                id = k+ny*femMesh->Nq;
                dfloat Grr = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + G00ID*femMesh->Np];

                val += Grr*femMesh->D[nx+k*femMesh->Nq]*femMesh->D[mx+k*femMesh->Nq];
              }
            }

            id = mx+ny*femMesh->Nq;
            dfloat Grs = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + G01ID*femMesh->Np];
            val += Grs*femMesh->D[nx+mx*femMesh->Nq]*femMesh->D[my+ny*femMesh->Nq];

            id = nx+my*femMesh->Nq;
            dfloat Gsr = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + G01ID*femMesh->Np];
            val += Gsr*femMesh->D[mx+nx*femMesh->Nq]*femMesh->D[ny+my*femMesh->Nq];

            if (nx==mx) {
              for (int k=0;k<femMesh->Nq;k++) {
                id = nx+k*femMesh->Nq;
                dfloat Gss = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + G11ID*femMesh->Np];

                val += Gss*femMesh->D[ny+k*femMesh->Nq]*femMesh->D[my+k*femMesh->Nq];
              }
            }

            if ((nx==mx)&&(ny==my)) {
              id = nx + ny*femMesh->Nq;
              dfloat JW = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + GWJID*femMesh->Np];
              val += JW*lambda;
            }

            dfloat nonZeroThreshold = 1e-7;
            if (fabs(val)>nonZeroThreshold) {
              #pragma omp critical
              {
                // pack non-zero
                sendNonZeros[cnt].val = val;
                sendNonZeros[cnt].row = globalNumbering[idn];
                sendNonZeros[cnt].col = globalNumbering[idm];
                sendNonZeros[cnt].ownerRank = mesh->globalOwners[idn];
                cnt++;
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
  dlong nnz = 0;
  for(int r=0;r<size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    nnz += ArecvCounts[r];
  }

  nonZero_t *A = (nonZero_t*) calloc(nnz, sizeof(nonZero_t));

  // determine number to receive
  MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_NONZERO_T,
                           A, ArecvCounts, ArecvOffsets, MPI_NONZERO_T,
                           MPI_COMM_WORLD);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort(A, nnz, sizeof(nonZero_t), parallelCompareRowColumn);

  // compress duplicates
  cnt = 0;
  for(dlong n=1;n<nnz;++n){
    if(A[n].row == A[cnt].row && A[n].col == A[cnt].col){
      A[cnt].val += A[n].val;
    } else{
      ++cnt;
      A[cnt] = A[n];
    }
  }
  if (nnz) cnt++;
  nnz = cnt;

  if(rank==0) printf("done.\n");

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Type_free(&MPI_NONZERO_T);

  hlong *Rows = (hlong *) calloc(nnz, sizeof(hlong));
  hlong *Cols = (hlong *) calloc(nnz, sizeof(hlong));
  dfloat *Vals = (dfloat*) calloc(nnz,sizeof(dfloat));

  for (dlong n=0;n<nnz;n++) {
    Rows[n] = A[n].row;
    Cols[n] = A[n].col;
    Vals[n] = A[n].val;
  }

  precon->parAlmond = parAlmondInit(mesh, parAlmondOptions);
  parAlmondAgmgSetup(precon->parAlmond,
                     globalStarts,
                     nnz,
                     Rows,
                     Cols,
                     Vals,
                     solver->allNeumann,
                     solver->allNeumannPenalty);
  free(A); free(Rows); free(Cols); free(Vals);

  //tell parAlmond to gather this level
  agmgLevel *baseLevel = precon->parAlmond->levels[0];

  baseLevel->gatherLevel = true;
  baseLevel->Srhs = (dfloat*) calloc(mesh->Np*mesh->Nelements,sizeof(dfloat));
  baseLevel->Sx   = (dfloat*) calloc(mesh->Np*mesh->Nelements,sizeof(dfloat));
  baseLevel->o_Srhs = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat));
  baseLevel->o_Sx   = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat));

  baseLevel->weightedInnerProds = false;

  baseLevel->gatherArgs = (void **) calloc(4,sizeof(void*));  
  baseLevel->gatherArgs[0] = (void *) solver;
  baseLevel->gatherArgs[1] = (void *) precon->FEMogs;  //use the gs made from the partial gathered femgrid 
  baseLevel->gatherArgs[2] = (void *) &(baseLevel->o_Sx);
  baseLevel->gatherArgs[3] = (void *) options;
  baseLevel->scatterArgs = baseLevel->gatherArgs;

  baseLevel->device_gather  = ellipticGather;
  baseLevel->device_scatter = ellipticScatter;  
}