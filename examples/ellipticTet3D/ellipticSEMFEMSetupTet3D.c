#include "ellipticTet3D.h"

typedef struct{

  dfloat VX;
  dfloat VY;
  dfloat VZ;

  dlong localId;
  hlong globalId;

}FEMverts_t;

typedef struct {

  dlong localId;
  hlong globalId;
  int ownerRank;

}parallelNode_t;

// compare on global owners
int parallelCompareOwnersAndGlobalId(const void *a, const void *b);

// compare on global indices
int parallelCompareGlobalId(const void *a, const void *b);

// compare xy coordinates
int parallelCompareFEMvertsLocation(const void *a, const void *b){
  dfloat NODETOL = 1e-6;

  FEMverts_t *fa = (FEMverts_t*) a;
  FEMverts_t *fb = (FEMverts_t*) b;

  if(fa->VX < fb->VX - NODETOL) return -1;
  if(fa->VX > fb->VX + NODETOL) return +1;

  if(fa->VY < fb->VY - NODETOL) return -1;
  if(fa->VY > fb->VY + NODETOL) return +1;

  if(fa->VZ < fb->VZ - NODETOL) return -1;
  if(fa->VZ > fb->VZ + NODETOL) return +1;

  return 0;
}

// compare local id
int parallelCompareFEMvertsLocalId(const void *a, const void *b){

  FEMverts_t *fa = (FEMverts_t*) a;
  FEMverts_t *fb = (FEMverts_t*) b;

  if(fa->localId < fb->localId) return -1;
  if(fa->localId > fb->localId) return +1;

  return 0;
}

int parallelCompareRowColumn(const void *a, const void *b);

void ellipticSEMFEMSetupTet3D(solver_t *solver, precon_t* precon,
                              dfloat tau, dfloat lambda, int *BCType,
                              const char *options, const char *parAlmondOptions) {

  if (!strstr(options, "CONTINUOUS")) {
    printf("SEMFEM is supported for CONTINUOUS only\n");
    exit(-1);
  }

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh3D* mesh = solver->mesh; //original mesh

  mesh3D* pmesh = (mesh3D*) calloc (1,sizeof(mesh3D)); //partially assembled fem mesh (result of projecting sem element to larger space)

  precon->femMesh = (mesh3D*) calloc (1,sizeof(mesh3D)); //full fem mesh
  mesh3D *femMesh = precon->femMesh;

  memcpy(pmesh  ,mesh,sizeof(mesh3D));
  memcpy(femMesh,mesh,sizeof(mesh3D));

  //set semfem nodes as the grid points
  pmesh->Np = mesh->NpFEM;
  pmesh->r  = mesh->rFEM;
  pmesh->s  = mesh->sFEM;
  pmesh->t  = mesh->tFEM;

  //count number of face nodes in the semfem element
  dfloat NODETOL = 1e-6;
  pmesh->Nfp=0;
  for (int n=0;n<pmesh->Np;n++)
    if (fabs(pmesh->t[n]+1)<NODETOL) pmesh->Nfp++;

  //remake the faceNodes array
  pmesh->faceNodes = (int *) calloc(pmesh->Nfaces*pmesh->Nfp,sizeof(int));
  int f0=0, f1=0, f2=0, f3=0;
  for (int n=0;n<pmesh->Np;n++) {
    if (fabs(pmesh->t[n]+1)<NODETOL)           pmesh->faceNodes[0*pmesh->Nfp+f0++] = n;
    if (fabs(pmesh->s[n]+1)<NODETOL)           pmesh->faceNodes[1*pmesh->Nfp+f1++] = n;
    if (fabs(pmesh->r[n]+pmesh->s[n]+
                     pmesh->t[n]+1.0)<NODETOL) pmesh->faceNodes[2*pmesh->Nfp+f2++] = n;
    if (fabs(pmesh->r[n]+1)<NODETOL)           pmesh->faceNodes[3*pmesh->Nfp+f3++] = n;
  }

  //remake vertexNodes array
  pmesh->vertexNodes = (int*) calloc(pmesh->Nverts, sizeof(int));
  for(int n=0;n<pmesh->Np;++n){
    if( (pmesh->r[n]+1)*(pmesh->r[n]+1)+(pmesh->s[n]+1)*(pmesh->s[n]+1)+(pmesh->t[n]+1)*(pmesh->t[n]+1)<NODETOL)
      pmesh->vertexNodes[0] = n;
    if( (pmesh->r[n]-1)*(pmesh->r[n]-1)+(pmesh->s[n]+1)*(pmesh->s[n]+1)+(pmesh->t[n]+1)*(pmesh->t[n]+1)<NODETOL)
      pmesh->vertexNodes[1] = n;
    if( (pmesh->r[n]+1)*(pmesh->r[n]+1)+(pmesh->s[n]-1)*(pmesh->s[n]-1)+(pmesh->t[n]+1)*(pmesh->t[n]+1)<NODETOL)
      pmesh->vertexNodes[2] = n;
    if( (pmesh->r[n]+1)*(pmesh->r[n]+1)+(pmesh->s[n]+1)*(pmesh->s[n]+1)+(pmesh->t[n]-1)*(pmesh->t[n]-1)<NODETOL)
      pmesh->vertexNodes[3] = n;
  }

  // connect elements using parallel sort
  meshParallelConnect(pmesh);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesTet3D(pmesh);

  // free(sendBuffer);
  meshHaloSetup(pmesh);

  // connect face nodes (find trace indices)
  meshConnectFaceNodes3D(pmesh);

  // global nodes
  meshParallelConnectNodes(pmesh);
  //pmesh->globalIds and pmesh->globalOwners are now populated

  //now build the full degree 1 fem mesh
  int femN = 1; //degree of fem approximation

  /* allocate space for node coordinates */
  femMesh->Nelements = mesh->NelFEM*mesh->Nelements;
  femMesh->EToV = (hlong*) calloc(femMesh->Nelements*femMesh->Nverts, sizeof(hlong));
  femMesh->EX = (dfloat*) calloc(femMesh->Nverts*femMesh->Nelements, sizeof(dfloat));
  femMesh->EY = (dfloat*) calloc(femMesh->Nverts*femMesh->Nelements, sizeof(dfloat));
  femMesh->EZ = (dfloat*) calloc(femMesh->Nverts*femMesh->Nelements, sizeof(dfloat));

  dlong *localIds = (dlong *) calloc(femMesh->Nverts*femMesh->Nelements,sizeof(dlong));

  // dlong NFEMverts = mesh->Nelements*mesh->NpFEM;
  for(dlong e=0;e<mesh->Nelements;++e){
    for (int n=0;n<mesh->NelFEM;n++) {
      //local ids in the subelement fem grid
      dlong id1 = e*mesh->NpFEM + mesh->FEMEToV[n*mesh->Nverts+0];
      dlong id2 = e*mesh->NpFEM + mesh->FEMEToV[n*mesh->Nverts+1];
      dlong id3 = e*mesh->NpFEM + mesh->FEMEToV[n*mesh->Nverts+2];
      dlong id4 = e*mesh->NpFEM + mesh->FEMEToV[n*mesh->Nverts+3];

      #if 1
      // check orientation
      dfloat xe1 = pmesh->x[id1], xe2 = pmesh->x[id2], xe3 = pmesh->x[id3], xe4 = pmesh->x[id4];
      dfloat ye1 = pmesh->y[id1], ye2 = pmesh->y[id2], ye3 = pmesh->y[id3], ye4 = pmesh->y[id4];
      dfloat ze1 = pmesh->z[id1], ze2 = pmesh->z[id2], ze3 = pmesh->z[id3], ze4 = pmesh->z[id4];
      
      /* Jacobian matrix */
      dfloat xr = 0.5*(xe2-xe1), xs = 0.5*(xe3-xe1), xt = 0.5*(xe4-xe1);
      dfloat yr = 0.5*(ye2-ye1), ys = 0.5*(ye3-ye1), yt = 0.5*(ye4-ye1);
      dfloat zr = 0.5*(ze2-ze1), zs = 0.5*(ze3-ze1), zt = 0.5*(ze4-ze1);

      /* compute geometric factors for affine coordinate transform*/
      dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
      if(J<0){
        printf("Uh oh\n");
        dlong id3tmp = id3;
        id3 = id4;
        id4 = id3tmp;
      }
      #endif

      /* read vertex triplet for triangle */
      dlong femId = e*mesh->NelFEM*mesh->Nverts+n*mesh->Nverts;
      femMesh->EToV[femId+0] = pmesh->globalIds[id1];
      femMesh->EToV[femId+1] = pmesh->globalIds[id2];
      femMesh->EToV[femId+2] = pmesh->globalIds[id3];
      femMesh->EToV[femId+3] = pmesh->globalIds[id4];

      femMesh->EX[femId+0] = pmesh->x[id1];
      femMesh->EX[femId+1] = pmesh->x[id2];
      femMesh->EX[femId+2] = pmesh->x[id3];
      femMesh->EX[femId+3] = pmesh->x[id4];

      femMesh->EY[femId+0] = pmesh->y[id1];
      femMesh->EY[femId+1] = pmesh->y[id2];
      femMesh->EY[femId+2] = pmesh->y[id3];
      femMesh->EY[femId+3] = pmesh->y[id4];

      femMesh->EZ[femId+0] = pmesh->z[id1];
      femMesh->EZ[femId+1] = pmesh->z[id2];
      femMesh->EZ[femId+2] = pmesh->z[id3];
      femMesh->EZ[femId+3] = pmesh->z[id4];

      localIds[femId+0] = id1;
      localIds[femId+1] = id2;
      localIds[femId+2] = id3;
      localIds[femId+3] = id4;
    }
  }

  // connect elements using parallel sort
  meshParallelConnect(femMesh);

  // load reference (r,s) element nodes
  meshLoadReferenceNodesTet3D(femMesh, femN);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesTet3D(femMesh);

  // compute geometric factors
  meshGeometricFactorsTet3D(femMesh);

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetup(femMesh);

  // connect face nodes (find trace indices)
  meshConnectFaceNodes3D(femMesh);

  // compute surface geofacs
  meshSurfaceGeometricFactorsTet3D(femMesh);

  // global nodes
  meshParallelConnectNodes(femMesh);


  //on-host version of gather-scatter
  int verbose = strstr(options,"VERBOSE") ? 1:0;
  pmesh->hostGsh = gsParallelGatherScatterSetup(pmesh->Nelements*pmesh->Np, pmesh->globalIds,verbose);

  //make a node-wise bc flag using the gsop (prioritize Dirichlet boundaries over Neumann)
  int *mapB = (int *) calloc(pmesh->Nelements*pmesh->Np,sizeof(int));
  for (dlong e=0;e<pmesh->Nelements;e++) {
    for (int n=0;n<pmesh->Np;n++) mapB[n+e*pmesh->Np] = 1E9;
    for (int f=0;f<pmesh->Nfaces;f++) {
      int bc = pmesh->EToB[f+e*pmesh->Nfaces];
      if (bc>0) {
        for (int n=0;n<pmesh->Nfp;n++) {
          int BCFlag = BCType[bc];
          int fid = pmesh->faceNodes[n+f*pmesh->Nfp];
          mapB[fid+e*pmesh->Np] = mymin(BCFlag,mapB[fid+e*pmesh->Np]);
        }
      }
    }
  }
  gsParallelGatherScatter(pmesh->hostGsh, mapB, "int", "min");

  //use the bc flags to find masked ids
  for (dlong n=0;n<pmesh->Nelements*pmesh->Np;n++) {
    if (mapB[n] == 1) { //Dirichlet boundary
      pmesh->globalIds[n] = -1;
    }
  }

  // squeeze node numbering
  hlong *globalStarts = (hlong*) calloc(size+1, sizeof(hlong));
  meshParallelConsecutiveGlobalNumbering(pmesh, pmesh->Np*pmesh->Nelements, pmesh->globalIds, pmesh->globalOwners, globalStarts);

  for (dlong n=0;n<pmesh->Np*pmesh->Nelements;n++) {
    dlong id = pmesh->gatherLocalIds[n];
    pmesh->gatherBaseIds[n] = pmesh->globalIds[id];
  }

  //build gather scatter with masked nodes
  precon->FEMogs = meshParallelGatherScatterSetup(pmesh, pmesh->Nelements*pmesh->Np, 
                                        pmesh->gatherLocalIds,  pmesh->gatherBaseIds, 
                                        pmesh->gatherBaseRanks, pmesh->gatherHaloFlags,verbose);



  //dont need these anymore
  free(pmesh->vmapM);
  free(pmesh->vmapP);
  free(pmesh->mapP);
  //maybe more cleanup can go here

  //build stiffness matrices
  femMesh->Srr = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
  femMesh->Srs = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
  femMesh->Srt = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
  femMesh->Ssr = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
  femMesh->Sss = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
  femMesh->Sst = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
  femMesh->Str = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
  femMesh->Sts = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
  femMesh->Stt = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
  for (int n=0;n<femMesh->Np;n++) {
    for (int m=0;m<femMesh->Np;m++) {
      for (int k=0;k<femMesh->Np;k++) {
        for (int l=0;l<femMesh->Np;l++) {
          femMesh->Srr[m+n*femMesh->Np] += femMesh->Dr[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dr[m+k*femMesh->Np];
          femMesh->Srs[m+n*femMesh->Np] += femMesh->Dr[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Ds[m+k*femMesh->Np];
          femMesh->Srt[m+n*femMesh->Np] += femMesh->Dr[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dt[m+k*femMesh->Np];
          femMesh->Ssr[m+n*femMesh->Np] += femMesh->Ds[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dr[m+k*femMesh->Np];
          femMesh->Sss[m+n*femMesh->Np] += femMesh->Ds[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Ds[m+k*femMesh->Np];
          femMesh->Sst[m+n*femMesh->Np] += femMesh->Ds[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dt[m+k*femMesh->Np];
          femMesh->Str[m+n*femMesh->Np] += femMesh->Dt[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dr[m+k*femMesh->Np];
          femMesh->Sts[m+n*femMesh->Np] += femMesh->Dt[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Ds[m+k*femMesh->Np];
          femMesh->Stt[m+n*femMesh->Np] += femMesh->Dt[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dt[m+k*femMesh->Np];
        }
      }
    }
  }

  printf("Building full SEMFEM matrix..."); fflush(stdout);

  // Build non-zeros of stiffness matrix (unassembled)
  long long int nnzLocal = femMesh->Np*femMesh->Np*femMesh->Nelements;

  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocal, sizeof(nonZero_t));
  int *AsendCounts  = (int*) calloc(size, sizeof(int));
  int *ArecvCounts  = (int*) calloc(size, sizeof(int));
  int *AsendOffsets = (int*) calloc(size+1, sizeof(int));
  int *ArecvOffsets = (int*) calloc(size+1, sizeof(int));

  //Build unassembed non-zeros
  dlong cnt =0;
  #pragma omp parallel for
  for (dlong e=0;e<femMesh->Nelements;e++) {

    dfloat Grr = femMesh->ggeo[e*femMesh->Nggeo + G00ID];
    dfloat Grs = femMesh->ggeo[e*femMesh->Nggeo + G01ID];
    dfloat Grt = femMesh->ggeo[e*femMesh->Nggeo + G02ID];
    dfloat Gss = femMesh->ggeo[e*femMesh->Nggeo + G11ID];
    dfloat Gst = femMesh->ggeo[e*femMesh->Nggeo + G12ID];
    dfloat Gtt = femMesh->ggeo[e*femMesh->Nggeo + G22ID];
    dfloat J   = femMesh->ggeo[e*femMesh->Nggeo + GWJID];

    for (int n=0;n<femMesh->Np;n++) {
      dlong idn = localIds[e*femMesh->Np + n];
      if (pmesh->globalIds[idn]<0) continue; //skip masked nodes
      for (int m=0;m<femMesh->Np;m++) {
        dlong idm = localIds[e*femMesh->Np + m];
        if (pmesh->globalIds[idm]<0) continue; //skip masked nodes

        dfloat val = 0.;
        val += Grr*femMesh->Srr[m+n*femMesh->Np];
        val += Grs*femMesh->Srs[m+n*femMesh->Np];
        val += Grt*femMesh->Srt[m+n*femMesh->Np];
        val += Grs*femMesh->Ssr[m+n*femMesh->Np];
        val += Gss*femMesh->Sss[m+n*femMesh->Np];
        val += Gst*femMesh->Sst[m+n*femMesh->Np];
        val += Grt*femMesh->Str[m+n*femMesh->Np];
        val += Gst*femMesh->Sts[m+n*femMesh->Np];
        val += Gtt*femMesh->Stt[m+n*femMesh->Np];
        val += J*lambda*femMesh->MM[m+n*femMesh->Np];

        dfloat nonZeroThreshold = 1e-7;
        if (fabs(val)>nonZeroThreshold) {
          #pragma omp critical
          {
            // pack non-zero
            sendNonZeros[cnt].val = val;
            sendNonZeros[cnt].row = pmesh->globalIds[idn];
            sendNonZeros[cnt].col = pmesh->globalIds[idm];
            sendNonZeros[cnt].ownerRank = pmesh->globalOwners[idn];
            cnt++;
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

  //tell parAlmond not to gather this level (its done manually)
  agmgLevel *baseLevel = precon->parAlmond->levels[0];
  baseLevel->gatherLevel = false;
  baseLevel->weightedInnerProds = false;

  // build interp and anterp
  dfloat *SEMFEMAnterp = (dfloat*) calloc(mesh->NpFEM*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->NpFEM;++n){
    for(int m=0;m<mesh->Np;++m){
      SEMFEMAnterp[n+m*mesh->NpFEM] = mesh->SEMFEMInterp[n*mesh->Np+m];
    }
  }

  mesh->o_SEMFEMInterp = mesh->device.malloc(mesh->NpFEM*mesh->Np*sizeof(dfloat),mesh->SEMFEMInterp);
  mesh->o_SEMFEMAnterp = mesh->device.malloc(mesh->NpFEM*mesh->Np*sizeof(dfloat),SEMFEMAnterp);

  free(SEMFEMAnterp);

  precon->o_rFEM = mesh->device.malloc(mesh->Nelements*mesh->NpFEM*sizeof(dfloat));
  precon->o_zFEM = mesh->device.malloc(mesh->Nelements*mesh->NpFEM*sizeof(dfloat));

  precon->o_GrFEM = mesh->device.malloc(precon->FEMogs->Ngather*sizeof(dfloat));
  precon->o_GzFEM = mesh->device.malloc(precon->FEMogs->Ngather*sizeof(dfloat));
}