#include "ellipticTri2D.h"

typedef struct{

  dfloat VX;
  dfloat VY;

  iint localId;
  iint globalId;

}FEMverts_t;

typedef struct {

  iint localId;
  iint globalId;
  iint ownerRank;

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

void ellipticSEMFEMSetupTri2D(solver_t *solver, precon_t* precon,
                              dfloat tau, dfloat lambda, iint *BCType,
                              const char *options, const char *parAlmondOptions) {

  if (!strstr(options, "CONTINUOUS")) {
    printf("SEMFEM is supported for CONTINUOUS only\n");
    exit(-1);
  }

  printf("Setting up FEM problem...");

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D* mesh = solver->mesh;

  precon->femMesh = (mesh2D*) calloc (1,sizeof(mesh2D));
  mesh2D *femMesh = precon->femMesh;
  memcpy(femMesh,mesh,sizeof(mesh2D));

  int femN = 1; //degree of fem approximation

  //making the fem grid on each element will make duplicate vertices,
  // We need a global ordering
  iint NFEMverts = mesh->Nelements*mesh->NpFEM;
  FEMverts_t *FEMverts = (FEMverts_t *) calloc(NFEMverts,sizeof(FEMverts_t));
  for(int e=0;e<mesh->Nelements;++e){
    iint id = e*mesh->Nverts;

    dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];

    dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];

    for (int n=0;n<mesh->NpFEM;n++) {
      dfloat rn = mesh->rFEM[n];
      dfloat sn = mesh->sFEM[n];

      iint fid = e*mesh->NpFEM+n;
      FEMverts[fid].localId = fid;
      FEMverts[fid].VX = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
      FEMverts[fid].VY = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
    }
  }

  // sort based on x,y coordinates
  qsort(FEMverts, NFEMverts, sizeof(FEMverts_t), parallelCompareFEMvertsLocation);

  //compress the duplicates to make a global numbering
  iint NFEMnodes = 0;
  for (iint n=1;n<NFEMverts;n++) {
    if (parallelCompareFEMvertsLocation(FEMverts+n-1, FEMverts+n)<0) NFEMnodes++;

    FEMverts[n].globalId = NFEMnodes;
  }
  if (NFEMverts) NFEMnodes++;

  qsort(FEMverts, NFEMverts, sizeof(FEMverts_t), parallelCompareFEMvertsLocalId);


  /* allocate space for node coordinates */
  femMesh->Nelements = mesh->NelFEM*mesh->Nelements;
  femMesh->EToV = (iint*) calloc(femMesh->Nelements*femMesh->Nverts, sizeof(iint));
  femMesh->EX = (dfloat*) calloc(femMesh->Nverts*femMesh->Nelements, sizeof(dfloat));
  femMesh->EY = (dfloat*) calloc(femMesh->Nverts*femMesh->Nelements, sizeof(dfloat));

  int localNodeCount = femMesh->Nelements*femMesh->Np;
  femMesh->gatherLocalIds  = (iint*) calloc(localNodeCount, sizeof(iint));
  femMesh->gatherBaseIds   = (iint*) calloc(localNodeCount, sizeof(iint));

  for(int e=0;e<mesh->Nelements;++e){
    for (int n=0;n<mesh->NelFEM;n++) {
      int id1 = e*mesh->NpFEM + mesh->FEMEToV[n*mesh->Nverts+0];
      int id2 = e*mesh->NpFEM + mesh->FEMEToV[n*mesh->Nverts+1];
      int id3 = e*mesh->NpFEM + mesh->FEMEToV[n*mesh->Nverts+2];

      // check orientation
      dfloat xe1 = FEMverts[id1].VX, xe2 = FEMverts[id2].VX, xe3 = FEMverts[id3].VX;
      dfloat ye1 = FEMverts[id1].VY, ye2 = FEMverts[id2].VY, ye3 = FEMverts[id3].VY;
      dfloat J = 0.25*((xe2-xe1)*(ye3-ye1) - (xe3-xe1)*(ye2-ye1));
      if(J<0){
        iint id3tmp = id3;
        id3 = id2;
        id2 = id3tmp;
      }

      /* read vertex triplet for triangle */
      iint femId = e*mesh->NelFEM*mesh->Nverts+n*mesh->Nverts;
      femMesh->EToV[femId+0] = FEMverts[id1].globalId;
      femMesh->EToV[femId+1] = FEMverts[id2].globalId;
      femMesh->EToV[femId+2] = FEMverts[id3].globalId;

      femMesh->EX[femId+0] = FEMverts[id1].VX;
      femMesh->EX[femId+1] = FEMverts[id2].VX;
      femMesh->EX[femId+2] = FEMverts[id3].VX;

      femMesh->EY[femId+0] = FEMverts[id1].VY;
      femMesh->EY[femId+1] = FEMverts[id2].VY;
      femMesh->EY[femId+2] = FEMverts[id3].VY;

      //propagate the global node numbering
      iint id = e*mesh->NelFEM*mesh->Nverts+n*mesh->Nverts;
      femMesh->gatherLocalIds[id+0] = id+0;
      femMesh->gatherLocalIds[id+1] = id+1;
      femMesh->gatherLocalIds[id+2] = id+2;

      femMesh->gatherBaseIds[id+0] = FEMverts[id1].globalId;
      femMesh->gatherBaseIds[id+1] = FEMverts[id2].globalId;
      femMesh->gatherBaseIds[id+2] = FEMverts[id3].globalId;
    }
  }

  //make a mapping to the original element face number to propagate boundary info
  dfloat NODETOL = 1e-6;
  int *faceMap = (int *) calloc(mesh->NelFEM*mesh->Nfaces,sizeof(int));
  for (int n=0;n<mesh->NelFEM;n++) {
    for (int f=0;f<mesh->Nfaces;f++) {
      int id1 = f;
      int id2 = (f < mesh->Nfaces-1) ? f+1 : 0;

      int v1 = mesh->FEMEToV[n*mesh->Nverts+id1];
      int v2 = mesh->FEMEToV[n*mesh->Nverts+id2];

      dfloat r1 = mesh->rFEM[v1];
      dfloat s1 = mesh->sFEM[v1];

      dfloat r2 = mesh->rFEM[v2];
      dfloat s2 = mesh->sFEM[v2];

      if ((fabs(s1+1)<NODETOL)&&(fabs(s2+1)<NODETOL)) {
        faceMap[n*mesh->Nfaces+f] = 0;
      } else if ((fabs(s1+r1)<NODETOL)&&(fabs(s2+r2)<NODETOL)) {
        faceMap[n*mesh->Nfaces+f] = 1;
      } else if ((fabs(r1+1)<NODETOL)&&(fabs(r2+1)<NODETOL)) {
        faceMap[n*mesh->Nfaces+f] = 2;
      } else {
        faceMap[n*mesh->Nfaces+f] = -1;
      }
    }
  }

  // connect elements using parallel sort
  meshParallelConnect(femMesh);

  //propagate boundary flags
  femMesh->EToB = (int*) calloc(femMesh->Nelements*mesh->Nfaces, sizeof(int));
  for(int n=0;n<femMesh->Nelements*mesh->Nfaces;++n) femMesh->EToB[n] = -1;

  for (iint e=0;e<mesh->Nelements;e++) {
    for (int f=0;f<mesh->Nfaces;f++) {
      if (mesh->EToB[e*mesh->Nfaces+f]>0) {
        for (int n=0;n<mesh->NelFEM*mesh->Nfaces;n++) {
          if (faceMap[n]==f)
              femMesh->EToB[e*mesh->NelFEM*mesh->Nfaces + n] = mesh->EToB[e*mesh->Nfaces+f];
        }
      }
    }
  }

  // load reference (r,s) element nodes
  meshLoadReferenceNodesTri2D(femMesh, femN);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesTri2D(femMesh);

  // compute geometric factors
  meshGeometricFactorsTri2D(femMesh);

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetup(femMesh);

  // connect face nodes (find trace indices)
  meshConnectFaceNodes2D(femMesh);

  // compute surface geofacs
  meshSurfaceGeometricFactorsTri2D(femMesh);

  //build stiffness matrices
  femMesh->Srr = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
  femMesh->Srs = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
  femMesh->Ssr = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
  femMesh->Sss = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
  for (iint n=0;n<femMesh->Np;n++) {
    for (iint m=0;m<femMesh->Np;m++) {
      for (iint k=0;k<femMesh->Np;k++) {
        for (iint l=0;l<femMesh->Np;l++) {
          femMesh->Srr[m+n*femMesh->Np] += femMesh->Dr[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dr[m+k*femMesh->Np];
          femMesh->Srs[m+n*femMesh->Np] += femMesh->Dr[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Ds[m+k*femMesh->Np];
          femMesh->Ssr[m+n*femMesh->Np] += femMesh->Ds[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dr[m+k*femMesh->Np];
          femMesh->Sss[m+n*femMesh->Np] += femMesh->Ds[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Ds[m+k*femMesh->Np];
        }
      }
    }
  }

  iint nnz;
  nonZero_t *A;

  iint Nnum = femMesh->Np*femMesh->Nelements;
  iint *globalStarts = (iint*) calloc(size+1, sizeof(iint));

  printf("done. \n");

  ellipticBuildContinuousTri2D(femMesh,lambda,&A,&nnz,&(precon->hgs),globalStarts,options);

  iint *Rows = (iint *) calloc(nnz, sizeof(iint));
  iint *Cols = (iint *) calloc(nnz, sizeof(iint));
  dfloat *Vals = (dfloat*) calloc(nnz,sizeof(dfloat));

  for (iint n=0;n<nnz;n++) {
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

  // build interp and anterp
  dfloat *SEMFEMAnterp = (dfloat*) calloc(mesh->NpFEM*mesh->Np, sizeof(dfloat));
  for(iint n=0;n<mesh->NpFEM;++n){
    for(iint m=0;m<mesh->Np;++m){
      SEMFEMAnterp[n+m*mesh->NpFEM] = mesh->SEMFEMInterp[n*mesh->Np+m];
    }
  }

  mesh->o_SEMFEMInterp = mesh->device.malloc(mesh->NpFEM*mesh->Np*sizeof(dfloat),mesh->SEMFEMInterp);
  mesh->o_SEMFEMAnterp = mesh->device.malloc(mesh->NpFEM*mesh->Np*sizeof(dfloat),SEMFEMAnterp);

  free(SEMFEMAnterp);

  //correct the gather operation (since the fem grid is partially gathered after interpolation)
  Nnum = mesh->NpFEM*mesh->Nelements;
  iint *globalOwners = (iint*) calloc(Nnum, sizeof(iint));
  iint *globalNumbering = (iint*) calloc(Nnum, sizeof(iint));

  for (iint n=0;n<Nnum;n++) {
    globalNumbering[n] = FEMverts[n].globalId;
  }

  dfloat *femMask; //TODO need to fill this

  // squeeze node numbering
  meshParallelConsecutiveGlobalNumbering(Nnum, globalNumbering, globalOwners, globalStarts,femMask);

  //use the ordering to define a gather+scatter for assembly
  precon->hgs = meshParallelGatherSetup(mesh, Nnum, globalNumbering, globalOwners);

  precon->o_rFEM = mesh->device.malloc(mesh->Nelements*mesh->NpFEM*sizeof(dfloat));
  precon->o_zFEM = mesh->device.malloc(mesh->Nelements*mesh->NpFEM*sizeof(dfloat));

  precon->o_GrFEM = mesh->device.malloc(precon->hgs->Ngather*sizeof(dfloat));
  precon->o_GzFEM = mesh->device.malloc(precon->hgs->Ngather*sizeof(dfloat));

  free(FEMverts);
  free(faceMap);
  free(globalOwners);
  free(globalNumbering);
}