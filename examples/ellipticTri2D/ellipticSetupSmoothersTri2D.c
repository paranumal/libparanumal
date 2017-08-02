#include "ellipticTri2D.h"

typedef struct{

  iint row;
  iint col;
  iint ownerRank;
  dfloat val;

}nonZero_t;

typedef struct{

  iint localId;
  iint baseId;
  iint haloFlag;

} preconGatherInfo_t;

int parallelCompareBaseId(const void *a, const void *b){

  preconGatherInfo_t *fa = (preconGatherInfo_t*) a;
  preconGatherInfo_t *fb = (preconGatherInfo_t*) b;

  if(fa->baseId < fb->baseId) return -1;
  if(fa->baseId > fb->baseId) return +1;

  return 0;
}

void ellipticBuildExactPatchesIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, iint *BCType,
                                   dfloat **patchesInvA, const char *options);

void ellipticBuildApproxPatchesIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, iint *BCType,
                                   iint *Npataches, iint **patchesIndex, dfloat **patchesInvA,
                                   const char *options);

void ellipticBuildJacobiIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda,
                                   iint *BCType, dfloat **invDiagA,
                                   const char *options);


void ellipticSetupSmootherTri2D(solver_t *solver, precon_t *precon,
                                dfloat tau, dfloat lambda, int* BCType,
                                const char *options) {

  
  
}


void ellipticSetupSmootherOverlappingPatchIpdg(solver_t *solver, precon_t *precon,
                                              dfloat tau, dfloat lambda, int* BCType,
                                              dfloat weight, const char *options) {

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D *mesh = solver->mesh;

  iint Nlocal = mesh->Np*mesh->Nelements;
  iint Nhalo  = mesh->Np*mesh->totalHaloPairs;
  iint Ntrace = mesh->Nfp*mesh->Nfaces*mesh->Nelements;

  // build gather-scatter
  iint NpP = mesh->Np + mesh->Nfaces*mesh->Nfp;

  // build gather-scatter for overlapping patches
  iint *allNelements = (iint*) calloc(size, sizeof(iint));
  MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT,
              allNelements, 1, MPI_IINT, MPI_COMM_WORLD);

  // offsets
  iint *startElement = (iint*) calloc(size, sizeof(iint));
  for(iint r=1;r<size;++r){
    startElement[r] = startElement[r-1]+allNelements[r-1];
  }

  // 1-indexed numbering of nodes on this process
  iint *localNums = (iint*) calloc((Nlocal+Nhalo), sizeof(iint));
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      localNums[e*mesh->Np+n] = 1 + e*mesh->Np + n + startElement[rank]*mesh->Np;
    }
  }

  if(Nhalo){
    // send buffer for outgoing halo
    iint *sendBuffer = (iint*) calloc(Nhalo, sizeof(iint));

    // exchange node numbers with neighbors
    meshHaloExchange(mesh,
                   mesh->Np*sizeof(iint),
                   localNums,
                   sendBuffer,
                   localNums+Nlocal);
  }

  preconGatherInfo_t *preconGatherInfoDg =
    (preconGatherInfo_t*) calloc(NpP*mesh->Nelements,sizeof(preconGatherInfo_t));

  // set local ids
  for(iint n=0;n<mesh->Nelements*NpP;++n)
    preconGatherInfoDg[n].localId = n;

  // numbering of patch interior nodes
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      iint id  = n + e*mesh->Np;
      iint pid = n + e*NpP;

      // all patch interior nodes are local
      preconGatherInfoDg[pid].baseId = localNums[id];
    }
  }

  // add patch boundary nodes
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint f=0;f<mesh->Nfaces;++f){
      // mark halo nodes
      iint rP = mesh->EToP[e*mesh->Nfaces+f];
      iint eP = mesh->EToE[e*mesh->Nfaces+f];
      iint fP = mesh->EToF[e*mesh->Nfaces+f];
      iint bc = mesh->EToB[e*mesh->Nfaces+f];

      for(iint n=0;n<mesh->Nfp;++n){
      iint id = n + f*mesh->Nfp+e*mesh->Nfp*mesh->Nfaces;
      iint idP = mesh->vmapP[id];

      // local numbers
      iint pidM = e*NpP + mesh->faceNodes[f*mesh->Nfp+n];
      iint pidP = e*NpP + mesh->Np + f*mesh->Nfp+n;
      preconGatherInfoDg[pidP].baseId = localNums[idP];

      if(rP!=-1){
        preconGatherInfoDg[pidM].haloFlag = 1;
        preconGatherInfoDg[pidP].haloFlag = 1;
      }
      }
    }
  }

  // sort by rank then base index
  qsort(preconGatherInfoDg, NpP*mesh->Nelements, sizeof(preconGatherInfo_t),
      parallelCompareBaseId);

  // do not gather-scatter nodes labelled zero
  int skip = 0;

  while(preconGatherInfoDg[skip].baseId==0 && skip<NpP*mesh->Nelements){
    ++skip;
  }

  // reset local ids
  iint NlocalDg = NpP*mesh->Nelements - skip;
  iint *gatherLocalIdsDg  = (iint*) calloc(NlocalDg, sizeof(iint));
  iint *gatherBaseIdsDg   = (iint*) calloc(NlocalDg, sizeof(iint));
  iint *gatherHaloFlagsDg = (iint*) calloc(NlocalDg, sizeof(iint));
  for(iint n=0;n<NlocalDg;++n){
    gatherLocalIdsDg[n]  = preconGatherInfoDg[n+skip].localId;
    gatherBaseIdsDg[n]   = preconGatherInfoDg[n+skip].baseId;
    gatherHaloFlagsDg[n] = preconGatherInfoDg[n+skip].haloFlag;
  }

  // make preconBaseIds => preconNumbering
  precon->ogsDg = meshParallelGatherScatterSetup(mesh,
                                               NlocalDg,
                                               sizeof(dfloat),
                                               gatherLocalIdsDg,
                                               gatherBaseIdsDg,
                                               gatherHaloFlagsDg);

  //correction for full patch
  NpP = mesh->NpP;

  // build degree vector
  iint NtotalDGP = NpP*mesh->Nelements;
  dfloat *invDegree = (dfloat*) calloc(NtotalDGP, sizeof(dfloat));
  dfloat *degree    = (dfloat*) calloc(NtotalDGP, sizeof(dfloat));
  precon->o_invDegreeDGP = mesh->device.malloc(NtotalDGP*sizeof(dfloat), invDegree);

  for(iint n=0;n<NtotalDGP;++n)
    degree[n] = 1;

  occa::memory o_deg = mesh->device.malloc(NtotalDGP*sizeof(dfloat), degree);
  meshParallelGatherScatter(mesh, precon->ogsDg, o_deg, o_deg, dfloatString, "add");
  o_deg.copyTo(degree);
  mesh->device.finish();
  o_deg.free();

  for(iint n=0;n<NtotalDGP;++n){ // need to weight inner products{
    if(degree[n] == 0) printf("WARNING!!!!\n");
    invDegree[n] = 1./degree[n];
  }

  precon->o_invDegreeDGP.copyFrom(invDegree);
  free(degree);
  free(invDegree);

  // -------------------------------------------------------------------------------------------
  // load prebuilt transform matrices
  precon->o_oasForwardDg = mesh->device.malloc(NpP*NpP*sizeof(dfloat), mesh->oasForwardDg);
  precon->o_oasBackDg    = mesh->device.malloc(NpP*NpP*sizeof(dfloat), mesh->oasBackDg);

  dfloat *oasForwardDgT = (dfloat*) calloc(NpP*NpP, sizeof(dfloat));
  dfloat *oasBackDgT = (dfloat*) calloc(NpP*NpP, sizeof(dfloat));
  for(iint n=0;n<NpP;++n){
    for(iint m=0;m<NpP;++m){
      oasForwardDgT[n+m*NpP] = mesh->oasForwardDg[m+n*NpP];
      oasBackDgT[n+m*NpP] = mesh->oasBackDg[m+n*NpP];
    }
  }

  precon->o_oasForwardDgT = mesh->device.malloc(NpP*NpP*sizeof(dfloat), oasForwardDgT);
  precon->o_oasBackDgT    = mesh->device.malloc(NpP*NpP*sizeof(dfloat), oasBackDgT);


  /// ---------------------------------------------------------------------------

  // hack estimate for Jacobian scaling

  dfloat *diagInvOp = (dfloat*) calloc(NpP*mesh->Nelements, sizeof(dfloat));
  dfloat *diagInvOpDg = (dfloat*) calloc(NpP*mesh->Nelements, sizeof(dfloat));
  for(iint e=0;e<mesh->Nelements;++e){

    // S = Jabc*(wa*wb*wc*lambda + wb*wc*Da'*wa*Da + wa*wc*Db'*wb*Db + wa*wb*Dc'*wc*Dc)
    // S = Jabc*wa*wb*wc*(lambda*I+1/wa*Da'*wa*Da + 1/wb*Db'*wb*Db + 1/wc*Dc'*wc*Dc)

    dfloat J = mesh->vgeo[e*mesh->Nvgeo + JID];
    dfloat rx = mesh->vgeo[e*mesh->Nvgeo + RXID];
    dfloat sx = mesh->vgeo[e*mesh->Nvgeo + SXID];
    dfloat ry = mesh->vgeo[e*mesh->Nvgeo + RYID];
    dfloat sy = mesh->vgeo[e*mesh->Nvgeo + SYID];

    //metric tensor on this element
    dfloat grr = rx*rx+ry*ry;
    dfloat grs = rx*sx+ry*sy;
    dfloat gss = sx*sx+sy*sy;

    //eigenvalues of the metric
    dfloat eigG1 = 0.5*(grr + gss) - 0.5*sqrt((grr-gss)*(grr-gss) + 4*grs*grs);
    dfloat eigG2 = 0.5*(grr + gss) + 0.5*sqrt((grr-gss)*(grr-gss) + 4*grs*grs);

    //hinv scaling?
    dfloat invhMin = 1e9;
    dfloat invhMax = 0;
    for (iint f=0;f<mesh->Nfaces;f++) {
      iint base = e*mesh->Nfaces + f;
      //invh = mymax(mesh->sgeo[base*mesh->Nsgeo+IHID],invh);
      invhMax =  mymax(mesh->sgeo[base*mesh->Nsgeo + SJID],invhMax);
      invhMin =  mymin(mesh->sgeo[base*mesh->Nsgeo + SJID],invhMin);
    }
    dfloat invh = invhMax/invhMin +invhMin/invhMax;

    //TODO Average? Min/Max/Avg Eigenvalue? What works best for the scaling?
    //dfloat Jhinv2 = J*(eigG1+eigG2);
    //dfloat Jhinv2 = (2*mesh->N+1)*J*mymax(eigG1,eigG2);
    dfloat Jhinv2 = J*mymin(eigG1,eigG2);
    //dfloat Jhinv2 = mesh->N*invh;

    for(iint n=0;n<NpP;++n){
      iint pid = n + e*NpP;

      diagInvOpDg[pid] = weight/(J*lambda + Jhinv2*mesh->oasDiagOpDg[n]);
    }
  }

  precon->o_oasDiagInvOpDg =
    mesh->device.malloc(NpP*mesh->Nelements*sizeof(dfloat), diagInvOpDg);

  //storage buffer for patches
  iint NtotalP = mesh->NpP*mesh->Nelements;
  precon->zP  = (dfloat*) calloc(NtotalP,  sizeof(dfloat));
  precon->o_zP  = mesh->device.malloc(NtotalP*sizeof(dfloat),precon->zP);
}

void ellipticSetupSmootherExactFullPatchIpdg(solver_t *solver, precon_t *precon,
                                            dfloat tau, dfloat lambda, int* BCType,
                                            dfloat weight, const char *options) {

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  dfloat *invAP;
  iint Npatches;
  mesh_t *mesh = solver->mesh;

  int NpP = mesh->Np*(mesh->Nfaces+1);

  //initialize the full inverse operators on each 4 element patch
  ellipticBuildExactPatchesIpdgTri2D(mesh, mesh->Np, NULL, tau, lambda, BCType, &invAP, options);

  precon->o_invAP = mesh->device.malloc(mesh->Nelements*NpP*NpP*sizeof(dfloat),invAP);

  dfloat *invDegree = (dfloat*) calloc(mesh->Nelements,sizeof(dfloat));
  for (iint e=0;e<mesh->Nelements;e++) {
    for (int f=0;f<mesh->Nfaces;f++)
        invDegree[e] += (mesh->EToE[e*mesh->Nfaces +f]<0) ? 0 : 1; //overlap degree = # of neighbours
    invDegree[e] = weight/invDegree[e]; //build in weight
  }
  precon->o_invDegreeAP = mesh->device.malloc(mesh->Nelements*sizeof(dfloat),invDegree);
  free(invDegree);

  mesh->o_EToE  = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),mesh->EToE);
  mesh->o_EToF  = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),mesh->EToF);
  mesh->o_rmapP = mesh->device.malloc(mesh->Np*mesh->Nfaces*sizeof(iint),mesh->rmapP);

  //set storage for larger patch
  precon->zP = (dfloat*) calloc(mesh->Nelements*NpP,  sizeof(dfloat));
  precon->o_zP = mesh->device.malloc(mesh->Nelements*NpP*sizeof(dfloat), precon->zP);
}

void ellipticSetupSmootherApproxFullPatchIpdg(solver_t *solver, precon_t *precon,
                                            dfloat tau, dfloat lambda, int* BCType,
                                            dfloat weight, const char *options) {

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  dfloat *invAP;
  iint Npatches;
  iint *patchesIndex;
  mesh_t *mesh = solver->mesh;

  int NpP = mesh->Np*(mesh->Nfaces+1);

  //initialize the full inverse operators on each 4 element patch
  ellipticBuildApproxPatchesIpdgTri2D(mesh, mesh->Np, NULL, tau, lambda, BCType,
                                      &Npatches, &patchesIndex, &invAP, options);

  precon->o_invAP = mesh->device.malloc(Npatches*NpP*NpP*sizeof(dfloat),invAP);
  precon->o_patchesIndex = mesh->device.malloc(mesh->Nelements*sizeof(iint), patchesIndex);

  dfloat *invDegree = (dfloat*) calloc(mesh->Nelements,sizeof(dfloat));
  for (iint e=0;e<mesh->Nelements;e++) {
    for (int f=0;f<mesh->Nfaces;f++)
        invDegree[e] += (mesh->EToE[e*mesh->Nfaces +f]<0) ? 0 : 1; //overlap degree = # of neighbours
    invDegree[e] = weight/invDegree[e]; //build in weight
  }
  precon->o_invDegreeAP = mesh->device.malloc(mesh->Nelements*sizeof(dfloat),invDegree);
  free(invDegree);

  mesh->o_EToE = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),mesh->EToE);
  mesh->o_EToF = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),mesh->EToF);

  //set storage for larger patch
  precon->zP = (dfloat*) calloc(mesh->Nelements*NpP,  sizeof(dfloat));
  precon->o_zP = mesh->device.malloc(mesh->Nelements*NpP*sizeof(dfloat), precon->zP);

}

void ellipticSetupSmootherDampedJacobiIpdg(solver_t *solver, precon_t *precon,
                                            dfloat tau, dfloat lambda, int* BCType,
                                            dfloat weight, const char *options) {

  dfloat *invDiagA;
  mesh_t *mesh = solver->mesh;

  ellipticBuildJacobiIpdgTri2D(mesh,mesh->Np,NULL,tau, lambda, BCType, &invDiagA,options);

  dfloat weight = 0.5; //dampening factor (may need to be tuned)
  for (iint n=0;n<mesh->Np*mesh->Nelements;n++)
    invDiagA[n] *= weight;

  precon->o_invDiagA = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), invDiagA);
}