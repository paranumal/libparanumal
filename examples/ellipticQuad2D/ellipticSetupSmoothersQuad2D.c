#include "ellipticQuad2D.h"

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

void ellipticSetupSmootherOverlappingPatchIpdg(solver_t *solver, precon_t *precon, agmgLevel *level,
                                              dfloat tau, dfloat lambda, int* BCType, const char *options) {

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D *mesh = solver->mesh;
  precon_t *precon= solver->precon;

  // offsets to extract second layer
  iint *offset = (iint*) calloc(mesh->Nfaces, sizeof(iint));
  offset[0] = +mesh->Nq;
  offset[1] = -1;
  offset[2] = -mesh->Nq;
  offset[3] = +1;

  // build gather-scatter
  iint NqP = mesh->Nq+2;
  iint NpP = NqP*NqP;

  iint *offsetP = (iint*) calloc(mesh->Nfaces, sizeof(iint));
  offsetP[0] = +NqP;
  offsetP[1] = -1;
  offsetP[2] = -NqP;
  offsetP[3] = +1;

  iint *faceNodesPrecon = (iint*) calloc(mesh->Nfp*mesh->Nfaces, sizeof(iint));

  for(iint i=0;i<mesh->Nq;++i) faceNodesPrecon[i+0*mesh->Nfp] = i+1 + 0*NqP;
  for(iint j=0;j<mesh->Nq;++j) faceNodesPrecon[j+1*mesh->Nfp] = NqP-1 + (j+1)*NqP;
  for(iint i=0;i<mesh->Nq;++i) faceNodesPrecon[i+2*mesh->Nfp] = i+1 + (NqP-1)*NqP;
  for(iint j=0;j<mesh->Nq;++j) faceNodesPrecon[j+3*mesh->Nfp] = 0 + (j+1)*NqP;

  iint *vmapMP = (iint*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(iint));
  iint *vmapPP = (iint*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(iint));

  // take node info from positive trace and put in overlap region on parallel gather info
  for(iint e=0;e<mesh->Nelements;++e){

    for(iint n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      iint fid = e*mesh->Nfp*mesh->Nfaces + n;

      // find face index
      iint fM = n/mesh->Nfp;
      iint fP = mesh->EToF[e*mesh->Nfaces+fM];
      if(fP<0) fP = fM;

      // find location of "+" trace in regular mesh
      iint idM = mesh->vmapM[fid] + offset[fM];
      iint idP = mesh->vmapP[fid] + offset[fP];

      vmapMP[fid] = idM;
      vmapPP[fid] = idP;
    }
  }

  // -------------------------------------------------------------------------------------------

  // space for gather base indices with halo
  preconGatherInfo_t *gatherInfo =
    (preconGatherInfo_t*) calloc(Nlocal+Nhalo, sizeof(preconGatherInfo_t));

  // rearrange in node order
  for(iint n=0;n<Nlocal;++n){
    iint id = mesh->gatherLocalIds[n];
    gatherInfo[id].baseId   = mesh->gatherBaseIds[n] + 1;
    gatherInfo[id].haloFlag = mesh->gatherHaloFlags[n];
  }

  if(Nhalo){
    // send buffer for outgoing halo
    preconGatherInfo_t *sendBuffer = (preconGatherInfo_t*) calloc(Nhalo, sizeof(preconGatherInfo_t));

    meshHaloExchange(mesh,
         mesh->Np*sizeof(preconGatherInfo_t),
         gatherInfo,
         sendBuffer,
         gatherInfo+Nlocal);
  }
  // now create padded version

  // now find info about halo nodes
  preconGatherInfo_t *preconGatherInfo =
    (preconGatherInfo_t*) calloc(NpP*mesh->Nelements, sizeof(preconGatherInfo_t));

  // push to non-overlap nodes
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint j=0;j<mesh->Nq;++j){
      for(iint i=0;i<mesh->Nq;++i){
        iint id  = i + j*mesh->Nq + e*mesh->Np;
        iint pid = i + 1 + (j+1)*NqP + e*NpP;

        preconGatherInfo[pid] = gatherInfo[id];
      }
    }
  }

  // add overlap region
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      iint id = n + e*mesh->Nfp*mesh->Nfaces;

      iint f = n/mesh->Nfp;
      iint idM = mesh->vmapM[id];
      iint idP = vmapPP[id];
      iint idMP = e*NpP + faceNodesPrecon[n];

      preconGatherInfo[idMP] = gatherInfo[idP];

#if 1
      if(gatherInfo[idM].haloFlag){
        preconGatherInfo[idMP].haloFlag = 1;
        preconGatherInfo[idMP+offsetP[f]].haloFlag = 1;
        preconGatherInfo[idMP+2*offsetP[f]].haloFlag = 1;
      }
#endif
    }
  }

  // reset local ids
  for(iint n=0;n<mesh->Nelements*NpP;++n)
    preconGatherInfo[n].localId = n;

  // sort by rank then base index
  qsort(preconGatherInfo, NpP*mesh->Nelements, sizeof(preconGatherInfo_t), parallelCompareBaseId);

  // do not gather-scatter nodes labelled zero
  iint skip = 0;
  while(preconGatherInfo[skip].baseId==0 && skip<NpP*mesh->Nelements){
    ++skip;
  }

  // reset local ids
  iint NlocalP = NpP*mesh->Nelements - skip;
  iint *gatherLocalIdsP  = (iint*) calloc(NlocalP, sizeof(iint));
  iint *gatherBaseIdsP   = (iint*) calloc(NlocalP, sizeof(iint));
  iint *gatherHaloFlagsP = (iint*) calloc(NlocalP, sizeof(iint));
  for(iint n=0;n<NlocalP;++n){
    gatherLocalIdsP[n]  = preconGatherInfo[n+skip].localId;
    gatherBaseIdsP[n]   = preconGatherInfo[n+skip].baseId;
    gatherHaloFlagsP[n] = preconGatherInfo[n+skip].haloFlag;
  }

  // make preconBaseIds => preconNumbering
  precon->ogsP = meshParallelGatherScatterSetup(mesh,
            NlocalP,
            sizeof(dfloat),
            gatherLocalIdsP,
            gatherBaseIdsP,
            gatherHaloFlagsP);

  // build degree vector
  iint NtotalP = mesh->NqP*mesh->NqP*mesh->Nelements;
  dfloat *invDegree = (dfloat*) calloc(NtotalP, sizeof(dfloat));
  dfloat *degree    = (dfloat*) calloc(NtotalP, sizeof(dfloat));
  precon->o_invDegreeP = mesh->device.malloc(NtotalP*sizeof(dfloat), invDegree);

  for(iint n=0;n<NtotalP;++n)
    degree[n] = 1;

  occa::memory o_deg = mesh->device.malloc(NtotalP*sizeof(dfloat), degree);
  meshParallelGatherScatter(mesh, precon->ogsP, o_deg, o_deg, dfloatString, "add");
  o_deg.copyTo(degree);
  mesh->device.finish();
  o_deg.free();

  for(iint n=0;n<NtotalP;++n){ // need to weight inner products{
    if(degree[n] == 0) printf("WARNING!!!!\n");
    invDegree[n] = 1./degree[n];
  }

  precon->o_invDegreeP.copyFrom(invDegree);
  free(degree);
  free(invDegree);

  // -------------------------------------------------------------------------------------------
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
    (preconGatherInfo_t*) calloc(NpP*mesh->Nelements,
         sizeof(preconGatherInfo_t));

  // set local ids
  for(iint n=0;n<mesh->Nelements*NpP;++n)
    preconGatherInfoDg[n].localId = n;

  // numbering of patch interior nodes
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint j=0;j<mesh->Nq;++j){
      for(iint i=0;i<mesh->Nq;++i){
        iint id  = i + j*mesh->Nq + e*mesh->Np;
        iint pid = (i+1) + (j+1)*NqP + e*NpP;

        // all patch interior nodes are local
        preconGatherInfoDg[pid].baseId = localNums[id];
      }
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
        iint pidM = e*NpP + faceNodesPrecon[f*mesh->Nfp+n] + offsetP[f];
        iint pidP = e*NpP + faceNodesPrecon[f*mesh->Nfp+n];
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
  skip = 0;

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

  // build degree vector
  iint NtotalDGP = NpP*mesh->Nelements;
  invDegree = (dfloat*) calloc(NtotalDGP, sizeof(dfloat));
  degree    = (dfloat*) calloc(NtotalDGP, sizeof(dfloat));
  precon->o_invDegreeDGP = mesh->device.malloc(NtotalDGP*sizeof(dfloat), invDegree);

  for(iint n=0;n<NtotalDGP;++n)
    degree[n] = 1;

  o_deg = mesh->device.malloc(NtotalDGP*sizeof(dfloat), degree);
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

  precon->o_faceNodesP = mesh->device.malloc(mesh->Nfp*mesh->Nfaces*sizeof(iint), faceNodesPrecon);
  precon->o_vmapPP     = mesh->device.malloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements*sizeof(iint), vmapPP);

  // load prebuilt transform matrices
  precon->o_oasForward = mesh->device.malloc(NqP*NqP*sizeof(dfloat), mesh->oasForward);
  precon->o_oasBack    = mesh->device.malloc(NqP*NqP*sizeof(dfloat), mesh->oasBack);

  precon->o_oasForwardDg = mesh->device.malloc(NqP*NqP*sizeof(dfloat), mesh->oasForwardDg);
  precon->o_oasBackDg    = mesh->device.malloc(NqP*NqP*sizeof(dfloat), mesh->oasBackDg);

  /// ---------------------------------------------------------------------------

  // hack estimate for Jacobian scaling

  dfloat *diagInvOp = (dfloat*) calloc(NpP*mesh->Nelements, sizeof(dfloat));
  dfloat *diagInvOpDg = (dfloat*) calloc(NpP*mesh->Nelements, sizeof(dfloat));
  for(iint e=0;e<mesh->Nelements;++e){

    // S = Jabc*(wa*wb*wc*lambda + wb*wc*Da'*wa*Da + wa*wc*Db'*wb*Db + wa*wb*Dc'*wc*Dc)
    // S = Jabc*wa*wb*wc*(lambda*I+1/wa*Da'*wa*Da + 1/wb*Db'*wb*Db + 1/wc*Dc'*wc*Dc)

    dfloat Jhrinv2 = 0, Jhsinv2 = 0, J = 0;
    for(iint n=0;n<mesh->Np;++n){

      dfloat W = mesh->gllw[n%mesh->Nq]*mesh->gllw[n/mesh->Nq];

      iint base = mesh->Nggeo*mesh->Np*e + n;

      J = mymax(J, mesh->ggeo[base + mesh->Np*GWJID]/W);
      Jhrinv2 = mymax(Jhrinv2, mesh->ggeo[base + mesh->Np*G00ID]/W);
      Jhsinv2 = mymax(Jhsinv2, mesh->ggeo[base + mesh->Np*G11ID]/W);
    }

    for(iint j=0;j<NqP;++j){
      for(iint i=0;i<NqP;++i){
        iint pid = i + j*NqP + e*NpP;

          diagInvOp[pid] =
            1./(J*lambda +
          Jhrinv2*mesh->oasDiagOp[i] +
          Jhsinv2*mesh->oasDiagOp[j]);

          diagInvOpDg[pid] =
            1./(J*lambda +
          Jhrinv2*mesh->oasDiagOpDg[i] +
          Jhsinv2*mesh->oasDiagOpDg[j]);

      }
    }
  }

  precon->o_oasDiagInvOp = mesh->device.malloc(NpP*mesh->Nelements*sizeof(dfloat), diagInvOp);
  precon->o_oasDiagInvOpDg = mesh->device.malloc(NpP*mesh->Nelements*sizeof(dfloat), diagInvOpDg);

  /// ---------------------------------------------------------------------------
  // compute diagonal of stiffness matrix for Jacobi

  iint Ntotal = mesh->Np*mesh->Nelements;
  dfloat *diagA = (dfloat*) calloc(Ntotal, sizeof(dfloat));

  for(iint e=0;e<mesh->Nelements;++e){
    iint cnt = 0;
    for(iint j=0;j<mesh->Nq;++j){
      for(iint i=0;i<mesh->Nq;++i){

        dfloat JW = mesh->ggeo[e*mesh->Np*mesh->Nggeo+ cnt + mesh->Np*GWJID];
        // (D_{ii}^2 + D_{jj}^2 + lambda)*w_i*w_j*w_k*J_{ijke}
        diagA[e*mesh->Np+cnt] = (pow(mesh->D[i*mesh->Nq+i],2) +
               pow(mesh->D[j*mesh->Nq+j],2) +
               lambda)*JW;
        ++cnt;
      }
    }
  }

  precon->o_diagA = mesh->device.malloc(Ntotal*sizeof(dfloat), diagA);

  // sum up
  meshParallelGatherScatter(mesh, ogs, precon->o_diagA, precon->o_diagA, dfloatString, "add");

  
  //storage buffer for patches
  iint NtotalP = mesh->NpP*mesh->Nelements;
  precon->zP  = (dfloat*) calloc(NtotalP,  sizeof(dfloat));
  precon->o_zP  = mesh->device.malloc(NtotalP*sizeof(dfloat),precon->zP);


  level->device_smoother = overlappingPatchIpdg;

  //check if stabilization is needed
  if (strstr(options,"MULTIGRID")||strstr(options,"FULLALMOND")) {
    //estimate the max eigenvalue of S*A
    dfloat rho = maxEigSmoothAx(solver, level);

    if (strstr(options,"CHEBYSHEV")) {

      level->smoother_params = (dfloat *) calloc(2,sizeof(dfloat));

      level->ChebyshevIterations = 2;
      level->smoother_params[0] = rho;
      level->smoother_params[1] = rho/10.;

    } else {

      //set the stabilty weight (jacobi-type interation)
      dfloat weight = (4./3.)/rho;

      printf("weight = %g \n", weight);

      for (iint n=0;n<NpP*mesh->Nelements;n++)
        diagInvOpDg[n] *= weight;

      //update diagonal with weight
      precon->o_oasDiagInvOpDg.copyFrom(diagInvOpDg);
    }
  }
}

void ellipticSetupSmootherFullPatchIpdg(solver_t *solver, precon_t *precon, agmgLevel *level,
                                              dfloat tau, dfloat lambda, int* BCType, dfloat rateTolerance, const char *options) {

  dfloat *invAP;
  iint Npatches;
  iint *patchesIndex;
  mesh_t *mesh = solver->mesh;

  int NpP = mesh->Np*(mesh->Nfaces+1);

  //initialize the full inverse operators on each 4 element patch
  ellipticBuildFullPatchesIpdgQuad2D(mesh, mesh->Np, NULL, tau, lambda, BCType, rateTolerance,
                                      &Npatches, &patchesIndex, &invAP, options);

  precon->o_invAP = mesh->device.malloc(Npatches*NpP*NpP*sizeof(dfloat),invAP);
  precon->o_patchesIndex = mesh->device.malloc(mesh->Nelements*sizeof(iint), patchesIndex);

  dfloat *invDegree = (dfloat*) calloc(mesh->Nelements,sizeof(dfloat));
  for (iint e=0;e<mesh->Nelements;e++) {
    for (int f=0;f<mesh->Nfaces;f++)
        invDegree[e] += (mesh->EToE[e*mesh->Nfaces +f]<0) ? 0 : 1; //overlap degree = # of neighbours
    invDegree[e] = 1.0/invDegree[e]; //build in weight
  }
  precon->o_invDegreeAP = mesh->device.malloc(mesh->Nelements*sizeof(dfloat),invDegree);

  mesh->o_EToE = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),mesh->EToE);
  mesh->o_EToF = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),mesh->EToF);

  //set storage for larger patch
  precon->zP = (dfloat*) calloc(mesh->Nelements*NpP,  sizeof(dfloat));
  precon->o_zP = mesh->device.malloc(mesh->Nelements*NpP*sizeof(dfloat), precon->zP);

  level->device_smoother = FullPatchIpdg;

  //check if stabilization is needed
  if (strstr(options,"MULTIGRID")||strstr(options,"FULLALMOND")) {
    //estimate the max eigenvalue of S*A
    dfloat rho = maxEigSmoothAx(solver, level);

    if (strstr(options,"CHEBYSHEV")) {

      level->smoother_params = (dfloat *) calloc(2,sizeof(dfloat));

      level->ChebyshevIterations = 2;
      level->smoother_params[0] = rho;
      level->smoother_params[1] = rho/10.;

    } else {

      //set the stabilty weight (jacobi-type interation)
      dfloat weight = (4./3.)/rho;

      printf("weight = %g \n", weight);

      for (iint e=0;e<mesh->Nelements;e++)
        invDegree[e] *= weight;

      //update with weight
      precon->o_invDegreeAP.copyFrom(invDegree);
    }
  }
  free(invDegree);
}

void ellipticSetupSmootherFacePatchIpdg(solver_t *solver, precon_t *precon, agmgLevel *level,
                                              dfloat tau, dfloat lambda, int* BCType, dfloat rateTolerance, const char *options) {

  dfloat *invAP;
  iint Npatches;
  iint *patchesIndex;
  mesh_t *mesh = solver->mesh;

  //initialize the full inverse operators on each 4 element patch
  ellipticBuildFacePatchesIpdgQuad2D(mesh, mesh->Np, NULL, tau, lambda, BCType, rateTolerance,
                                      &Npatches, &patchesIndex, &invAP, options);

  int NpP = 2*mesh->Np;

  precon->o_invAP = mesh->device.malloc(Npatches*NpP*NpP*sizeof(dfloat),invAP);
  precon->o_patchesIndex = mesh->device.malloc(mesh->NfacePairs*sizeof(iint), patchesIndex);

  dfloat *invDegree = (dfloat*) calloc(mesh->Nelements+mesh->totalHaloPairs,sizeof(dfloat));
  for (iint face=0;face<mesh->NfacePairs;face++) {
    iint eM = mesh->FPairsToE[2*face+0];
    iint eP = mesh->FPairsToE[2*face+1];

    invDegree[eM]++; //overlap degree = # of patches
    if (eP>=0) invDegree[eP]++; //overlap degree = # of patches
  }
  for (iint e=0;e<mesh->Nelements+mesh->totalHaloPairs;e++) {
    invDegree[e] = 1.0/invDegree[e];
  }

  precon->o_invDegreeAP = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*sizeof(dfloat),invDegree);

  mesh->o_FPairsToE = mesh->device.malloc(2*mesh->NfacePairs*sizeof(iint),mesh->FPairsToE);
  mesh->o_FPairsToF = mesh->device.malloc(2*mesh->NfacePairs*sizeof(iint),mesh->FPairsToF);
  mesh->o_EToFPairs = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),mesh->EToFPairs);

  //set storage for larger patch
  precon->zP = (dfloat*) calloc(mesh->NfacePairs*NpP,  sizeof(dfloat));
  precon->o_zP = mesh->device.malloc(mesh->NfacePairs*NpP*sizeof(dfloat), precon->zP);


  level->device_smoother = FacePatchIpdg;

  //check if stabilization is needed
  if (strstr(options,"MULTIGRID")||strstr(options,"FULLALMOND")) {
    //estimate the max eigenvalue of S*A
    dfloat rho = maxEigSmoothAx(solver, level);

    if (strstr(options,"CHEBYSHEV")) {

      level->smoother_params = (dfloat *) calloc(2,sizeof(dfloat));

      level->ChebyshevIterations = 2;
      level->smoother_params[0] = rho;
      level->smoother_params[1] = rho/10.;

    } else {

      //set the stabilty weight (jacobi-type interation)
      dfloat weight = (4./3.)/rho;

      printf("weight = %g \n", weight);

      for (iint e=0;e<mesh->Nelements;e++)
        invDegree[e] *= weight;

      //update with weight
      precon->o_invDegreeAP.copyFrom(invDegree);
    }
  }
  free(invDegree);
}

void ellipticSetupSmootherLocalPatchIpdg(solver_t *solver, precon_t *precon, agmgLevel *level,
                                              dfloat tau, dfloat lambda, int* BCType, dfloat rateTolerance, const char *options) {

  dfloat *invAP;
  iint Npatches;
  iint *patchesIndex;
  mesh_t *mesh = solver->mesh;

  int NpP = mesh->Np;

  //initialize the full inverse operators on each 4 element patch
  ellipticBuildLocalPatchesIpdgQuad2D(mesh, mesh->Np, NULL, tau, lambda, BCType, rateTolerance,
                                      &Npatches, &patchesIndex, &invAP, options);

  precon->o_invAP = mesh->device.malloc(Npatches*NpP*NpP*sizeof(dfloat),invAP);
  precon->o_patchesIndex = mesh->device.malloc(mesh->Nelements*sizeof(iint), patchesIndex);

  dfloat *invDegree = (dfloat*) calloc(mesh->Nelements,sizeof(dfloat));
  for (iint e=0;e<mesh->Nelements;e++) {
    invDegree[e] = 1.0;
  }
  precon->o_invDegreeAP = mesh->device.malloc(mesh->Nelements*sizeof(dfloat),invDegree);

  level->device_smoother = LocalPatchIpdg;

  //check if stabilization is needed
  if (strstr(options,"MULTIGRID")||strstr(options,"FULLALMOND")) {
    //estimate the max eigenvalue of S*A
    dfloat rho = maxEigSmoothAx(solver, level);

    if (strstr(options,"CHEBYSHEV")) {

      level->smoother_params = (dfloat *) calloc(2,sizeof(dfloat));

      level->ChebyshevIterations = 2;
      level->smoother_params[0] = rho;
      level->smoother_params[1] = rho/10.;

    } else {

      //set the stabilty weight (jacobi-type interation)
      dfloat weight = (4./3.)/rho;

      printf("weight = %g \n", weight);

      for (iint e=0;e<mesh->Nelements;e++)
        invDegree[e] *= weight;

      //update with weight
      precon->o_invDegreeAP.copyFrom(invDegree);
    }
  }
  free(invDegree);
}

void ellipticSetupSmootherDampedJacobiIpdg(solver_t *solver, precon_t *precon, agmgLevel *level,
                                              dfloat tau, dfloat lambda, int* BCType, const char *options) {

  dfloat *invDiagA;
  mesh_t *mesh = solver->mesh;

  ellipticBuildJacobiIpdgQuad2D(mesh,mesh->Np,NULL,tau, lambda, BCType, &invDiagA,options);

  precon->o_invDiagA = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), invDiagA);

  level->device_smoother = dampedJacobi;

  //check if stabilization is needed
  if (strstr(options,"MULTIGRID")||strstr(options,"FULLALMOND")) {
    //estimate the max eigenvalue of S*A
    dfloat rho = maxEigSmoothAx(solver, level);

    if (strstr(options,"CHEBYSHEV")) {

      level->smoother_params = (dfloat *) calloc(2,sizeof(dfloat));

      level->ChebyshevIterations = 2;
      level->smoother_params[0] = rho;
      level->smoother_params[1] = rho/10.;

    } else {

      //set the stabilty weight (jacobi-type interation)
      dfloat weight = (4./3.)/rho;

      printf("weight = %g \n", weight);

      for (iint n=0;n<mesh->Np*mesh->Nelements;n++)
        invDiagA[n] *= weight;

      //update diagonal with weight
      precon->o_invDiagA.copyFrom(invDiagA);
    }
  }

  free(invDiagA);
}

static void eig(const int Nrows, double *A, double *WR, double *WI){

  int NB  = 256;
  char JOBVL  = 'V';
  char JOBVR  = 'V';
  int     N = Nrows;
  int   LDA = Nrows;
  int  LWORK  = (NB+2)*N;

  double *WORK  = new double[LWORK];
  double *VL  = new double[Nrows*Nrows];
  double *VR  = new double[Nrows*Nrows];

  int INFO = -999;

  dgeev_ (&JOBVL, &JOBVR, &N, A, &LDA, WR, WI,
    VL, &LDA, VR, &LDA, WORK, &LWORK, &INFO);


  assert(INFO == 0);

  delete [] VL;
  delete [] VR;
  delete [] WORK;
}

dfloat maxEigSmoothAx(solver_t* solver, agmgLevel *level){

  mesh_t *mesh = solver->mesh;

  const iint N = level->Nrows;
  const iint M = level->Ncols;

  int k = 10;

  iint rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iint Ntotal=0;
  MPI_Allreduce(&N, &Ntotal, 1, MPI_IINT, MPI_SUM, MPI_COMM_WORLD);
  if(k > Ntotal)
    k = Ntotal;

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  double *H = (double *) calloc(k*k,sizeof(double));

  // allocate memory for basis
  dfloat **V = (dfloat **) calloc(k+1, sizeof(dfloat *));
  dfloat *Vx = (dfloat *) calloc(M, sizeof(dfloat));

  occa::memory o_Vx  = mesh->device.malloc(M*sizeof(dfloat));
  occa::memory o_AVx = mesh->device.malloc(M*sizeof(dfloat));

  for(int i=0; i<=k; i++)
    V[i] = (dfloat *) calloc(N, sizeof(dfloat));

  // generate a random vector for initial basis vector
  for (iint i=0;i<N;i++)
    Vx[i] = (dfloat) drand48();

  dfloat norm_vo = 0.;
  for (iint i=0;i<N;i++)
    norm_vo += Vx[i]*Vx[i];

  dfloat gNorm_vo = 0;
  MPI_Allreduce(&norm_vo, &gNorm_vo, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
  gNorm_vo = sqrt(gNorm_vo);

  for (iint i=0;i<N;i++)
    Vx[i] /= gNorm_vo;

  for (iint i=0;i<N;i++)
    V[0][i] = Vx[i];

  for(int j=0; j<k; j++){

    for (iint i=0;i<N;i++)
      Vx[i] = V[j][i];

    o_Vx.copyFrom(Vx); //send to device

    // v[j+1] = invD*(A*v[j])
    level->device_Ax(level->AxArgs,o_Vx,o_AVx);
    level->device_smoother(level->smootherArgs, o_AVx, o_AVx);

    o_AVx.copyTo(V[j+1],N*sizeof(dfloat)); //send result to host

    // modified Gram-Schmidth
    for(int i=0; i<=j; i++){
      // H(i,j) = v[i]'*A*v[j]
      dfloat hij = 0.;
      for (iint n=0;n<N;n++)
        hij += V[i][n]*V[j+1][n];

      dfloat ghij = 0;
      MPI_Allreduce(&hij, &ghij, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

      // v[j+1] = v[j+1] - hij*v[i]
      for (iint n=0;n<N;n++)
        V[j+1][n] -= ghij*V[i][n];

      H[i + j*k] = (double) ghij;
    }

    if(j+1 < k){

      dfloat norm_vj = 0.;
      for (iint n=0;n<N;n++)
        norm_vj += V[j+1][n]*V[j+1][n];

      dfloat gNorm_vj;
      MPI_Allreduce(&norm_vj, &gNorm_vj, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
      gNorm_vj = sqrt(gNorm_vj);

      H[j+1+ j*k] = (double) gNorm_vj;

      for (iint n=0;n<N;n++)
        V[j+1][n] *= 1./H[j+1 + j*k];
    }
  }

  double *WR = (double *) calloc(k,sizeof(double));
  double *WI = (double *) calloc(k,sizeof(double));

  eig(k, H, WR, WI);

  double rho = 0.;

  for(int i=0; i<k; i++){
    double rho_i  = sqrt(WR[i]*WR[i] + WI[i]*WI[i]);

    if(rho < rho_i) {
      rho = rho_i;
    }
  }

  free(H);
  free(WR);
  free(WI);

  // free memory
  for(int i=0; i<=k; i++){
    free(V[i]);
  }

  o_Vx.free();
  o_AVx.free();

  return rho;
}