#include "ellipticTet3D.h"

typedef struct{

  dlong localId;
  hlong baseId;
  int haloFlag;

} preconGatherInfo_t;

int parallelCompareBaseId(const void *a, const void *b){

  preconGatherInfo_t *fa = (preconGatherInfo_t*) a;
  preconGatherInfo_t *fb = (preconGatherInfo_t*) b;

  if(fa->baseId < fb->baseId) return -1;
  if(fa->baseId > fb->baseId) return +1;

  return 0;
}

void ellipticSetupSmootherOverlappingPatch(solver_t *solver, precon_t *precon, agmgLevel *level,
                                              dfloat tau, dfloat lambda, int* BCType, const char *options) {

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh3D *mesh = solver->mesh;

  dlong Nlocal = mesh->Np*mesh->Nelements;
  dlong Nhalo  = mesh->Np*mesh->totalHaloPairs;
  // dlong Ntrace = mesh->Nfp*mesh->Nfaces*mesh->Nelements;

  // build gather-scatter
  int NpP = mesh->Np + mesh->Nfaces*mesh->Nfp;

  // build gather-scatter for overlapping patches
  dlong *allNelements = (dlong*) calloc(size, sizeof(dlong));
  MPI_Allgather(&(mesh->Nelements), 1, MPI_DLONG,
                      allNelements, 1, MPI_DLONG, MPI_COMM_WORLD);

  // offsets
  hlong *startElement = (hlong*) calloc(size, sizeof(hlong));
  for(int r=1;r<size;++r){
    startElement[r] = startElement[r-1]+allNelements[r-1];
  }

  // 1-indexed numbering of nodes on this process
  hlong *localNums = (hlong*) calloc((Nlocal+Nhalo), sizeof(hlong));
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      localNums[e*mesh->Np+n] = 1 + e*mesh->Np + n + startElement[rank]*mesh->Np;
    }
  }

  if(Nhalo){
    // send buffer for outgoing halo
    hlong *sendBuffer = (hlong*) calloc(Nhalo, sizeof(hlong));

    // exchange node numbers with neighbors
    meshHaloExchange(mesh,
         mesh->Np*sizeof(hlong),
         localNums,
         sendBuffer,
         localNums+Nlocal);
  }

  preconGatherInfo_t *preconGatherInfoDg =
    (preconGatherInfo_t*) calloc(NpP*mesh->Nelements,
         sizeof(preconGatherInfo_t));

  // set local ids
  for(dlong n=0;n<mesh->Nelements*NpP;++n)
    preconGatherInfoDg[n].localId = n;

  // numbering of patch interior nodes
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dlong id  = n + e*mesh->Np;
      dlong pid = n + e*NpP;

      // all patch interior nodes are local
      preconGatherInfoDg[pid].baseId = localNums[id];
    }
  }

  // add patch boundary nodes
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      // mark halo nodes
      int rP = mesh->EToP[e*mesh->Nfaces+f];
      // dlong eP = mesh->EToE[e*mesh->Nfaces+f];
      // int fP = mesh->EToF[e*mesh->Nfaces+f];
      // int bc = mesh->EToB[e*mesh->Nfaces+f];

      for(int n=0;n<mesh->Nfp;++n){
        dlong id = n + f*mesh->Nfp+e*mesh->Nfp*mesh->Nfaces;
        dlong idP = mesh->vmapP[id];

        // local numbers
        dlong pidM = e*NpP + mesh->faceNodes[f*mesh->Nfp+n];
        dlong pidP = e*NpP + mesh->Np + f*mesh->Nfp+n;
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
  dlong skip = 0;

  while(preconGatherInfoDg[skip].baseId==0 && skip<NpP*mesh->Nelements){
    ++skip;
  }

  // reset local ids
  dlong NlocalDg = NpP*mesh->Nelements - skip;
  dlong *gatherLocalIdsDg  = (dlong*) calloc(NlocalDg, sizeof(dlong));
  hlong *gatherBaseIdsDg   = (hlong*) calloc(NlocalDg, sizeof(hlong));
  int *gatherHaloFlagsDg = (int*) calloc(NlocalDg, sizeof(int));
  for(dlong n=0;n<NlocalDg;++n){
    gatherLocalIdsDg[n]  = preconGatherInfoDg[n+skip].localId;
    gatherBaseIdsDg[n]   = preconGatherInfoDg[n+skip].baseId;
    gatherHaloFlagsDg[n] = preconGatherInfoDg[n+skip].haloFlag;
  }

  // make preconBaseIds => preconNumbering

  /* depreciated */
  // precon->ogsDg = meshParallelGatherScatterSetup(mesh,
  //            NlocalDg,
  //            sizeof(dfloat),
  //            gatherLocalIdsDg,
  //            gatherBaseIdsDg,
  //            gatherHaloFlagsDg);

  //correction for full patch
  NpP = mesh->NpP;

  // build degree vector
  dlong NtotalDGP = NpP*mesh->Nelements;
  dfloat *invDegree = (dfloat*) calloc(NtotalDGP, sizeof(dfloat));
  dfloat *degree    = (dfloat*) calloc(NtotalDGP, sizeof(dfloat));
  precon->o_invDegreeDGP = mesh->device.malloc(NtotalDGP*sizeof(dfloat), invDegree);

  for(dlong n=0;n<NtotalDGP;++n)
    degree[n] = 1;

  occa::memory o_deg = mesh->device.malloc(NtotalDGP*sizeof(dfloat), degree);
  meshParallelGatherScatter(mesh, precon->ogsDg, o_deg);
  o_deg.copyTo(degree);
  mesh->device.finish();
  o_deg.free();

  for(dlong n=0;n<NtotalDGP;++n){ // need to weight inner products{
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
  for(int n=0;n<NpP;++n){
    for(int m=0;m<NpP;++m){
      oasForwardDgT[n+m*NpP] = mesh->oasForwardDg[m+n*NpP];
      oasBackDgT[n+m*NpP] = mesh->oasBackDg[m+n*NpP];
    }
  }

  precon->o_oasForwardDgT = mesh->device.malloc(NpP*NpP*sizeof(dfloat), oasForwardDgT);
  precon->o_oasBackDgT    = mesh->device.malloc(NpP*NpP*sizeof(dfloat), oasBackDgT);


  /// ---------------------------------------------------------------------------

  // hack estimate for Jacobian scaling

  // dfloat *diagInvOp = (dfloat*) calloc(NpP*mesh->Nelements, sizeof(dfloat));
  dfloat *diagInvOpDg = (dfloat*) calloc(NpP*mesh->Nelements, sizeof(dfloat));
  for(dlong e=0;e<mesh->Nelements;++e){

    dfloat J = mesh->vgeo[e*mesh->Nvgeo + JID];
    dfloat rx = mesh->vgeo[e*mesh->Nvgeo + RXID];
    dfloat sx = mesh->vgeo[e*mesh->Nvgeo + SXID];
    dfloat tx = mesh->vgeo[e*mesh->Nvgeo + TXID];
    dfloat ry = mesh->vgeo[e*mesh->Nvgeo + RYID];
    dfloat sy = mesh->vgeo[e*mesh->Nvgeo + SYID];
    dfloat ty = mesh->vgeo[e*mesh->Nvgeo + TYID];
    dfloat rz = mesh->vgeo[e*mesh->Nvgeo + RZID];
    dfloat sz = mesh->vgeo[e*mesh->Nvgeo + SZID];
    dfloat tz = mesh->vgeo[e*mesh->Nvgeo + TZID];

    //TODO change this to the eigenvalues of the element metric tensor
    dfloat Jhrinv2 = J*(rx*rx+ry*ry+rz*rz);
    dfloat Jhsinv2 = J*(sx*sx+sy*sy+sz*sz);
    dfloat Jhtinv2 = J*(tx*tx+ty*ty+tz*tz);
    dfloat Jhinv2 = (Jhrinv2 + Jhsinv2 + Jhtinv2)/3;
    //dfloat Jhinv2 = mymax(mymax(Jhrinv2, Jhsinv2),Jhtinv2);
    for(int n=0;n<NpP;++n){
      dlong pid = n + e*NpP;

      diagInvOpDg[pid] = 1./(J*lambda + Jhinv2*mesh->oasDiagOpDg[n]);
    }
  }

  precon->o_oasDiagInvOpDg =
    mesh->device.malloc(NpP*mesh->Nelements*sizeof(dfloat), diagInvOpDg);

  //storage buffer for patches
  dlong NtotalP = mesh->NpP*mesh->Nelements;
  precon->zP  = (dfloat*) calloc(NtotalP,  sizeof(dfloat));
  precon->o_zP  = mesh->device.malloc(NtotalP*sizeof(dfloat),precon->zP);


  level->device_smoother = overlappingPatchIpdg;

  //check if stabilization is needed
  if (strstr(options,"MULTIGRID")||strstr(options,"FULLALMOND")) {
    //estimate the max eigenvalue of S*A
    dfloat rho = maxEigSmoothAx(solver, level, options);

    if (strstr(options,"CHEBYSHEV")) {

      level->smoother_params = (dfloat *) calloc(2,sizeof(dfloat));

      level->ChebyshevIterations = 2;
      level->smoother_params[0] = rho;
      level->smoother_params[1] = rho/10.;

    } else {

      //set the stabilty weight (jacobi-type interation)
      dfloat weight = (4./3.)/rho;

      // printf("weight = %g \n", weight);

      for (int n=0;n<NpP*mesh->Nelements;n++)
        diagInvOpDg[n] *= weight;

      //update diagonal with weight
      precon->o_oasDiagInvOpDg.copyFrom(diagInvOpDg);
    }
  }
}

void ellipticSetupSmootherFullPatch(solver_t *solver, precon_t *precon, agmgLevel *level,
                                              dfloat tau, dfloat lambda, int* BCType, dfloat rateTolerance, const char *options) {

  dfloat *invAP;
  dlong Npatches;
  dlong *patchesIndex;
  mesh_t *mesh = solver->mesh;

  int NpP = mesh->Np*(mesh->Nfaces+1);

  //initialize the full inverse operators on each 4 element patch
  ellipticBuildFullPatchesTet3D(solver, mesh, mesh->Np, NULL, tau, lambda, BCType, rateTolerance,
                                      &Npatches, &patchesIndex, &invAP, options);

  precon->o_invAP = mesh->device.malloc(Npatches*NpP*NpP*sizeof(dfloat),invAP);
  precon->o_patchesIndex = mesh->device.malloc(mesh->Nelements*sizeof(int), patchesIndex);

  dfloat *invDegree = (dfloat*) calloc(mesh->Nelements,sizeof(dfloat));
  for (dlong e=0;e<mesh->Nelements;e++) {
    for (int f=0;f<mesh->Nfaces;f++)
        invDegree[e] += (mesh->EToE[e*mesh->Nfaces +f]<0) ? 0 : 1; //overlap degree = # of neighbours
    invDegree[e] = 1.0/invDegree[e]; //build in weight
  }
  precon->o_invDegreeAP = mesh->device.malloc(mesh->Nelements*sizeof(dfloat),invDegree);

  mesh->o_EToE = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(int),mesh->EToE);
  mesh->o_EToF = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(int),mesh->EToF);

  //set storage for larger patch
  precon->zP = (dfloat*) calloc(mesh->Nelements*NpP,  sizeof(dfloat));
  precon->o_zP = mesh->device.malloc(mesh->Nelements*NpP*sizeof(dfloat), precon->zP);

  level->device_smoother = FullPatchIpdg;

  //check if stabilization is needed
  if (strstr(options,"MULTIGRID")||strstr(options,"FULLALMOND")) {
    //estimate the max eigenvalue of S*A
    dfloat rho = maxEigSmoothAx(solver, level, options);

    if (strstr(options,"CHEBYSHEV")) {

      level->smoother_params = (dfloat *) calloc(2,sizeof(dfloat));

      level->ChebyshevIterations = 2;
      level->smoother_params[0] = rho;
      level->smoother_params[1] = rho/10.;

    } else {

      //set the stabilty weight (jacobi-type interation)
      dfloat weight = (4./3.)/rho;

      for (int e=0;e<mesh->Nelements;e++)
        invDegree[e] *= weight;

      //update with weight
      precon->o_invDegreeAP.copyFrom(invDegree);
    }
  }
  free(invDegree);
}

void ellipticSetupSmootherFacePatch(solver_t *solver, precon_t *precon, agmgLevel *level,
                                              dfloat tau, dfloat lambda, int* BCType, dfloat rateTolerance, const char *options) {

  dfloat *invAP;
  dlong Npatches;
  dlong *patchesIndex;
  mesh_t *mesh = solver->mesh;

  //initialize the full inverse operators on each 4 element patch
  ellipticBuildFacePatchesTet3D(solver, mesh, mesh->Np, NULL, tau, lambda, BCType, rateTolerance,
                                      &Npatches, &patchesIndex, &invAP, options);

  int NpP = 2*mesh->Np;

  precon->o_invAP = mesh->device.malloc(Npatches*NpP*NpP*sizeof(dfloat),invAP);
  precon->o_patchesIndex = mesh->device.malloc(mesh->NfacePairs*sizeof(dlong), patchesIndex);

  dfloat *invDegree = (dfloat*) calloc(mesh->Nelements+mesh->totalHaloPairs,sizeof(dfloat));
  for (dlong face=0;face<mesh->NfacePairs;face++) {
    dlong eM = mesh->FPairsToE[2*face+0];
    dlong eP = mesh->FPairsToE[2*face+1];

    invDegree[eM]++; //overlap degree = # of patches
    if (eP>=0) invDegree[eP]++; //overlap degree = # of patches
  }
  for (dlong e=0;e<mesh->Nelements+mesh->totalHaloPairs;e++) {
    invDegree[e] = 1.0/invDegree[e];
  }

  precon->o_invDegreeAP = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*sizeof(dfloat),invDegree);

  mesh->o_FPairsToE = mesh->device.malloc(2*mesh->NfacePairs*sizeof(dlong),mesh->FPairsToE);
  mesh->o_FPairsToF = mesh->device.malloc(2*mesh->NfacePairs*sizeof(int),mesh->FPairsToF);
  mesh->o_EToFPairs = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(dlong),mesh->EToFPairs);

  //set storage for larger patch
  precon->zP = (dfloat*) calloc(mesh->NfacePairs*NpP,  sizeof(dfloat));
  precon->o_zP = mesh->device.malloc(mesh->NfacePairs*NpP*sizeof(dfloat), precon->zP);


  level->device_smoother = FacePatchIpdg;

  //check if stabilization is needed
  if (strstr(options,"MULTIGRID")||strstr(options,"FULLALMOND")) {
    //estimate the max eigenvalue of S*A
    dfloat rho = maxEigSmoothAx(solver, level, options);

    if (strstr(options,"CHEBYSHEV")) {

      level->smoother_params = (dfloat *) calloc(2,sizeof(dfloat));

      level->ChebyshevIterations = 2;
      level->smoother_params[0] = rho;
      level->smoother_params[1] = rho/10.;

    } else {

      //set the stabilty weight (jacobi-type interation)
      dfloat weight = (4./3.)/rho;

      for (int e=0;e<mesh->Nelements;e++)
        invDegree[e] *= weight;

      //update with weight
      precon->o_invDegreeAP.copyFrom(invDegree);
    }
  }
  free(invDegree);
}

void ellipticSetupSmootherLocalPatch(solver_t *solver, precon_t *precon, agmgLevel *level,
                                              dfloat tau, dfloat lambda, int* BCType, dfloat rateTolerance, const char *options) {

  dfloat *invAP;
  dlong Npatches;
  dlong *patchesIndex;
  mesh_t *mesh = solver->mesh;

  int NpP = mesh->Np;

  //initialize the full inverse operators on each 4 element patch
  ellipticBuildLocalPatchesTet3D(solver, mesh, mesh->Np, NULL, tau, lambda, BCType, rateTolerance,
                                      &Npatches, &patchesIndex, &invAP, options);

  precon->o_invAP = mesh->device.malloc(Npatches*NpP*NpP*sizeof(dfloat),invAP);
  precon->o_patchesIndex = mesh->device.malloc(mesh->Nelements*sizeof(dlong), patchesIndex);

  dfloat *invDegree = (dfloat*) calloc(mesh->Nelements,sizeof(dfloat));
  for (dlong e=0;e<mesh->Nelements;e++) {
    invDegree[e] = 1.0;
  }
  precon->o_invDegreeAP = mesh->device.malloc(mesh->Nelements*sizeof(dfloat),invDegree);

  level->device_smoother = LocalPatchIpdg;

  //check if stabilization is needed
  if (strstr(options,"MULTIGRID")||strstr(options,"FULLALMOND")) {
    //estimate the max eigenvalue of S*A
    dfloat rho = maxEigSmoothAx(solver, level, options);

    if (strstr(options,"CHEBYSHEV")) {

      level->smoother_params = (dfloat *) calloc(2,sizeof(dfloat));

      level->ChebyshevIterations = 2;
      level->smoother_params[0] = rho;
      level->smoother_params[1] = rho/10.;

    } else {

      //set the stabilty weight (jacobi-type interation)
      dfloat weight = (4./3.)/rho;

      for (int e=0;e<mesh->Nelements;e++)
        invDegree[e] *= weight;

      //update with weight
      precon->o_invDegreeAP.copyFrom(invDegree);
    }
  }
  free(invDegree);
}

void ellipticSetupSmootherDampedJacobi(solver_t *solver, precon_t *precon, agmgLevel *level,
                                              dfloat tau, dfloat lambda, int* BCType, const char *options) {

  dfloat *invDiagA;
  mesh_t *mesh = solver->mesh;

  ellipticBuildJacobiTet3D(solver, mesh,mesh->Np,NULL,tau, lambda, BCType, &invDiagA,options);

  precon->o_invDiagA = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), invDiagA);

  level->device_smoother = dampedJacobi;

  //check if stabilization is needed
  if (strstr(options,"MULTIGRID")||strstr(options,"FULLALMOND")) {
    //estimate the max eigenvalue of S*A
    dfloat rho = maxEigSmoothAx(solver, level, options);

    if (strstr(options,"CHEBYSHEV")) {

      level->smoother_params = (dfloat *) calloc(2,sizeof(dfloat));

      level->ChebyshevIterations = 2;
      level->smoother_params[0] = rho;
      level->smoother_params[1] = rho/10.;

    } else {

      //set the stabilty weight (jacobi-type interation)
      dfloat weight = (4./3.)/rho;

      for (dlong n=0;n<mesh->Np*mesh->Nelements;n++)
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

dfloat maxEigSmoothAx(solver_t* solver, agmgLevel *level, const char* options){

  mesh_t *mesh = solver->mesh;

  const dlong N = level->Nrows;
  const dlong M = level->Ncols;

  int k = 10;

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  hlong Nlocal = (hlong) level->Nrows;
  hlong Ntotal = 0;
  MPI_Allreduce(&Nlocal, &Ntotal, 1, MPI_HLONG, MPI_SUM, MPI_COMM_WORLD);
  if(k > Ntotal) k = (int) Ntotal;

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  double *H = (double *) calloc(k*k,sizeof(double));

  // allocate memory for basis
  dfloat *Vx = (dfloat*) calloc(M, sizeof(dfloat));
  occa::memory *o_V = (occa::memory *) calloc(k+1, sizeof(occa::memory));
  
  occa::memory o_Vx  = mesh->device.malloc(M*sizeof(dfloat),Vx);
  occa::memory o_AVx = mesh->device.malloc(M*sizeof(dfloat),Vx);

  for(int i=0; i<=k; i++)
    o_V[i] = mesh->device.malloc(M*sizeof(dfloat),Vx);

  // generate a random vector for initial basis vector
  for (dlong i=0;i<N;i++) Vx[i] = (dfloat) drand48(); 

  //gather-scatter 
  if (strstr(options,"CONTINUOUS")) {
    gsParallelGatherScatter(mesh->hostGsh, Vx, dfloatString, "add");   
    for (dlong i=0;i<solver->Nmasked;i++) Vx[solver->maskIds[i]] = 0.;
  }

  o_Vx.copyFrom(Vx); //copy to device
  dfloat norm_vo = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_Vx, o_Vx, options);
  norm_vo = sqrt(norm_vo);

  ellipticScaledAdd(solver, 1./norm_vo, o_Vx, 0. , o_V[0]);

  for(int j=0; j<k; j++){
    // v[j+1] = invD*(A*v[j])
    level->device_Ax(level->AxArgs,o_V[j],o_AVx);
    level->device_smoother(level->smootherArgs, o_AVx, o_V[j+1]);

    // modified Gram-Schmidth
    for(int i=0; i<=j; i++){
      // H(i,j) = v[i]'*A*v[j]
      dfloat hij = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_V[i], o_V[j+1], options);

      // v[j+1] = v[j+1] - hij*v[i]
      ellipticScaledAdd(solver, -hij, o_V[i], 1., o_V[j+1]);

      H[i + j*k] = (double) hij;
    }

    if(j+1 < k){
      // v[j+1] = v[j+1]/||v[j+1]||
      dfloat norm_vj = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_V[j+1], o_V[j+1], options);
      norm_vj = sqrt(norm_vj);
      ellipticScaledAdd(solver, 1/norm_vj, o_V[j+1], 0., o_V[j+1]);
      
      H[j+1+ j*k] = (double) norm_vj;
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

  // free memory
  free(H);
  free(WR);
  free(WI);

  free(Vx);
  o_Vx.free();
  o_AVx.free();
  for(int i=0; i<=k; i++) o_V[i].free();
  free((void*)o_V);

  if((rank==0)&&strstr(options,"VERBOSE")) printf("weight = %g \n", rho);

  return rho;
}