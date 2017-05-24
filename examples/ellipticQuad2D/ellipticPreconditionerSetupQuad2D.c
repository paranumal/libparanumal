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

typedef struct{

  iint row;
  iint col;
  iint ownerRank;
  dfloat val;

}nonZero_t;

// compare on global indices 
int parallelCompareRowColumn(const void *a, const void *b);

void ellipticBuildIpdgQuad2D(mesh2D *mesh, dfloat lambda, nonZero_t **A, iint *nnzA, const char *options);

void ellipticBuildContinuousQuad2D(mesh2D *mesh, dfloat lambda, nonZero_t **A, iint *nnz, hgs_t **hgs, iint *globalStarts, const char* options);

precon_t *ellipticPreconditionerSetupQuad2D(mesh2D *mesh, ogs_t *ogs, dfloat lambda, const char *options){

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  iint Nlocal = mesh->Np*mesh->Nelements;
  iint Nhalo  = mesh->Np*mesh->totalHaloPairs;
  iint Ntrace = mesh->Nfp*mesh->Nfaces*mesh->Nelements;

  precon_t *precon = (precon_t*) calloc(1, sizeof(precon_t));

  if (strstr(options,"OAS")) {
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

    free(diagA);

    // coarse grid preconditioner (only continous elements)
    occaTimerTic(mesh->device,"CoarsePreconditionerSetup");
    ellipticCoarsePreconditionerSetupQuad2D(mesh, precon, lambda, options);
    occaTimerToc(mesh->device,"CoarsePreconditionerSetup");

  } else if (strstr(options,"FULLALMOND")) {
    iint nnz;
    nonZero_t *A;
    hgs_t *hgs;

    iint Nnum = mesh->Np*mesh->Nelements;
    iint *globalStarts = (iint*) calloc(size+1, sizeof(iint));

    if (strstr(options,"IPDG")) {

      MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT, globalStarts+1, 1, MPI_IINT, MPI_COMM_WORLD);
      for(iint r=0;r<size;++r)
        globalStarts[r+1] = globalStarts[r]+globalStarts[r+1]*mesh->Np;

      ellipticBuildIpdgQuad2D(mesh, lambda, &A, &nnz,options);
    
      qsort(A, nnz, sizeof(nonZero_t), parallelCompareRowColumn);

    } else if (strstr(options,"CONTINUOUS")) {
      
      ellipticBuildContinuousQuad2D(mesh,lambda,&A,&nnz,&hgs,globalStarts, options);
    }

    //collect global assembled matrix
    iint *globalnnz       = (iint *) calloc(size  ,sizeof(iint));
    iint *globalnnzOffset = (iint *) calloc(size+1,sizeof(iint));
    MPI_Allgather(&nnz, 1, MPI_IINT, 
                  globalnnz, 1, MPI_IINT, MPI_COMM_WORLD);
    globalnnzOffset[0] = 0;
    for (iint n=0;n<size;n++)
      globalnnzOffset[n+1] = globalnnzOffset[n]+globalnnz[n];

    iint globalnnzTotal = globalnnzOffset[size];

    iint *globalRecvCounts  = (iint *) calloc(size,sizeof(iint));
    iint *globalRecvOffsets = (iint *) calloc(size,sizeof(iint));
    for (iint n=0;n<size;n++){
      globalRecvCounts[n] = globalnnz[n]*sizeof(nonZero_t);
      globalRecvOffsets[n] = globalnnzOffset[n]*sizeof(nonZero_t);
    }
    nonZero_t *globalNonZero = (nonZero_t*) calloc(globalnnzTotal, sizeof(nonZero_t));

    MPI_Allgatherv(A, nnz*sizeof(nonZero_t), MPI_CHAR, 
                  globalNonZero, globalRecvCounts, globalRecvOffsets, MPI_CHAR, MPI_COMM_WORLD);
    

    iint *globalIndex = (iint *) calloc(globalnnzTotal, sizeof(iint));
    iint *globalRows = (iint *) calloc(globalnnzTotal, sizeof(iint));
    iint *globalCols = (iint *) calloc(globalnnzTotal, sizeof(iint));
    dfloat *globalVals = (dfloat*) calloc(globalnnzTotal,sizeof(dfloat));

    for (iint n=0;n<globalnnzTotal;n++) {
      globalRows[n] = globalNonZero[n].row;
      globalCols[n] = globalNonZero[n].col;
      globalVals[n] = globalNonZero[n].val;
    }

    precon->parAlmond = parAlmondSetup(mesh, 
                                   Nnum, 
                                   globalStarts, 
                                   globalnnzTotal,      
                                   globalRows,        
                                   globalCols,        
                                   globalVals,    
                                   0,             // 0 if no null space
                                   hgs,
                                   options);

    precon->o_r1 = mesh->device.malloc(Nnum*sizeof(dfloat));
    precon->o_z1 = mesh->device.malloc(Nnum*sizeof(dfloat));
    precon->r1 = (dfloat*) malloc(Nnum*sizeof(dfloat));
    precon->z1 = (dfloat*) malloc(Nnum*sizeof(dfloat));
  } 

  return precon;
}
