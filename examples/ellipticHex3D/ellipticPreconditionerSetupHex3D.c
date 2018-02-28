#include "ellipticHex3D.h"

typedef struct{

  int localId;
  int baseId;
  int haloFlag;
  
} preconGatherInfo_t;

int parallelCompareBaseId(const void *a, const void *b){

  preconGatherInfo_t *fa = (preconGatherInfo_t*) a;
  preconGatherInfo_t *fb = (preconGatherInfo_t*) b;

  if(fa->baseId < fb->baseId) return -1;
  if(fa->baseId > fb->baseId) return +1;

  return 0;

}

typedef struct{

  int row;
  int col;
  int ownerRank;
  dfloat val;

}nonZero_t;

// compare on global indices 
int parallelCompareRowColumn(const void *a, const void *b);

void ellipticBuildIpdgHex3D(mesh3D *mesh, dfloat lambda, nonZero_t **A, int *nnzA, const char *options);

void ellipticBuildContinuousHex3D(mesh3D *mesh, dfloat lambda, nonZero_t **A, int *nnz, hgs_t **hgs, int *globalStarts, const char* options);

precon_t *ellipticPreconditionerSetupHex3D(mesh3D *mesh, ogs_t *ogs, dfloat lambda, const char *options){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // assumes meshParallelGatherScatterSetup3D has been called
  
  // ????? need to extend storage for halo ?????
  int Nlocal = mesh->Np*mesh->Nelements;
  int Nhalo  = mesh->Np*mesh->totalHaloPairs;
  int Ntrace = mesh->Nfp*mesh->Nfaces*mesh->Nelements;

  precon_t *precon = (precon_t*) calloc(1, sizeof(precon_t));

  if (strstr(options,"OAS")) {

    // offsets to extract second layer
    int *offset = (int*) calloc(mesh->Nfaces, sizeof(int));
    offset[0] = +mesh->Nfp;
    offset[1] = +mesh->Nq;
    offset[2] = -1;
    offset[3] = -mesh->Nq;
    offset[4] = +1;
    offset[5] = -mesh->Nfp;

    // build gather-scatter
    int NqP = mesh->Nq+2;
    int NpP = NqP*NqP*NqP;

    int *offsetP = (int*) calloc(mesh->Nfaces, sizeof(int));
    offsetP[0] = +NqP*NqP;
    offsetP[1] = +NqP;
    offsetP[2] = -1;
    offsetP[3] = -NqP;
    offsetP[4] = +1;
    offsetP[5] = -NqP*NqP;

    int *faceNodesPrecon = (int*) calloc(mesh->Nfp*mesh->Nfaces, sizeof(int));
    for(int j=0;j<mesh->Nq;++j)
      for(int i=0;i<mesh->Nq;++i)
        faceNodesPrecon[i+j*mesh->Nq+0*mesh->Nfp] = i+1 + (j+1)*NqP + 0*NqP*NqP;
    
    for(int k=0;k<mesh->Nq;++k)
      for(int i=0;i<mesh->Nq;++i)
        faceNodesPrecon[i+k*mesh->Nq+1*mesh->Nfp] = i+1 + 0*NqP + (k+1)*NqP*NqP;
    
    for(int k=0;k<mesh->Nq;++k)
      for(int j=0;j<mesh->Nq;++j)
        faceNodesPrecon[j+k*mesh->Nq+2*mesh->Nfp] = NqP-1 + (j+1)*NqP + (k+1)*NqP*NqP;

    for(int k=0;k<mesh->Nq;++k)
      for(int i=0;i<mesh->Nq;++i)
        faceNodesPrecon[i+k*mesh->Nq+3*mesh->Nfp] = i+1 + (NqP-1)*NqP + (k+1)*NqP*NqP;

    for(int k=0;k<mesh->Nq;++k)
      for(int j=0;j<mesh->Nq;++j)
        faceNodesPrecon[j+k*mesh->Nq+4*mesh->Nfp] = 0 + (j+1)*NqP + (k+1)*NqP*NqP;
    
    for(int j=0;j<mesh->Nq;++j)
      for(int i=0;i<mesh->Nq;++i)
        faceNodesPrecon[i+j*mesh->Nq+5*mesh->Nfp] = i+1 + (j+1)*NqP + (NqP-1)*NqP*NqP;

    int *vmapMP = (int*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(int));
    int *vmapPP = (int*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(int));
    
    // take node info from positive trace and put in overlap region on parallel gather info
    for(int e=0;e<mesh->Nelements;++e){

      for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
        int fid = e*mesh->Nfp*mesh->Nfaces + n;

        // find face index
        int fM = n/mesh->Nfp;
        int fP = mesh->EToF[e*mesh->Nfaces+fM];
        if(fP<0) fP = fM;

        // find location of "+" trace in regular mesh
        int idM = mesh->vmapM[fid] + offset[fM];
        int idP = mesh->vmapP[fid] + offset[fP];

        vmapMP[fid] = idM;
        vmapPP[fid] = idP;
      }
    }
    
    // space for gather base indices with halo
    preconGatherInfo_t *gatherInfo =
      (preconGatherInfo_t*) calloc(Nlocal+Nhalo, sizeof(preconGatherInfo_t));

    // rearrange in node order
    for(int n=0;n<Nlocal;++n){
      int id = mesh->gatherLocalIds[n];
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
    for(int e=0;e<mesh->Nelements;++e){
      for(int k=0;k<mesh->Nq;++k){
        for(int j=0;j<mesh->Nq;++j){
  	for(int i=0;i<mesh->Nq;++i){
  	  int id  = i + j*mesh->Nq + k*mesh->Nq*mesh->Nq + e*mesh->Np;
  	  int pid = i + 1 + (j+1)*NqP + (k+1)*NqP*NqP + e*NpP;

  	  preconGatherInfo[pid] = gatherInfo[id];
  	}
        }
      }
    }

    // add overlap region
    for(int e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
        int id = n + e*mesh->Nfp*mesh->Nfaces;

        int f = n/mesh->Nfp;
        int idM = mesh->vmapM[id];
        int idP = vmapPP[id];
        int idMP = e*NpP + faceNodesPrecon[n];
  	    
        preconGatherInfo[idMP] = gatherInfo[idP];

        if(gatherInfo[idM].haloFlag){
  	preconGatherInfo[idMP].haloFlag = 1;
  	preconGatherInfo[idMP+offsetP[f]].haloFlag = 1;
  	preconGatherInfo[idMP+2*offsetP[f]].haloFlag = 1;
        }
      }
    }
    
    // reset local ids
    for(int n=0;n<mesh->Nelements*NpP;++n)
      preconGatherInfo[n].localId = n;

    char fname[BUFSIZ];
    sprintf(fname, "haloFlag%05d.dat", rank);
    FILE *fp = fopen(fname, "w");
    
    for(int p=0;p<size;++p){
      if(p==rank){
        for(int e=0;e<mesh->Nelements;++e){
  	fprintf(fp,"e=%d: \n", e);
  	for(int k=0;k<mesh->NqP;++k){
  	  for(int j=0;j<mesh->NqP;++j){
  	    for(int i=0;i<mesh->NqP;++i){
  	      int id = i + mesh->NqP*j + mesh->NqP*mesh->NqP*k+ e*NpP;
  	      fprintf(fp,"%d ", preconGatherInfo[id].haloFlag);
  	    }
  	    fprintf(fp,"\n");
  	  }
  	  fprintf(fp,"\n");
  	}
  	fprintf(fp,"\n");
        }
      }
    }
    fclose(fp);
    
    // sort by rank then base index
    qsort(preconGatherInfo, NpP*mesh->Nelements, sizeof(preconGatherInfo_t), parallelCompareBaseId);
    
    // do not gather-scatter nodes labelled zero
    int skip = 0;
    while(preconGatherInfo[skip].baseId==0 && skip<NpP*mesh->Nelements){
      ++skip;
    }
    printf("skip = %d out of %d\n", skip, NpP*mesh->Nelements);

    // reset local ids
    int NlocalP = NpP*mesh->Nelements - skip;
    int *gatherLocalIdsP  = (int*) calloc(NlocalP, sizeof(int));
    int *gatherBaseIdsP   = (int*) calloc(NlocalP, sizeof(int));
    int *gatherHaloFlagsP = (int*) calloc(NlocalP, sizeof(int));
    for(int n=0;n<NlocalP;++n){
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
    int NtotalP = mesh->NqP*mesh->NqP*mesh->NqP*mesh->Nelements;
    dfloat *invDegree = (dfloat*) calloc(NtotalP, sizeof(dfloat));
    dfloat *degree    = (dfloat*) calloc(NtotalP, sizeof(dfloat));
    precon->o_invDegreeP = mesh->device.malloc(NtotalP*sizeof(dfloat), invDegree);
    
    for(int n=0;n<NtotalP;++n)
      degree[n] = 1;

    occa::memory o_deg = mesh->device.malloc(NtotalP*sizeof(dfloat), degree);
    meshParallelGatherScatter(mesh, precon->ogsP, o_deg, o_deg, dfloatString, "add");
    o_deg.copyTo(degree);
    mesh->device.finish();
    o_deg.free();

    for(int n=0;n<NtotalP;++n){ // need to weight inner products{
      if(degree[n] == 0) printf("WARNING!!!!\n");
      invDegree[n] = 1./degree[n];
    }
    
    precon->o_invDegreeP.copyFrom(invDegree);
    free(degree);
    free(invDegree);

    // -------------------------------------------------------------------------------------------
    // build gather-scatter for overlapping patches
    int *allNelements = (int*) calloc(size, sizeof(int));
    MPI_Allgather(&(mesh->Nelements), 1, MPI_INT,
  		allNelements, 1, MPI_INT, MPI_COMM_WORLD);

    // offsets
    int *startElement = (int*) calloc(size, sizeof(int));
    for(int r=1;r<size;++r){
      startElement[r] = startElement[r-1]+allNelements[r-1];
    }

    // 1-indexed numbering of nodes on this process
    int *localNums = (int*) calloc((Nlocal+Nhalo), sizeof(int));
    for(int e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        localNums[e*mesh->Np+n] = 1 + e*mesh->Np + n + startElement[rank]*mesh->Np;
      }
    }
    
    if(Nhalo){
      // send buffer for outgoing halo
      int *sendBuffer = (int*) calloc(Nhalo, sizeof(int));

      // exchange node numbers with neighbors
      meshHaloExchange(mesh,
  		     mesh->Np*sizeof(int),
  		     localNums,
  		     sendBuffer,
  		     localNums+Nlocal);
    }
    
    preconGatherInfo_t *preconGatherInfoDg = 
      (preconGatherInfo_t*) calloc(NpP*mesh->Nelements,
  				 sizeof(preconGatherInfo_t));

    // set local ids
    for(int n=0;n<mesh->Nelements*NpP;++n)
      preconGatherInfoDg[n].localId = n;

    // numbering of patch interior nodes
    for(int e=0;e<mesh->Nelements;++e){
      for(int k=0;k<mesh->Nq;++k){
        for(int j=0;j<mesh->Nq;++j){
  	for(int i=0;i<mesh->Nq;++i){
  	  int id  = i + j*mesh->Nq + k*mesh->Nq*mesh->Nq + e*mesh->Np;
  	  int pid = (i+1) + (j+1)*NqP + (k+1)*NqP*NqP + e*NpP;

  	  // all patch interior nodes are local
  	  preconGatherInfoDg[pid].baseId = localNums[id];
  	}
        }
      }
    }
    // add patch boundary nodes
    for(int e=0;e<mesh->Nelements;++e){
      for(int f=0;f<mesh->Nfaces;++f){
        // mark halo nodes
        int rP = mesh->EToP[e*mesh->Nfaces+f];
        int eP = mesh->EToE[e*mesh->Nfaces+f];
        int fP = mesh->EToF[e*mesh->Nfaces+f];
        int bc = mesh->EToB[e*mesh->Nfaces+f];
        
        for(int n=0;n<mesh->Nfp;++n){
  	int id = n + f*mesh->Nfp+e*mesh->Nfp*mesh->Nfaces;
  	int idP = mesh->vmapP[id];
  	
  	// local numbers
  	int pidM = e*NpP + faceNodesPrecon[f*mesh->Nfp+n] + offsetP[f]; 
  	int pidP = e*NpP + faceNodesPrecon[f*mesh->Nfp+n];
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
    int NlocalDg = NpP*mesh->Nelements - skip;
    int *gatherLocalIdsDg  = (int*) calloc(NlocalDg, sizeof(int));
    int *gatherBaseIdsDg   = (int*) calloc(NlocalDg, sizeof(int));
    int *gatherHaloFlagsDg = (int*) calloc(NlocalDg, sizeof(int));
    for(int n=0;n<NlocalDg;++n){
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
    int NtotalDGP = NpP*mesh->Nelements;
    invDegree = (dfloat*) calloc(NtotalDGP, sizeof(dfloat));
    degree    = (dfloat*) calloc(NtotalDGP, sizeof(dfloat));
    precon->o_invDegreeDGP = mesh->device.malloc(NtotalDGP*sizeof(dfloat), invDegree);
    
    for(int n=0;n<NtotalDGP;++n)
      degree[n] = 1;

    o_deg = mesh->device.malloc(NtotalDGP*sizeof(dfloat), degree);
    meshParallelGatherScatter(mesh, precon->ogsDg, o_deg, o_deg, dfloatString, "add");
    o_deg.copyTo(degree);
    mesh->device.finish();
    o_deg.free();

    for(int n=0;n<NtotalDGP;++n){ // need to weight inner products{
      if(degree[n] == 0) printf("WARNING!!!!\n");
      invDegree[n] = 1./degree[n];
    }
    
    precon->o_invDegreeDGP.copyFrom(invDegree);
    free(degree);
    free(invDegree);

    // -------------------------------------------------------------------------------------------
    

    
    precon->o_faceNodesP = mesh->device.malloc(mesh->Nfp*mesh->Nfaces*sizeof(int), faceNodesPrecon);
    precon->o_vmapPP     = mesh->device.malloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements*sizeof(int), vmapPP);

    precon->o_oasForward = mesh->device.malloc(NqP*NqP*sizeof(dfloat), mesh->oasForward);
    precon->o_oasBack    = mesh->device.malloc(NqP*NqP*sizeof(dfloat), mesh->oasBack);

    precon->o_oasForwardDg = mesh->device.malloc(NqP*NqP*sizeof(dfloat), mesh->oasForwardDg);
    precon->o_oasBackDg    = mesh->device.malloc(NqP*NqP*sizeof(dfloat), mesh->oasBackDg);
    
    /// ---------------------------------------------------------------------------
    
    // hack estimate for Jacobian scaling
    dfloat *diagInvOp = (dfloat*) calloc(NpP*mesh->Nelements, sizeof(dfloat));
    dfloat *diagInvOpDg = (dfloat*) calloc(NpP*mesh->Nelements, sizeof(dfloat));
  					  
    for(int e=0;e<mesh->Nelements;++e){

      // S = Jabc*(wa*wb*wc*lambda + wb*wc*Da'*wa*Da + wa*wc*Db'*wb*Db + wa*wb*Dc'*wc*Dc)
      // S = Jabc*wa*wb*wc*(lambda*I+1/wa*Da'*wa*Da + 1/wb*Db'*wb*Db + 1/wc*Dc'*wc*Dc)
      
      dfloat Jhrinv2 = 0, Jhsinv2 = 0, Jhtinv2 = 0, J = 0;
      for(int n=0;n<mesh->Np;++n){
        dfloat W = mesh->gllw[n%mesh->Nq]*
  	mesh->gllw[(n/mesh->Nq)%mesh->Nq]*
  	mesh->gllw[n/(mesh->Nq*mesh->Nq)];
        int base = mesh->Nggeo*mesh->Np*e + n;

        J = mymax(J, mesh->ggeo[base + mesh->Np*GWJID]/W);
        Jhrinv2 = mymax(Jhrinv2, mesh->ggeo[base + mesh->Np*G00ID]/W);
        Jhsinv2 = mymax(Jhsinv2, mesh->ggeo[base + mesh->Np*G11ID]/W);
        Jhtinv2 = mymax(Jhtinv2, mesh->ggeo[base + mesh->Np*G22ID]/W);
        
      }
      
      for(int k=0;k<NqP;++k){
        for(int j=0;j<NqP;++j){
  	for(int i=0;i<NqP;++i){
  	  int pid = i + j*NqP + k*NqP*NqP + e*NpP;
  	  
  	  diagInvOp[pid] =
  	    1./(J*lambda +
  		Jhrinv2*mesh->oasDiagOp[i] +
  		Jhsinv2*mesh->oasDiagOp[j] +
  		Jhtinv2*mesh->oasDiagOp[k]);


  	  diagInvOpDg[pid] =
  	    1./(J*lambda +
  		Jhrinv2*mesh->oasDiagOpDg[i] +
  		Jhsinv2*mesh->oasDiagOpDg[j] +
  		Jhtinv2*mesh->oasDiagOpDg[k]);

  	}
        }
      }
    }
    
    precon->o_oasDiagInvOp = mesh->device.malloc(NpP*mesh->Nelements*sizeof(dfloat), diagInvOp);
    precon->o_oasDiagInvOpDg = mesh->device.malloc(NpP*mesh->Nelements*sizeof(dfloat), diagInvOpDg);
    
    /// ---------------------------------------------------------------------------
    // compute diagonal of stiffness matrix for Jacobi

    int Ntotal = mesh->Np*mesh->Nelements;
    dfloat *diagA = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  				   
    for(int e=0;e<mesh->Nelements;++e){
      int cnt = 0;
      for(int k=0;k<mesh->Nq;++k){
        for(int j=0;j<mesh->Nq;++j){
  	for(int i=0;i<mesh->Nq;++i){
  	  
  	  dfloat JW = mesh->ggeo[e*mesh->Np*mesh->Nggeo+ cnt + mesh->Np*GWJID];
  	  // (D_{ii}^2 + D_{jj}^2 + D_{kk}^2 + lambda)*w_i*w_j*w_k*J_{ijke}
  	  diagA[e*mesh->Np+cnt] = (pow(mesh->D[i*mesh->Nq+i],2) +
  				   pow(mesh->D[j*mesh->Nq+j],2) +
  				   pow(mesh->D[k*mesh->Nq+k],2) +
  				   lambda)*JW;
  	  ++cnt;
  	}
        }
      }
    }

    precon->o_diagA = mesh->device.malloc(Ntotal*sizeof(dfloat), diagA);

    // sum up
    meshParallelGatherScatter(mesh, ogs, precon->o_diagA, precon->o_diagA, dfloatString, "add");

    occaTimerTic(mesh->device,"CoarsePreconditionerSetup");
    // do this here since we need A*x
    ellipticCoarsePreconditionerSetupHex3D(mesh, precon, ogs, lambda, options);
    occaTimerToc(mesh->device,"CoarsePreconditionerSetup");

    
    free(diagA);
  } else if (strstr(options,"FULLALMOND")) {
    int nnz;
    nonZero_t *A;
    hgs_t *hgs;

    int Nnum = mesh->Np*mesh->Nelements;
    int *globalStarts = (int*) calloc(size+1, sizeof(int));

    if (strstr(options,"IPDG")) {

      MPI_Allgather(&(mesh->Nelements), 1, MPI_INT, globalStarts+1, 1, MPI_INT, MPI_COMM_WORLD);
      for(int r=0;r<size;++r)
        globalStarts[r+1] = globalStarts[r]+globalStarts[r+1]*mesh->Np;

      ellipticBuildIpdgHex3D(mesh, lambda, &A, &nnz,options);
    
      qsort(A, nnz, sizeof(nonZero_t), parallelCompareRowColumn);

    } else if (strstr(options,"CONTINUOUS")) {
      
      ellipticBuildContinuousHex3D(mesh,lambda,&A,&nnz,&hgs,globalStarts, options);
    }


    //collect global assembled matrix
    int *globalnnz       = (int *) calloc(size  ,sizeof(int));
    int *globalnnzOffset = (int *) calloc(size+1,sizeof(int));
    MPI_Allgather(&nnz, 1, MPI_INT, 
                  globalnnz, 1, MPI_INT, MPI_COMM_WORLD);
    globalnnzOffset[0] = 0;
    for (int n=0;n<size;n++)
      globalnnzOffset[n+1] = globalnnzOffset[n]+globalnnz[n];

    int globalnnzTotal = globalnnzOffset[size];

    int *globalRecvCounts  = (int *) calloc(size,sizeof(int));
    int *globalRecvOffsets = (int *) calloc(size,sizeof(int));
    for (int n=0;n<size;n++){
      globalRecvCounts[n] = globalnnz[n]*sizeof(nonZero_t);
      globalRecvOffsets[n] = globalnnzOffset[n]*sizeof(nonZero_t);
    }
    nonZero_t *globalA = (nonZero_t*) calloc(globalnnzTotal, sizeof(nonZero_t));

    MPI_Allgatherv(A, nnz*sizeof(nonZero_t), MPI_CHAR, 
                  globalA, globalRecvCounts, globalRecvOffsets, MPI_CHAR, MPI_COMM_WORLD);
    
    int *globalIndex = (int *) calloc(globalnnzTotal, sizeof(int));
    int *globalRows = (int *) calloc(globalnnzTotal, sizeof(int));
    int *globalCols = (int *) calloc(globalnnzTotal, sizeof(int));
    dfloat *globalVals = (dfloat*) calloc(globalnnzTotal,sizeof(dfloat));

    for (int n=0;n<globalnnzTotal;n++) {
      globalRows[n] = globalA[n].row;
      globalCols[n] = globalA[n].col;
      globalVals[n] = globalA[n].val;
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
