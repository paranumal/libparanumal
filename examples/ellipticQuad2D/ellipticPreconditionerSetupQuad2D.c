#include "ellipticQuad2D.h"

typedef struct{

  iint localId;
  iint baseId;
  iint baseRank;
  iint haloFlag;
  
} preconGatherInfo_t;

int parallelCompareBaseRank(const void *a, const void *b){

  preconGatherInfo_t *fa = (preconGatherInfo_t*) a;
  preconGatherInfo_t *fb = (preconGatherInfo_t*) b;

  if(fa->baseRank < fb->baseRank) return -1;
  if(fa->baseRank > fb->baseRank) return +1;

  if(fa->baseId < fb->baseId) return -1;
  if(fa->baseId > fb->baseId) return +1;

  return 0;

}

precon_t *ellipticPreconditionerSetupQuad2D(mesh2D *mesh, ogs_t *ogs, dfloat lambda){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // assumes meshParallelGatherScatterSetup2D has been called
  
  // ????? need to extend storage for halo ?????
  iint Nlocal = mesh->Np*mesh->Nelements;
  iint Nhalo  = mesh->Np*mesh->totalHaloPairs;
  iint Ntrace = mesh->Nfp*mesh->Nfaces*mesh->Nelements;

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
  for(iint i=0;i<mesh->Nq;++i)
    faceNodesPrecon[i+0*mesh->Nfp] = i+1 + 0*NqP;
  
  for(iint j=0;j<mesh->Nq;++j)
    faceNodesPrecon[j+1*mesh->Nfp] = NqP-1 + (j+1)*NqP;
  
  for(iint i=0;i<mesh->Nq;++i)
    faceNodesPrecon[i+2*mesh->Nfp] = i+1 + (NqP-1)*NqP;

  for(iint j=0;j<mesh->Nq;++j)
    faceNodesPrecon[j+3*mesh->Nfp] = 0 + (j+1)*NqP;

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
  
  // space for gather base indices with halo
  preconGatherInfo_t *gatherInfo =
    (preconGatherInfo_t*) calloc(Nlocal+Nhalo, sizeof(preconGatherInfo_t));

  // rearrange in node order
  for(iint n=0;n<Nlocal;++n){
    iint id = mesh->gatherLocalIds[n];
    gatherInfo[id].baseId   = mesh->gatherBaseIds[n] + 1;
    gatherInfo[id].baseRank = mesh->gatherBaseRanks[n];
    gatherInfo[id].haloFlag = mesh->gatherHaloFlags[n];
  }

  if(Nhalo){
    // send buffer for outgoing halo
    preconGatherInfo_t *sendBuffer = (preconGatherInfo_t*) calloc(Nhalo, sizeof(preconGatherInfo_t));
    
    meshHaloExchange2D(mesh,
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

      if(gatherInfo[idM].haloFlag){
	preconGatherInfo[idMP].haloFlag = 1;
	preconGatherInfo[idMP+offsetP[f]].haloFlag = 1;
	preconGatherInfo[idMP+2*offsetP[f]].haloFlag = 1;
      }
    }
  }
  
  // reset local ids
  for(iint n=0;n<mesh->Nelements*NpP;++n)
    preconGatherInfo[n].localId = n;

  char fname[BUFSIZ];
  sprintf(fname, "haloFlag%05d.dat", rank);
  FILE *fp = fopen(fname, "w");
  
  for(iint p=0;p<size;++p){
    if(p==rank){
      for(iint e=0;e<mesh->Nelements;++e){
	fprintf(fp,"e=%d: \n", e);
	for(iint j=0;j<mesh->NqP;++j){
	  for(iint i=0;i<mesh->NqP;++i){
	    iint id = i + mesh->NqP*j + e*NpP;
	    fprintf(fp,"%d ", preconGatherInfo[id].haloFlag);
	  }
	  fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
      }
    }
  }
  fclose(fp);
  
  // sort by rank then base index
  qsort(preconGatherInfo, NpP*mesh->Nelements, sizeof(preconGatherInfo_t), parallelCompareBaseRank);
  
  // do not gather-scatter nodes labelled zero
  iint skip = 0;
  while(preconGatherInfo[skip].baseId==0 && skip<NpP*mesh->Nelements){
    ++skip;
  }
  printf("skip = %d out of %d\n", skip, NpP*mesh->Nelements);

  // reset local ids
  iint NlocalP = NpP*mesh->Nelements - skip;
  iint *gatherLocalIdsP  = (iint*) calloc(NlocalP, sizeof(iint));
  iint *gatherBaseIdsP   = (iint*) calloc(NlocalP, sizeof(iint));
  iint *gatherBaseRanksP = (iint*) calloc(NlocalP, sizeof(iint));
  iint *gatherHaloFlagsP = (iint*) calloc(NlocalP, sizeof(iint));
  for(iint n=0;n<NlocalP;++n){
    gatherLocalIdsP[n]  = preconGatherInfo[n+skip].localId;
    gatherBaseIdsP[n]   = preconGatherInfo[n+skip].baseId;
    gatherBaseRanksP[n] = preconGatherInfo[n+skip].baseRank;
    gatherHaloFlagsP[n] = preconGatherInfo[n+skip].haloFlag;
  }

  // make preconBaseIds => preconNumbering
  precon_t *precon = (precon_t*) calloc(1, sizeof(precon_t));
  precon->ogsP = meshParallelGatherScatterSetup2D(mesh,
						  NlocalP,
						  sizeof(dfloat),
						  gatherLocalIdsP,
						  gatherBaseIdsP,
						  gatherBaseRanksP,
						  gatherHaloFlagsP);

  precon->o_faceNodesP = mesh->device.malloc(mesh->Nfp*mesh->Nfaces*sizeof(iint), faceNodesPrecon);
  precon->o_vmapPP     = mesh->device.malloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements*sizeof(iint), vmapPP);
  precon->o_oasForward = mesh->device.malloc(NqP*NqP*sizeof(dfloat), mesh->oasForward);
  precon->o_oasBack    = mesh->device.malloc(NqP*NqP*sizeof(dfloat), mesh->oasBack);

  /// ---------------------------------------------------------------------------
  
  // hack estimate for Jacobian scaling
  dfloat *diagInvOp = (dfloat*) calloc(NpP*mesh->Nelements, sizeof(dfloat));
  for(iint e=0;e<mesh->Nelements;++e){

    // S = Jabc*(wa*wb*wc*lambda + wb*wc*Da'*wa*Da + wa*wc*Db'*wb*Db + wa*wb*Dc'*wc*Dc)
    // S = Jabc*wa*wb*wc*(lambda*I+1/wa*Da'*wa*Da + 1/wb*Db'*wb*Db + 1/wc*Dc'*wc*Dc)
    
    dfloat Jhrinv2 = 0, Jhsinv2 = 0, J = 0;
    for(iint n=0;n<mesh->Np;++n){

      dfloat W = mesh->gllz[n%mesh->Nq]*mesh->gllz[n/mesh->Nq];
      
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
      }
    }
  }
  
  precon->o_oasDiagInvOp = mesh->device.malloc(NpP*mesh->Nelements*sizeof(dfloat), diagInvOp);

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
  meshParallelGatherScatter2D(mesh, ogs, precon->o_diagA, precon->o_diagA, dfloatString);

  free(diagA);
  
  return precon;
}
