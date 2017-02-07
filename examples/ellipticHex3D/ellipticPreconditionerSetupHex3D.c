#include "ellipticHex3D.h"

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

precon_t *ellipticPreconditionerSetupHex3D(mesh3D *mesh, ogs_t *ogs, dfloat lambda){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // assumes meshParallelGatherScatterSetup3D has been called
  
  // ????? need to extend storage for halo ?????
  iint Nlocal = mesh->Np*mesh->Nelements;
  iint Nhalo  = mesh->Np*mesh->totalHaloPairs;
  iint Ntrace = mesh->Nfp*mesh->Nfaces*mesh->Nelements;

  // offsets to extract second layer
  iint *offset = (iint*) calloc(mesh->Nfaces, sizeof(iint));
  offset[0] = +mesh->Nfp;
  offset[1] = +mesh->Nq;
  offset[2] = -1;
  offset[3] = -mesh->Nq;
  offset[4] = +1;
  offset[5] = -mesh->Nfp;

  // build gather-scatter
  iint NqP = mesh->Nq+2;
  iint NpP = NqP*NqP*NqP;

  iint *offsetP = (iint*) calloc(mesh->Nfaces, sizeof(iint));
  offsetP[0] = +NqP*NqP;
  offsetP[1] = +NqP;
  offsetP[2] = -1;
  offsetP[3] = -NqP;
  offsetP[4] = +1;
  offsetP[5] = -NqP*NqP;

  iint *faceNodesPrecon = (iint*) calloc(mesh->Nfp*mesh->Nfaces, sizeof(iint));
  for(iint j=0;j<mesh->Nq;++j)
    for(iint i=0;i<mesh->Nq;++i)
      faceNodesPrecon[i+j*mesh->Nq+0*mesh->Nfp] = i+1 + (j+1)*NqP + 0*NqP*NqP;
  
  for(iint k=0;k<mesh->Nq;++k)
    for(iint i=0;i<mesh->Nq;++i)
      faceNodesPrecon[i+k*mesh->Nq+1*mesh->Nfp] = i+1 + 0*NqP + (k+1)*NqP*NqP;
  
  for(iint k=0;k<mesh->Nq;++k)
    for(iint j=0;j<mesh->Nq;++j)
      faceNodesPrecon[j+k*mesh->Nq+2*mesh->Nfp] = NqP-1 + (j+1)*NqP + (k+1)*NqP*NqP;

  for(iint k=0;k<mesh->Nq;++k)
    for(iint i=0;i<mesh->Nq;++i)
      faceNodesPrecon[i+k*mesh->Nq+3*mesh->Nfp] = i+1 + (NqP-1)*NqP + (k+1)*NqP*NqP;

  for(iint k=0;k<mesh->Nq;++k)
    for(iint j=0;j<mesh->Nq;++j)
      faceNodesPrecon[j+k*mesh->Nq+4*mesh->Nfp] = 0 + (j+1)*NqP + (k+1)*NqP*NqP;
  
  for(iint j=0;j<mesh->Nq;++j)
    for(iint i=0;i<mesh->Nq;++i)
      faceNodesPrecon[i+j*mesh->Nq+5*mesh->Nfp] = i+1 + (j+1)*NqP + (NqP-1)*NqP*NqP;

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
    for(iint k=0;k<mesh->Nq;++k){
      for(iint j=0;j<mesh->Nq;++j){
	for(iint i=0;i<mesh->Nq;++i){
	  iint id  = i + j*mesh->Nq + k*mesh->Nq*mesh->Nq + e*mesh->Np;
	  iint pid = i + 1 + (j+1)*NqP + (k+1)*NqP*NqP + e*NpP;

	  preconGatherInfo[pid] = gatherInfo[id];
	}
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
	for(iint k=0;k<mesh->NqP;++k){
	  for(iint j=0;j<mesh->NqP;++j){
	    for(iint i=0;i<mesh->NqP;++i){
	      iint id = i + mesh->NqP*j + mesh->NqP*mesh->NqP*k+ e*NpP;
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
  iint skip = 0;
  while(preconGatherInfo[skip].baseId==0 && skip<NpP*mesh->Nelements){
    ++skip;
  }
  printf("skip = %d out of %d\n", skip, NpP*mesh->Nelements);

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
  precon_t *precon = (precon_t*) calloc(1, sizeof(precon_t));
  precon->ogsP = meshParallelGatherScatterSetup3D(mesh,
						  NlocalP,
						  sizeof(dfloat),
						  gatherLocalIdsP,
						  gatherBaseIdsP,
						  gatherHaloFlagsP);

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
    for(iint k=0;k<mesh->Nq;++k){
      for(iint j=0;j<mesh->Nq;++j){
	for(iint i=0;i<mesh->Nq;++i){
	  iint id  = i + j*mesh->Nq + k*mesh->Nq*mesh->Nq + e*mesh->Np;
	  iint pid = (i+1) + (j+1)*NqP + (k+1)*NqP*NqP + e*NpP;

	  // all patch interior nodes are local
	  preconGatherInfoDg[pid].baseId = localNums[id];
	}
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
  precon->ogsDg = meshParallelGatherScatterSetup3D(mesh,
						   NlocalDg,
						   sizeof(dfloat),
						   gatherLocalIdsDg,
						   gatherBaseIdsDg,
						   gatherHaloFlagsDg);
  
  // -------------------------------------------------------------------------------------------
  

  
  precon->o_faceNodesP = mesh->device.malloc(mesh->Nfp*mesh->Nfaces*sizeof(iint), faceNodesPrecon);
  precon->o_vmapPP     = mesh->device.malloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements*sizeof(iint), vmapPP);

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
    
    dfloat Jhrinv2 = 0, Jhsinv2 = 0, Jhtinv2 = 0, J = 0;
    for(iint n=0;n<mesh->Np;++n){
      dfloat W = mesh->gllw[n%mesh->Nq]*
	mesh->gllw[(n/mesh->Nq)%mesh->Nq]*
	mesh->gllw[n/(mesh->Nq*mesh->Nq)];
      iint base = mesh->Nggeo*mesh->Np*e + n;

      J = mymax(J, mesh->ggeo[base + mesh->Np*GWJID]/W);
      Jhrinv2 = mymax(Jhrinv2, mesh->ggeo[base + mesh->Np*G00ID]/W);
      Jhsinv2 = mymax(Jhsinv2, mesh->ggeo[base + mesh->Np*G11ID]/W);
      Jhtinv2 = mymax(Jhtinv2, mesh->ggeo[base + mesh->Np*G22ID]/W);
      
    }
    
    for(iint k=0;k<NqP;++k){
      for(iint j=0;j<NqP;++j){
	for(iint i=0;i<NqP;++i){
	  iint pid = i + j*NqP + k*NqP*NqP + e*NpP;
	  
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

  iint Ntotal = mesh->Np*mesh->Nelements;
  dfloat *diagA = (dfloat*) calloc(Ntotal, sizeof(dfloat));
				   
  for(iint e=0;e<mesh->Nelements;++e){
    iint cnt = 0;
    for(iint k=0;k<mesh->Nq;++k){
      for(iint j=0;j<mesh->Nq;++j){
	for(iint i=0;i<mesh->Nq;++i){
	  
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
  meshParallelGatherScatter3D(mesh, ogs, precon->o_diagA, precon->o_diagA, dfloatString, "add");

  if(Nhalo){
    dfloat *vgeoSendBuffer = (dfloat*) calloc(Nhalo*mesh->Nvgeo, sizeof(dfloat));
    
    // import geometric factors from halo elements
    mesh->vgeo = (dfloat*) realloc(mesh->vgeo, (Nlocal+Nhalo)*mesh->Nvgeo*sizeof(dfloat));
    
    meshHaloExchange(mesh,
		     mesh->Nvgeo*mesh->Np*sizeof(dfloat),
		     mesh->vgeo,
		     vgeoSendBuffer,
		     mesh->vgeo + Nlocal*mesh->Nvgeo);
    
    mesh->o_vgeo =
      mesh->device.malloc((Nlocal+Nhalo)*mesh->Nvgeo*sizeof(dfloat), mesh->vgeo);
  }

  
  free(diagA);
  
  return precon;
}
