#include "ellipticTet3D.h"

typedef struct{

  iint row;
  iint col;
  iint ownerRank;
  dfloat val;

}nonZero_t;

// compare on global indices 
int parallelCompareRowColumn(const void *a, const void *b){

  nonZero_t *fa = (nonZero_t*) a;
  nonZero_t *fb = (nonZero_t*) b;
  
  if(fa->row < fb->row) return -1;
  if(fa->row > fb->row) return +1;

  if(fa->col < fb->col) return -1;
  if(fa->col > fb->col) return +1;

  return 0;

}


void ellipticCoarsePreconditionerSetupTet3D(mesh_t *mesh, precon_t *precon, dfloat lambda, const char *options){

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  // ------------------------------------------------------------------------------------
  // 1. Create a contiguous numbering system, starting from the element-vertex connectivity
  iint Nnum = mesh->Nverts*mesh->Nelements;

  iint *globalNumbering = (iint*) calloc(Nnum, sizeof(iint));

  iint *globalOwners = (iint*) calloc(Nnum, sizeof(iint));
  iint *globalStarts = (iint*) calloc(size+1, sizeof(iint));

  // use original vertex numbering
  memcpy(globalNumbering, mesh->EToV, Nnum*sizeof(iint));

  // squeeze numbering
  meshParallelConsecutiveGlobalNumbering(Nnum, globalNumbering, globalOwners, globalStarts);
  
  //use the ordering to define a gather+scatter for assembly
  hgs_t *hgs = meshParallelGatherSetup(mesh, Nnum, globalNumbering, globalOwners);

  // temporary
  precon->o_ztmp = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat));
  
  // ------------------------------------------------------------------------------------
  // 2. Build coarse grid element basis functions

  dfloat *V1  = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  dfloat *Vr1 = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  dfloat *Vs1 = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  dfloat *Vt1 = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));

  for(iint n=0;n<mesh->Np;++n){   
    dfloat rn = mesh->r[n];
    dfloat sn = mesh->s[n];
    dfloat tn = mesh->t[n];
    
    V1[0*mesh->Np+n] = -0.5*(1.0+rn+sn+tn);
    V1[1*mesh->Np+n] = +0.5*(1.+rn);
    V1[2*mesh->Np+n] = +0.5*(1.+sn);
    V1[3*mesh->Np+n] = +0.5*(1.+tn);


    Vr1[0*mesh->Np+n] = 0.5*(-1);
    Vr1[1*mesh->Np+n] = 0.5*(+1);
    Vr1[2*mesh->Np+n] = 0;
    Vr1[3*mesh->Np+n] = 0;
      
    Vs1[0*mesh->Np+n] = 0.5*(-1);
    Vs1[1*mesh->Np+n] = 0.;
    Vs1[2*mesh->Np+n] = 0.5*(+1);
    Vs1[3*mesh->Np+n] = 0.;

    Vt1[0*mesh->Np+n] = 0.5*(-1);
    Vt1[1*mesh->Np+n] = 0.;
    Vt1[2*mesh->Np+n] = 0.;
    Vt1[3*mesh->Np+n] = 0.5*(+1);
  }
  precon->o_V1  = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), V1);
  precon->o_Vr1 = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), Vr1);
  precon->o_Vs1 = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), Vs1);
  precon->o_Vt1 = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), Vt1);

  // ------------------------------------------------------------------------------------
  // 3. Build non-zeros of stiffness matrix (unassembled)
  iint nnz = mesh->Nverts*mesh->Nverts*mesh->Nelements;
  iint   *rowsA = (iint*) calloc(nnz, sizeof(iint));
  iint   *colsA = (iint*) calloc(nnz, sizeof(iint));
  dfloat *valsA = (dfloat*) calloc(nnz, sizeof(dfloat));

  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnz, sizeof(nonZero_t));
  iint *AsendCounts  = (iint*) calloc(size, sizeof(iint));
  iint *ArecvCounts  = (iint*) calloc(size, sizeof(iint));
  iint *AsendOffsets = (iint*) calloc(size+1, sizeof(iint));
  iint *ArecvOffsets = (iint*) calloc(size+1, sizeof(iint));
    
  iint cnt = 0;

  dfloat *cV1  = (dfloat*) calloc(mesh->cubNp*mesh->Nverts, sizeof(dfloat));
  dfloat *cVr1 = (dfloat*) calloc(mesh->cubNp*mesh->Nverts, sizeof(dfloat));
  dfloat *cVs1 = (dfloat*) calloc(mesh->cubNp*mesh->Nverts, sizeof(dfloat));
  dfloat *cVt1 = (dfloat*) calloc(mesh->cubNp*mesh->Nverts, sizeof(dfloat));

  for(iint n=0;n<mesh->cubNp;++n){   
    dfloat rn = mesh->cubr[n];
    dfloat sn = mesh->cubs[n];
    dfloat tn = mesh->cubt[n];
    
    cV1[0*mesh->cubNp+n] = -0.5*(1.0+rn+sn+tn);
    cV1[1*mesh->cubNp+n] = +0.5*(1.+rn);
    cV1[2*mesh->cubNp+n] = +0.5*(1.+sn);
    cV1[3*mesh->cubNp+n] = +0.5*(1.+tn);

    cVr1[0*mesh->cubNp+n] = 0.5*(-1);
    cVr1[1*mesh->cubNp+n] = 0.5*(+1);
    cVr1[2*mesh->cubNp+n] = 0;
    cVr1[3*mesh->cubNp+n] = 0;
      
    cVs1[0*mesh->cubNp+n] = 0.5*(-1);
    cVs1[1*mesh->cubNp+n] = 0;
    cVs1[2*mesh->cubNp+n] = 0.5*(+1);
    cVs1[3*mesh->cubNp+n] = 0;

    cVt1[0*mesh->cubNp+n] = 0.5*(-1);
    cVt1[1*mesh->cubNp+n] = 0;
    cVt1[2*mesh->cubNp+n] = 0;
    cVt1[3*mesh->cubNp+n] = 0.5*(+1);
  }

  dfloat Srr1[4][4], Srs1[4][4], Srt1[4][4];
  dfloat Ssr1[4][4], Sss1[4][4], Sst1[4][4];
  dfloat Str1[4][4], Sts1[4][4], Stt1[4][4]; 
  dfloat MM1[4][4];
  for(iint n=0;n<mesh->Nverts;++n){
    for(iint m=0;m<mesh->Nverts;++m){
      Srr1[n][m] = 0;
      Srs1[n][m] = 0;
      Srt1[n][m] = 0;
      Ssr1[n][m] = 0;
      Sss1[n][m] = 0;
      Sst1[n][m] = 0;
      Str1[n][m] = 0;
      Sts1[n][m] = 0;
      Stt1[n][m] = 0;
      MM1[n][m] = 0;
      
      for(iint i=0;i<mesh->cubNp;++i){
      	iint idn = n*mesh->cubNp+i;
      	iint idm = m*mesh->cubNp+i;
      	dfloat cw = mesh->cubw[i];
      	Srr1[n][m] += cw*(cVr1[idn]*cVr1[idm]);
      	Srs1[n][m] += cw*(cVr1[idn]*cVs1[idm]);
        Srt1[n][m] += cw*(cVr1[idn]*cVt1[idm]);
      	Ssr1[n][m] += cw*(cVs1[idn]*cVr1[idm]);
        Sss1[n][m] += cw*(cVs1[idn]*cVs1[idm]);
        Sst1[n][m] += cw*(cVs1[idn]*cVt1[idm]);
        Str1[n][m] += cw*(cVt1[idn]*cVr1[idm]);
        Sts1[n][m] += cw*(cVt1[idn]*cVs1[idm]);
        Stt1[n][m] += cw*(cVt1[idn]*cVt1[idm]);
      	MM1[n][m] += cw*(cV1[idn]*cV1[idm]);
      }
    }
  }	

  
  printf("Building coarse matrix system\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Nverts;++n){
      for(iint m=0;m<mesh->Nverts;++m){
      	dfloat Snm = 0;

      	dfloat rx = mesh->vgeo[e*mesh->Nvgeo + RXID];
      	dfloat sx = mesh->vgeo[e*mesh->Nvgeo + SXID];
        dfloat tx = mesh->vgeo[e*mesh->Nvgeo + TXID];
      	dfloat ry = mesh->vgeo[e*mesh->Nvgeo + RYID];
      	dfloat sy = mesh->vgeo[e*mesh->Nvgeo + SYID];
        dfloat ty = mesh->vgeo[e*mesh->Nvgeo + TYID];
        dfloat rz = mesh->vgeo[e*mesh->Nvgeo + RZID];
        dfloat sz = mesh->vgeo[e*mesh->Nvgeo + SZID];
        dfloat tz = mesh->vgeo[e*mesh->Nvgeo + TZID];
      	dfloat J  = mesh->vgeo[e*mesh->Nvgeo +  JID];
      	
      	Snm  = J*(rx*rx+ry*ry+rz*rz)*Srr1[n][m];
      	Snm += J*(rx*sx+ry*sy+rz*sz)*Srs1[n][m];
        Snm += J*(rx*tx+ry*ty+rz*tz)*Srt1[n][m];
        Snm += J*(sx*rx+sy*ry+sz*rz)*Ssr1[n][m];
      	Snm += J*(sx*sx+sy*sy+sz*sz)*Sss1[n][m];
        Snm += J*(sx*tx+sy*ty+sz*tz)*Sst1[n][m];
        Snm += J*(tx*rx+ty*ry+tz*rz)*Str1[n][m];
        Snm += J*(tx*sx+ty*sy+tz*sz)*Sts1[n][m];
        Snm += J*(tx*tx+ty*ty+tz*tz)*Stt1[n][m];
      	Snm += J*lambda*MM1[n][m];

        dfloat nonZeroThreshold = 1e-7;
        if(fabs(Snm)>nonZeroThreshold) {
        	valsA[cnt] = Snm;
        	rowsA[cnt] = e*mesh->Nverts+n;
        	colsA[cnt] = e*mesh->Nverts+m;

        	// pack non-zero
        	sendNonZeros[cnt].val = Snm;
        	sendNonZeros[cnt].row = globalNumbering[e*mesh->Nverts+n];
        	sendNonZeros[cnt].col = globalNumbering[e*mesh->Nverts+m];
        	sendNonZeros[cnt].ownerRank = globalOwners[e*mesh->Nverts+n];
        	
        	++cnt;
        }
      }
    }
  }

  // count how many non-zeros to send to each process
  for(iint n=0;n<cnt;++n)
    AsendCounts[sendNonZeros[n].ownerRank] += sizeof(nonZero_t);

  // sort by row ordering
  qsort(sendNonZeros, cnt, sizeof(nonZero_t), parallelCompareRowColumn);
  
  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(AsendCounts, 1, MPI_IINT, ArecvCounts, 1, MPI_IINT, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  iint recvNtotal = 0;
  for(iint r=0;r<size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    recvNtotal += ArecvCounts[r]/sizeof(nonZero_t);
  }

  nonZero_t *recvNonZeros = (nonZero_t*) calloc(recvNtotal, sizeof(nonZero_t));
  
  // determine number to receive
  MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_CHAR,
		recvNonZeros, ArecvCounts, ArecvOffsets, MPI_CHAR,
		MPI_COMM_WORLD);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort(recvNonZeros, recvNtotal, sizeof(nonZero_t), parallelCompareRowColumn);

  // compress duplicates
  cnt = 0;
  for(iint n=1;n<recvNtotal;++n){
    if(recvNonZeros[n].row == recvNonZeros[cnt].row &&
       recvNonZeros[n].col == recvNonZeros[cnt].col){
      recvNonZeros[cnt].val += recvNonZeros[n].val;
    }
    else{
      ++cnt;
      recvNonZeros[cnt] = recvNonZeros[n];
    }
  }
  recvNtotal = cnt+1;

  iint *recvRows = (iint *) calloc(recvNtotal,sizeof(iint));
  iint *recvCols = (iint *) calloc(recvNtotal,sizeof(iint));
  dfloat *recvVals = (dfloat *) calloc(recvNtotal,sizeof(dfloat));
  
  for (iint n=0;n<recvNtotal;n++) {
    recvRows[n] = recvNonZeros[n].row;
    recvCols[n] = recvNonZeros[n].col;
    recvVals[n] = recvNonZeros[n].val;
  }
  
  //collect global assembled matrix
  iint globalnnz[size];
  iint globalnnzOffset[size+1];
  MPI_Allgather(&recvNtotal, 1, MPI_IINT, 
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

  MPI_Allgatherv(recvNonZeros, recvNtotal*sizeof(nonZero_t), MPI_CHAR, 
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

  
  printf("Done building coarse matrix system\n");
  if(strstr(options, "XXT")){
    precon->xxt = xxtSetup(Nnum,
			   globalNumbering,
			   nnz,
			   rowsA,
			   colsA,
			   valsA,
			   0,
			   iintString,
			   dfloatString); // 0 if no null space
  }

  if(strstr(options, "ALMOND")){
 
    precon->parAlmond = parAlmondSetup(mesh,
         Nnum, 
         globalStarts,
         recvNtotal,      
         recvRows,        
         recvCols,       
         recvVals,   
         0, // 0 if no null space
         hgs,
         options); 
    
  }

  precon->o_r1 = mesh->device.malloc(Nnum*sizeof(dfloat));
  precon->o_z1 = mesh->device.malloc(Nnum*sizeof(dfloat));
  precon->r1 = (dfloat*) malloc(Nnum*sizeof(dfloat));
  precon->z1 = (dfloat*) malloc(Nnum*sizeof(dfloat));
  
  free(AsendCounts);
  free(ArecvCounts);
  free(AsendOffsets);
  free(ArecvOffsets);

}
