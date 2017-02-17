#include "ellipticTri2D.h"

void ellipticCoarsePreconditionerSetupTri2D(mesh_t *mesh, precon_t *precon, dfloat lambda){

  // ------------------------------------------------------------------------------------
  // 1. Create a contiguous numbering system, starting from the element-vertex connectivity
  iint Nnum = mesh->Nverts*mesh->Nelements;
  
  iint *globalNumbering = (iint*) calloc(Nnum, sizeof(iint));

  // use original vertex numbering
  memcpy(globalNumbering, mesh->EToV, Nnum*sizeof(iint));

  // build gs
  void *gsh = gsParallelGatherScatterSetup(Nnum, globalNumbering);

  dfloat *degree = (dfloat*) calloc(Nnum, sizeof(dfloat));
  for(iint n=0;n<Nnum;++n)
    degree[n] = 1;
  
  gsParallelGatherScatter(gsh, degree, dfloatString, "add");
  
  dfloat *invDegree = (dfloat*) calloc(Nnum, sizeof(dfloat));
  for(iint n=0;n<Nnum;++n)
    invDegree[n] = 1./degree[n];

  precon->o_coarseInvDegree = mesh->device.malloc(Nnum*sizeof(dfloat), invDegree);

  // clean up
  gsParallelGatherScatterDestroy(gsh);

  // temporary
  precon->o_ztmp = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat));
  
  // ------------------------------------------------------------------------------------
  // 2. Build coarse grid element basis functions

  dfloat *V1  = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  dfloat *Vr1 = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  dfloat *Vs1 = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));

  for(iint n=0;n<mesh->Np;++n){   
    dfloat rn = mesh->r[n];
    dfloat sn = mesh->s[n];
    
    V1[0*mesh->Np+n] = -0.5*(rn+sn);
    V1[1*mesh->Np+n] = +0.5*(1.+rn);
    V1[2*mesh->Np+n] = +0.5*(1.+sn);

    Vr1[0*mesh->Np+n] = 0.5*(-1);
    Vr1[1*mesh->Np+n] = 0.5*(+1);
    Vr1[2*mesh->Np+n] = 0;
      
    Vs1[0*mesh->Np+n] = 0.5*(-1);
    Vs1[1*mesh->Np+n] = 0;
    Vs1[2*mesh->Np+n] = 0.5*(+1);
  }
  precon->o_V1  = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), V1);
  precon->o_Vr1 = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), Vr1);
  precon->o_Vs1 = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), Vs1);

  // ------------------------------------------------------------------------------------
  // 3. Build non-zeros of stiffness matrix (unassembled)
  iint nnz = mesh->Nverts*mesh->Nverts*mesh->Nelements;
  iint   *rowsA = (iint*) calloc(nnz, sizeof(iint));
  iint   *colsA = (iint*) calloc(nnz, sizeof(iint));
  dfloat *valsA = (dfloat*) calloc(nnz, sizeof(dfloat));
  
  iint cnt = 0;

  dfloat *cV1  = (dfloat*) calloc(mesh->cubNp*mesh->Nverts, sizeof(dfloat));
  dfloat *cVr1 = (dfloat*) calloc(mesh->cubNp*mesh->Nverts, sizeof(dfloat));
  dfloat *cVs1 = (dfloat*) calloc(mesh->cubNp*mesh->Nverts, sizeof(dfloat));

  for(iint n=0;n<mesh->cubNp;++n){   
    dfloat rn = mesh->cubr[n];
    dfloat sn = mesh->cubs[n];
    
    cV1[0*mesh->cubNp+n] = -0.5*(rn+sn);
    cV1[1*mesh->cubNp+n] = +0.5*(1.+rn);
    cV1[2*mesh->cubNp+n] = +0.5*(1.+sn);

    cVr1[0*mesh->cubNp+n] = 0.5*(-1);
    cVr1[1*mesh->cubNp+n] = 0.5*(+1);
    cVr1[2*mesh->cubNp+n] = 0;
      
    cVs1[0*mesh->cubNp+n] = 0.5*(-1);
    cVs1[1*mesh->cubNp+n] = 0;
    cVs1[2*mesh->cubNp+n] = 0.5*(+1);
  }

  dfloat Srr1[3][3], Srs1[3][3], Sss1[3][3], MM1[3][3];
  for(iint n=0;n<mesh->Nverts;++n){
    for(iint m=0;m<mesh->Nverts;++m){
      Srr1[n][m] = 0;
      Srs1[n][m] = 0;
      Sss1[n][m] = 0;
      MM1[n][m] = 0;
      
      for(iint i=0;i<mesh->cubNp;++i){
	iint idn = n*mesh->cubNp+i;
	iint idm = m*mesh->cubNp+i;
	dfloat cw = mesh->cubw[i];
	Srr1[n][m] += cw*(cVr1[idn]*cVr1[idm]);
	Srs1[n][m] += cw*(cVr1[idn]*cVs1[idm]);
	Sss1[n][m] += cw*(cVs1[idn]*cVs1[idm]);
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
	dfloat ry = mesh->vgeo[e*mesh->Nvgeo + RYID];
	dfloat sy = mesh->vgeo[e*mesh->Nvgeo + SYID];
	dfloat J  = mesh->vgeo[e*mesh->Nvgeo +  JID];
	
	Snm  = J*(rx*rx+ry*ry)*Srr1[n][m];
	Snm += J*(rx*sx+ry*sy)*Srs1[n][m];
	Snm += J*(sx*sx+sy*sy)*Sss1[n][m];
	Snm += J*lambda*MM1[n][m];

	valsA[cnt] = Snm;
	rowsA[cnt] = e*mesh->Nverts+n;
	colsA[cnt] = e*mesh->Nverts+m;
	++cnt;
      }
    }
  }

  printf("Done building coarse matrix system\n");
  precon->xxt = xxtSetup(Nnum,
			 globalNumbering,
			 nnz,
			 rowsA,
			 colsA,
			 valsA,
			 0,
			 iintString,
			 dfloatString); // 0 if no null space
  
  precon->o_r1 = mesh->device.malloc(Nnum*sizeof(dfloat));
  precon->o_z1 = mesh->device.malloc(Nnum*sizeof(dfloat));
  precon->r1 = (dfloat*) malloc(Nnum*sizeof(dfloat));
  precon->z1 = (dfloat*) malloc(Nnum*sizeof(dfloat));
  
}
