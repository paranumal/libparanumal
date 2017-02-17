#include "ellipticQuad2D.h"

void ellipticCoarsePreconditionerSetupQuad2D(mesh_t *mesh, precon_t *precon, dfloat lambda){

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

  for(iint j=0;j<mesh->Nq;++j){
    for(iint i=0;i<mesh->Nq;++i){
      iint n = i+j*mesh->Nq;
      /*
      dfloat rn = mesh->r[n];
      dfloat sn = mesh->s[n];
      */
      dfloat rn = mesh->gllz[i];
      dfloat sn = mesh->gllz[j];
      V1[0*mesh->Np+n] = 0.25*(1-rn)*(1-sn);
      V1[1*mesh->Np+n] = 0.25*(1+rn)*(1-sn);
      V1[2*mesh->Np+n] = 0.25*(1+rn)*(1+sn);
      V1[3*mesh->Np+n] = 0.25*(1-rn)*(1+sn);

      Vr1[0*mesh->Np+n] = 0.25*(-1)*(1-sn);
      Vr1[1*mesh->Np+n] = 0.25*(+1)*(1-sn);
      Vr1[2*mesh->Np+n] = 0.25*(+1)*(1+sn);
      Vr1[3*mesh->Np+n] = 0.25*(-1)*(1+sn);
      
      Vs1[0*mesh->Np+n] = 0.25*(1-rn)*(-1);
      Vs1[1*mesh->Np+n] = 0.25*(1+rn)*(-1);
      Vs1[2*mesh->Np+n] = 0.25*(1+rn)*(+1);
      Vs1[3*mesh->Np+n] = 0.25*(1-rn)*(+1);
    }
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

  printf("Building coarse matrix system\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Nverts;++n){
      for(iint m=0;m<mesh->Nverts;++m){
	dfloat Snm = 0;
 
	// use GLL nodes for integration
	// (since Jacobian is high order tensor-product polynomial)
	for(iint j=0;j<mesh->Nq;++j){
	  for(iint i=0;i<mesh->Nq;++i){
	    iint id = i+j*mesh->Nq;
	    
	    dfloat Vr1ni = Vr1[n*mesh->Np+id];
	    dfloat Vs1ni = Vs1[n*mesh->Np+id];
	    dfloat V1ni  = V1[n*mesh->Np+id];
	    
	    dfloat Vr1mi = Vr1[m*mesh->Np+id];
	    dfloat Vs1mi = Vs1[m*mesh->Np+id];
	    dfloat V1mi  = V1[m*mesh->Np+id];

	    dfloat rx = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + RXID*mesh->Np];
	    dfloat sx = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + SXID*mesh->Np];
	    dfloat ry = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + RYID*mesh->Np];
	    dfloat sy = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + SYID*mesh->Np];
	    dfloat JW = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + JWID*mesh->Np];

	    dfloat Vx1ni = rx*Vr1ni+sx*Vs1ni;
	    dfloat Vy1ni = ry*Vr1ni+sy*Vs1ni;
	    dfloat Vx1mi = rx*Vr1mi+sx*Vs1mi;
	    dfloat Vy1mi = ry*Vr1mi+sy*Vs1mi;
	    
	    Snm += (Vx1ni*Vx1mi+Vy1ni*Vy1mi)*JW;
	    Snm += (lambda*V1ni*V1mi)*JW;
	  }
	}
	//	Snm = (n==m) ? 1: 0;

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
