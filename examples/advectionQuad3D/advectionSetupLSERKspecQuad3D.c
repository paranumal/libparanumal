#include "advectionQuad3D.h"

void rk4_lserkspec_coeffs(solver_t *solver) {
  int Nrk = 5;
  
  solver->rka = (dfloat *) calloc(Nrk,sizeof(dfloat));
  solver->rkb = (dfloat *) calloc(Nrk,sizeof(dfloat));
  solver->rkc = (dfloat *) calloc(Nrk+1,sizeof(dfloat));
  
  dfloat rka[5] = {0.0,
		   -567301805773.0/1357537059087.0 ,
		   -2404267990393.0/2016746695238.0 ,
		   -3550918686646.0/2091501179385.0  ,
		   -1275806237668.0/842570457699.0};
  dfloat rkb[5] = { 1432997174477.0/9575080441755.0 ,
		    5161836677717.0/13612068292357.0 ,
		    1720146321549.0/2090206949498.0  ,
		    3134564353537.0/4481467310338.0  ,
		    2277821191437.0/14882151754819.0};
  dfloat rkc[6] = {0.0  ,
		   1432997174477.0/9575080441755.0 ,
		   2526269341429.0/6820363962896.0 ,
		   2006345519317.0/3224310063776.0 ,
		   2802321613138.0/2924317926251.0,
		   1.};
  solver->Nrk = Nrk;
  memcpy(solver->rka, rka, Nrk*sizeof(dfloat));
  memcpy(solver->rkb, rkb, Nrk*sizeof(dfloat));
  memcpy(solver->rkc, rkc, (Nrk+1)*sizeof(dfloat));
}

void advectionSetupLSERKspecQuad3D (solver_t *solver) {
  
  mesh_t *mesh = solver->mesh;
  
  solver->rhsq = (dfloat*) calloc(mesh->NgridElements*mesh->Np*solver->Nfields,
				  sizeof(dfloat));
  
  rk4_lserkspec_coeffs(solver);
  
  solver->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*solver->Nfields,
				  sizeof(dfloat));

    // set time step
    dfloat hmin = 1e9, hmax = 0;
    for(iint e=0;e<mesh->Nelements;++e){

      for(iint f=0;f<mesh->Nfaces;++f){
	for(iint n=0;n<mesh->Nfp;++n){
	  iint sid = mesh->Nsgeo*mesh->Nfp*mesh->Nfaces*e + mesh->Nsgeo*mesh->Nfp*f+n;

	  dfloat sJ   = mesh->sgeo[sid + mesh->Nq*SJID];
	  dfloat invJ = mesh->sgeo[sid + mesh->Nq*IJID];

	  // A = 0.5*h*L
	  // => J*2 = 0.5*h*sJ*2
	  // => h = 2*J/sJ

	  dfloat hest = 2./(sJ*invJ);

	  hmin = mymin(hmin, hest);
	  hmax = mymax(hmax, hest);
	}
      }
    }
    
    dfloat cfl = 1; // depends on the stability region size

    // dt ~ cfl (h/(N+1)^2)/(Lambda^2*fastest wave speed)
    solver->dt = cfl*hmin/((mesh->N+1.)*(mesh->N+1.));

    //dt = mymin(dt, cfl/mesh->tauInv);
    
    solver->finalTime = 5;
    solver->NtimeSteps = solver->finalTime/solver->dt;
    solver->dt = solver->finalTime/solver->NtimeSteps;
        
    printf("cfl = %g \n", cfl);
    printf("dt = %g\n", solver->dt);
    printf("max wave speed = %g\n", sqrt(3.)*solver->sqrtRT);
    
    // errorStep
    solver->errorStep = 100*mesh->Nq;
    
    printf("dt = %g\n", solver->dt);
    
    occa::kernelInfo kernelInfo;

    // resize some variables prior to loading
    solver->q = (dfloat *)realloc(solver->q,mesh->Np*mesh->NgridElements*solver->Nfields*sizeof(dfloat));
    for (iint i = 0; i <  mesh->NgridElements - mesh->Nelements; ++i) {
      solver->q[(mesh->Nelements + i)*mesh->Np*solver->Nfields] = 0;
    }
    
    advectionSetupOccaQuad3D(solver,&kernelInfo);

    dfloat *test = (dfloat *) calloc(mesh->Np*mesh->NgridElements*solver->Nfields,sizeof(dfloat));

    solver->o_weakD = solver->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D);

    solver->o_mass = solver->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->weakD);
    
    solver->o_q =
	solver->device.malloc(mesh->Np*mesh->NgridElements*solver->Nfields*sizeof(dfloat));
    solver->o_qs =
	solver->device.malloc(mesh->Np*mesh->NgridElements*solver->Nfields*sizeof(dfloat));
    solver->o_qw =
	solver->device.malloc(mesh->Np*mesh->NgridElements*solver->Nfields*sizeof(dfloat));

    solver->o_qpre =
      solver->device.malloc(mesh->Np*mesh->NgridElements*solver->Nfields*sizeof(dfloat), solver->q);

    solver->o_rhsqs =
	solver->device.malloc(mesh->Np*mesh->NgridElements*solver->Nfields*sizeof(dfloat));

    solver->o_rhsqw =
	solver->device.malloc(mesh->Np*mesh->NgridElements*solver->Nfields*sizeof(dfloat));

    solver->o_eInterp =
      solver->device.malloc(mesh->Np*(mesh->NgridElements-mesh->Nelements+1)*sizeof(iint),mesh->eInterp);

    solver->o_gridToE =
      solver->device.malloc(mesh->Nfaces*mesh->Nelements*sizeof(iint),mesh->gridToE);
    
    solver->o_overlapDirection =
      solver->device.malloc((mesh->NgridElements - mesh->Nelements + 1)*sizeof(char),mesh->overlapDirection);

    solver->o_rlocal =
      solver->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),mesh->rlocal);

    solver->o_slocal =
      solver->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),mesh->slocal);

    solver->o_par_loc =
      solver->device.malloc((mesh->NgridElements - mesh->Nelements + 1)*mesh->Np*sizeof(dfloat),mesh->par_loc);

    solver->o_perp_index = solver->device.malloc((mesh->NgridElements - mesh->Nelements + 1)*mesh->Np*sizeof(iint),mesh->perp_index);
    
    solver->o_qFilter =
	solver->device.malloc(mesh->NgridElements*solver->Nfields*mesh->Np*sizeof(dfloat));

    solver->o_qFilters =
	solver->device.malloc(mesh->NgridElements*solver->Nfields*mesh->Np*sizeof(dfloat));

    solver->o_qFilterw =
	solver->device.malloc(mesh->NgridElements*solver->Nfields*mesh->Np*sizeof(dfloat));
    
    solver->o_qCorr =
	solver->device.malloc(mesh->Nelements*solver->Nfields*mesh->Np*sizeof(dfloat));
    solver->o_resq =
      solver->device.malloc(mesh->Np*mesh->Nelements*solver->Nfields*sizeof(dfloat), solver->resq);

    solver->o_invmass =
	solver->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->inv_mass);
    solver->volumeKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/advectionVolumeQuad3D.okl",
					 "advectionVolumeLSERKmixedQuad3D",
					 kernelInfo);
    solver->volumeCorrectionKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolumeCorrectionQuad3D.okl",
					 "boltzmannVolumeCorrectionDOPRIQuad3D",
					 kernelInfo);
    solver->massMatrixKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolumeCorrectionQuad3D.okl",
					 "boltzmannMassMatrixQuad3D",
					 kernelInfo);
    solver->surfaceKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/advectionSurfaceQuad3D.okl",
					 "advectionSurfaceLSERKmixedQuad3D",
					 kernelInfo);
    solver->loadFilterGridKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannLoadFilterGridQuad3D.okl",
					   "boltzmannLoadFilterGridQuad3D",
					   kernelInfo);
    solver->updateKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdateQuad3D.okl",
					   "boltzmannLSERKsymUpdateQuad3D",
					   kernelInfo);
    solver->filterKernelH =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterHQuad3D.okl",
					 "boltzmannFilterHq0Quad3D",
					 kernelInfo);
    solver->filterKernelV =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterVQuad3D.okl",
					 "boltzmannFilterVq0Quad3D",
					 kernelInfo);
    solver->filterWeakKernelH =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterHQuad3D.okl",
					 "boltzmannFilterHtransQuad3D",
					 kernelInfo);
    solver->filterWeakKernelV =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterVQuad3D.okl",
					 "boltzmannFilterVtransQuad3D",
					 kernelInfo);
}
