#include "advectionQuad3D.h"

void rk6(solver_t *solver) {
  int Nrk = 7;
  
  solver->rkA = (dfloat *) calloc(Nrk*Nrk,sizeof(dfloat));
  solver->rkE = (dfloat *) calloc(Nrk,sizeof(dfloat));
  solver->rkC = (dfloat *) calloc(Nrk,sizeof(dfloat));    
  
  dfloat rkC[7] = {0.0,
		   0.2,
		   0.3,
		   0.8,
		   8.0/9.0,
		   1.0,
		   1.0};

  dfloat rkA[7*7] ={
                           0.0,             0.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                           0.2,             0.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                      3.0/40.0,        9.0/40.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                     44.0/45.0,      -56.0/15.0,       32.0/9.0,          0.0,             0.0,       0.0, 0.0,
                19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0,             0.0,       0.0, 0.0,
                 9017.0/3168.0,     -355.0/33.0, 46732.0/5247.0,   49.0/176.0, -5103.0/18656.0,       0.0, 0.0, 
                    35.0/384.0,             0.0,   500.0/1113.0,  125.0/192.0,  -2187.0/6784.0, 11.0/84.0, 0.0};

  dfloat rkE[7] = {     71.0/57600.0,
			         0.0,
		       -71.0/16695.0,
	           	 71.0/1920.0,
		   -17253.0/339200.0,
		 	  22.0/525.0,
			   -1.0/40.0  };
    solver->Nrk = Nrk;
    memcpy(solver->rkA, rkA, Nrk*Nrk*sizeof(dfloat));
    memcpy(solver->rkE, rkE, Nrk*sizeof(dfloat));
    memcpy(solver->rkC, rkC, Nrk*sizeof(dfloat));
}

void advectionSetupDOPRIQuad3D (solver_t *solver) {

  mesh_t *mesh = solver->mesh;
  
  rk6(solver);

  solver->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*solver->Nfields*solver->Nrk,
				  sizeof(dfloat));

  solver->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*solver->Nfields*solver->Nrk,
				  sizeof(dfloat));
  
  dfloat *q_zero = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*solver->Nfields,
				    sizeof(dfloat));
  
  dfloat hmin = 1e9;
  for(iint e=0;e<mesh->Nelements;++e){
    dfloat lmin = 1e9, lmax = 0;
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
      }
    }
  }

  // need to change cfl and defn of dt
  dfloat cfl = 0.5; // depends on the stability region size
  
  solver->dt = cfl*hmin/((mesh->N+1.)*(mesh->N+1.)*solver->sqrtRT);
  //
  solver->finalTime = 5;
  
  //create dopri blocks
  solver->blockSize = 256;
  solver->Ntotal = mesh->Nelements*mesh->Np*solver->Nfields;
  solver->Nblock = (solver->Ntotal+solver->blockSize-1)/solver->blockSize;
  
  solver->errtmp = (dfloat *) calloc(solver->Nblock,sizeof(dfloat));
  
  solver->dtmin = 1E-7; 
  solver->absTol = 1E-9;
  solver->relTol = 1E-7;
  solver->safety = 0.9;
  
  //error control parameters
  solver->beta = 0.05;
  solver->factor1 = 0.2;
  solver->factor2 = 10.0;
  
  solver->exp1 = 0.2 - 0.75*solver->beta;
  solver->invfactor1 = 1.0/solver->factor1;
  solver->invfactor2 = 1.0/solver->factor2;
  solver->oldFactor = 1E-4;
  
  // hard code this for the moment
  solver->outputInterval = 2.5;
  solver->nextOutputTime = solver->outputInterval;
  solver->outputNumber = 0;
  
  //initial time
  solver->time = 0.0;
  solver->tstep = 0;
  solver->allStep = 0;

  occa::kernelInfo kernelInfo;
  
  //common device initialization
  advectionSetupOccaQuad3D(solver, &kernelInfo);
    
  solver->o_q =
    solver->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*solver->Nfields*sizeof(dfloat), solver->q);

  solver->o_qCorr =
    solver->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*solver->Nfields*sizeof(dfloat), q_zero);

  solver->o_qFilter =
    solver->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*solver->Nfields*sizeof(dfloat), solver->q);

  solver->o_qFiltered =
    solver->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*solver->Nfields*sizeof(dfloat), solver->q);
  
  solver->o_rhsq =
    solver->device.malloc(mesh->Np*mesh->Nelements*solver->Nfields*sizeof(dfloat), solver->rhsq);

  solver->o_resq =
    solver->device.malloc(mesh->Np*mesh->Nelements*solver->Nfields*solver->Nrk*sizeof(dfloat), solver->resq);

  solver->o_rkA = solver->device.malloc(solver->Nrk*solver->Nrk*sizeof(dfloat), solver->rkA);
  solver->o_rkE = solver->device.malloc(solver->Nrk*sizeof(dfloat), solver->rkE);

  solver->o_rkq =
    solver->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*solver->Nfields*sizeof(dfloat),solver->q);
  
  solver->o_rkerr =
    solver->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*solver->Nfields*sizeof(dfloat),q_zero);
  
  solver->o_errtmp = solver->device.malloc(solver->Nblock*sizeof(dfloat), solver->errtmp);
  
  solver->volumeKernel =
    solver->device.buildKernelFromSource(DHOLMES "/okl/advectionVolumeQuad3D.okl",
					 "advectionVolumeLSERKQuad3D",
					 kernelInfo);
  
  solver->volumeCorrectionKernel =
    solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolumeCorrectionQuad3D.okl",
				       "boltzmannVolumeCorrectionDOPRIQuad3D",
				       kernelInfo);
  solver->surfaceKernel =
    solver->device.buildKernelFromSource(DHOLMES "/okl/advectionSurfaceQuad3D.okl",
				       "advectionSurfaceLSERKQuad3D",
				       kernelInfo);
  solver->updateKernel =
    solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdateQuad3D.okl",
				       "boltzmannDOPRIUpdateQuad3D",
				       kernelInfo);
  solver->rkStageKernel =
    solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdateQuad3D.okl",
					 "boltzmannDOPRIrkStageQuad3D",
					 kernelInfo);
  solver->rkErrorEstimateKernel = 
    solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdateQuad3D.okl",
					 "boltzmannDOPRIerrorEstimateQuad3D",
					 kernelInfo);
  solver->filterKernelH =
    solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterHQuad3D.okl",
				       "boltzmannFilterHq0Quad3D",
				       kernelInfo);
  solver->filterKernelV =
    solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterVQuad3D.okl",
				       "boltzmannFilterVq0Quad3D",
				       kernelInfo);
}
