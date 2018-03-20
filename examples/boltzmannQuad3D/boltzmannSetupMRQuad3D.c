#include "boltzmannQuad3D.h"

void brownMinion(dfloat bmRho, dfloat bmDelta, dfloat sphereRadius,
		 dfloat x, dfloat y, dfloat z,
		 dfloat *u, dfloat *v, dfloat *w){

  dfloat Utangential = 0.25*(1+tanh(bmRho*(-z+0.5)))*(1+tanh(bmRho*(0.5+z)));

  dfloat uout, vout;

  if(x*x+y*y>1e-4) {
    uout =  -y*Utangential/(x*x+y*y);
    vout =   x*Utangential/(x*x+y*y);
  }
  else{
    uout = 0;
    vout = 0;
  }

  dfloat wout = bmDelta*sin(2*atan2(y,x))*(1-z*z);

  dfloat udotx = uout*x+vout*y+wout*z;
  *u = uout - udotx*x/(sphereRadius*sphereRadius);
  *v = vout - udotx*y/(sphereRadius*sphereRadius);
  *w = wout - udotx*z/(sphereRadius*sphereRadius);

}


void rk_coeffs(mesh_t *mesh) {

  iint Nlevels = mesh->MRABNlevels;

  const iint Nr = 32;
  dfloat complex R[Nr];

  for(iint ind =1; ind <= Nr; ++ind){
    const dfloat theta = (dfloat) (ind - 0.5) / (dfloat) Nr;
    R[ind-1] = cexp(I*M_PI* theta);
  }

  mesh->MRSAAB_A = (dfloat *) calloc(3*3*Nlevels,sizeof(dfloat));
  mesh->MRSAAB_B = (dfloat *) calloc(3*3*Nlevels,sizeof(dfloat));
  mesh->MRSAAB_C = (dfloat *) calloc(    Nlevels,sizeof(dfloat));
  mesh->MRAB_A   = (dfloat *) calloc(3*3*Nlevels,sizeof(dfloat));
  mesh->MRAB_B   = (dfloat *) calloc(3*3*Nlevels,sizeof(dfloat));
  mesh->MRAB_C   = (dfloat *) calloc(    Nlevels,sizeof(dfloat));

  iint MRABorder = mesh->Nrhs;

  for(iint l = 0; l<Nlevels; ++l){
    // MRSAAB coefficients
    dfloat alpha = -mesh->tauInv*mesh->dt*pow(2,l);
    dfloat h  = mesh->dt * pow(2,l);
    //
    for (iint order=0; order<3; ++order){
      // computation of coefficients based on magnitude
      const iint id = order*Nlevels*3 + l*3;
      if(order==0){

	double complex a1 = 0. + 0.* I;
	double complex b1 = 0. + 0.* I;

	for(iint i = 0; i<Nr; ++i ){
	  double complex lr = alpha  + R[i];
	  a1 +=  h*(cexp(lr) - 1.)/lr;
	  b1 +=  h*(cexp(lr/2.) - 1.)/lr;
	}
	// Full dt coeeficients
	mesh->MRSAAB_A[id + 0] = creal(a1)/Nr;
	mesh->MRSAAB_A[id + 1] = 0.f;
	mesh->MRSAAB_A[id + 2] = 0.f;
	// Half coefficients
	mesh->MRSAAB_B[id + 0] = creal(b1)/Nr;
	mesh->MRSAAB_B[id + 1] = 0.f;
	mesh->MRSAAB_B[id + 2] = 0.f;

	// MRAB coefficients
	mesh->MRAB_A[id + 0]   =  h ;
	mesh->MRAB_A[id + 1]   =  0.f ;
	mesh->MRAB_A[id + 2]   =  0.f ;

	mesh->MRAB_B[id+0]     =  h/2. ;
	mesh->MRAB_B[id+1]     =  0.f ;
	mesh->MRAB_B[id+2]     =  0.f ;
      }

      else if(order==1){

	double complex a1 = 0. + 0.* I;
	double complex b1 = 0. + 0.* I;
	double complex a2 = 0. + 0.* I;
	double complex b2 = 0. + 0.* I;

	for(iint i = 0; i<Nr; ++i ){
	  double complex lr = alpha  + R[i];
	  a1 +=  h*(-2.*lr + (1.+lr)*cexp(lr) - 1.)/cpow(lr,2);
	  a2 +=  h*(lr - cexp(lr) + 1.)/cpow(lr,2);
	  b1 +=  h*(-1.5*lr + (1.+lr)*cexp(lr/2.) - 1.)/cpow(lr,2);
	  b2 +=  h*(0.5*lr - cexp(lr/2.) + 1.)/cpow(lr,2);
	}
	// Full dt coeeficients
	mesh->MRSAAB_A[id + 0] = creal(a1)/Nr;
	mesh->MRSAAB_A[id + 1] = creal(a2)/Nr;
	mesh->MRSAAB_A[id + 2] = 0.f;
	// Half coefficients
	mesh->MRSAAB_B[id + 0] = creal(b1)/Nr;
	mesh->MRSAAB_B[id + 1] = creal(b2)/Nr;
	mesh->MRSAAB_B[id + 2] = 0.f;


	// MRAB coefficients
	mesh->MRAB_A[id + 0]   =  3.*h/2. ;
	mesh->MRAB_A[id + 1]   = -1.*h/2. ;
	mesh->MRAB_A[id + 2]   =  0.f ;

	mesh->MRAB_B[id + 0]   =  5.*h/8. ;
	mesh->MRAB_B[id + 1]   = -1.*h/8. ;
	mesh->MRAB_B[id + 2]   =   0.f ;
      }

      else{
	double complex a1 = 0. + 0.* I;
	double complex b1 = 0. + 0.* I;
	double complex a2 = 0. + 0.* I;
	double complex b2 = 0. + 0.* I;
	double complex a3 = 0. + 0.* I;
	double complex b3 = 0. + 0.* I;

	for(iint i = 0; i<Nr; ++i ){
	  double complex lr = alpha  + R[i];
	  a1 += h*(-2.5*lr - 3.*cpow(lr,2) + (1.+cpow(lr,2)+1.5*lr)*cexp(lr) - 1.)/cpow(lr,3);
	  a2 += h*(4.*lr + 3.*cpow(lr,2)- (2.*lr + 2.0)*cexp(lr) + 2.)/cpow(lr,3);
	  a3 +=-h*(1.5*lr + cpow(lr,2)- (0.5*lr + 1.)*cexp(lr) + 1.)/cpow(lr,3);
	  b1 += h*(cexp(lr/2.)- 2.*lr - (15.*cpow(lr,2))/8.f + cpow(lr,2)*cexp(lr/2.) + 3.*lr*cexp(lr/2.)/2. - 1.)/cpow(lr,3);
	  b2 += h*(3.*lr - 2.*cexp(lr/2.0) + 1.25*cpow(lr,2) - 2.*lr*cexp(lr/2.) + 2.)/cpow(lr,3);
	  b3 +=-h*(lr - cexp(lr/2.) + 0.375*cpow(lr,2) - 0.5*lr*cexp(lr/2.) + 1.)/cpow(lr,3);
	}


	// Full dt coeeficients
	mesh->MRSAAB_A[id+0] = creal(a1)/Nr;
	mesh->MRSAAB_A[id+1] = creal(a2)/Nr;
	mesh->MRSAAB_A[id+2] = creal(a3)/Nr;
	// Half coefficients
	mesh->MRSAAB_B[id+0] = creal(b1)/Nr;
	mesh->MRSAAB_B[id+1] = creal(b2)/Nr;
	mesh->MRSAAB_B[id+2] = creal(b3)/Nr;

	// MRAB coefficients
	mesh->MRAB_A[id+0]   =  23.*h/12. ;
	mesh->MRAB_A[id+1]   = -16.*h/12. ;
	mesh->MRAB_A[id+2]   =  5. *h/12. ;

	mesh->MRAB_B[id+0]   =  17.*h/24. ;
	mesh->MRAB_B[id+1]   = - 7.*h/24. ;
	mesh->MRAB_B[id+2]   =   2.*h/24. ;


      }
    }

    // Exponential part
    mesh->MRSAAB_C[l]    = exp(alpha);
    mesh->MRAB_C[l]      =   h ;
  }
}

solver_t *boltzmannSetupMRQuad3D(mesh_t *mesh){
  
  solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));

  solver->mesh = mesh;
  
  mesh->Nrhs = 3; //hardcoded order of multirate solver  
  
  mesh->Nfields = 10;
  
  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Nrhs*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  mesh->MRABshiftIndex = (iint*) calloc(mesh->MRABNlevels,sizeof(iint));
  mesh->shift_init = (iint*) calloc(mesh->MRABNlevels,sizeof(iint));

  // set temperature, gas constant, wave speeds
  mesh->RT = 9.;
  mesh->sqrtRT = sqrt(mesh->RT);
  dfloat nu = 4e-4; 8.9e-5; 9.8e-4;  // was 6.e-3

  // initial conditions
  dfloat sR = mesh->sphereRadius;

  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){

      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];
      dfloat z = mesh->z[n + mesh->Np*e];

      // Brown Minion shear layer roll up
      dfloat bmRho = 40; // was 40
      dfloat bmDelta  = 0.05;

      dfloat rho = 1;

      dfloat umod, vmod, wmod;

      brownMinion(bmRho, bmDelta, sR, x, y, z, &umod, &vmod, &wmod);

      dfloat delta = 1e-5;

      dfloat uP, uM, vP, vM, wP, wM;

      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x+delta, y, z, &uP, &vP, &wP);
      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x-delta, y, z, &uM, &vM, &wM);

      dfloat dudx = (uP-uM)/(2*delta);
      dfloat dvdx = (vP-vM)/(2*delta);
      dfloat dwdx = (wP-wM)/(2*delta);

      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x, y+delta, z, &uP, &vP, &wP);
      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x, y-delta, z, &uM, &vM, &wM);

      dfloat dudy = (uP-uM)/(2*delta);
      dfloat dvdy = (vP-vM)/(2*delta);
      dfloat dwdy = (wP-wM)/(2*delta);

      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x, y, z+delta, &uP, &vP, &wP);
      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x, y, z-delta, &uM, &vM, &wM);

      dfloat dudz = (uP-uM)/(2*delta);
      dfloat dvdz = (vP-vM)/(2*delta);
      dfloat dwdz = (wP-wM)/(2*delta);

      dfloat divu = dudx + dvdy + dwdz;

#if 1
      dfloat sigma11 = nu*(dudx+dudx - (2*divu/3));
      dfloat sigma12 = nu*(dvdx+dudy);
      dfloat sigma13 = nu*(dwdx+dudz);
      dfloat sigma22 = nu*(dvdy+dvdy - (2*divu/3));
      dfloat sigma23 = nu*(dwdy+dvdz);
      dfloat sigma33 = nu*(dwdz+dwdz - (2*divu/3));
#else
      dfloat sigma11 = 0;
      dfloat sigma12 = 0;
      dfloat sigma13 = 0;
      dfloat sigma22 = 0;
      dfloat sigma23 = 0;
      dfloat sigma33 = 0;
#endif
      dfloat q1bar = rho;
      dfloat q2bar = rho*umod/mesh->sqrtRT;
      dfloat q3bar = rho*vmod/mesh->sqrtRT;
      dfloat q4bar = rho*wmod/mesh->sqrtRT;
      dfloat q5bar = (rho*umod*umod - sigma11)/(sqrt(2.)*mesh->RT);
      dfloat q6bar = (rho*vmod*vmod - sigma22)/(sqrt(2.)*mesh->RT);
      dfloat q7bar = (rho*wmod*wmod - sigma33)/(sqrt(2.)*mesh->RT);
      dfloat q8bar  = (rho*umod*vmod - sigma12)/mesh->RT;
      dfloat q9bar =  (rho*umod*wmod - sigma13)/mesh->RT;
      dfloat q10bar = (rho*vmod*wmod - sigma23)/mesh->RT;

      dfloat t = 0;

      int base = n + e*mesh->Np*mesh->Nfields;

#if 0
      // Gaussian pulse centered at X=(1,0,0)
      q1bar = 1 + .1*exp(-20*((x-1)*(x-1)+y*y+z*z));
#endif


      mesh->q[base+0*mesh->Np] = q1bar; // uniform density, zero flow

      mesh->q[base+1*mesh->Np] = q2bar;
      mesh->q[base+2*mesh->Np] = q3bar;
      mesh->q[base+3*mesh->Np] = q4bar;

      mesh->q[base+4*mesh->Np] = q5bar;
      mesh->q[base+5*mesh->Np] = q6bar;
      mesh->q[base+6*mesh->Np] = q7bar;

      mesh->q[base+7*mesh->Np] = q8bar;
      mesh->q[base+8*mesh->Np] = q9bar;
     mesh->q[base+9*mesh->Np] = q10bar;
    }
  }
  // set BGK collision relaxation rate
  // nu = R*T*tau
  // 1/tau = RT/nu
  //  dfloat nu = 1.e-2/.5;
  //  dfloat nu = 1.e-3/.5;
  //  dfloat nu = 5.e-4;
  //    dfloat nu = 1.e-2; TW works for start up fence
  dfloat cfl_small = 0.2; // depends on the stability region size (was .4, then 2)
  dfloat cfl_large = 4*cfl_small;
  
  mesh->localdt = (dfloat *) calloc(mesh->Nelements,sizeof(dfloat));
  
  mesh->tauInv = mesh->RT/nu; // TW


  dfloat glmin = 1e9, glmax = -1e9;
  // set time step
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

	lmin = mymin(lmin, hest);
	lmax = mymax(lmax, hest);
      }
    }
    if (mesh->cubeDistance[e] == 0) {
      mesh->localdt[e]  = cfl_small*lmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*mesh->sqrtRT);
    }
    else {
      mesh->localdt[e]  = cfl_large*lmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*mesh->sqrtRT);
    }

    glmin = mymin(glmin, lmin);
    glmax = mymax(glmax, lmax);
    
  }

  //dt = mymin(dt, cfl/mesh->tauInv);

  mesh->finalTime = 10;
  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  
  iint maxLevels=100;
  meshMRABSetupQuad3D(mesh,mesh->localdt,maxLevels);

  dfloat dt = mesh->dt;

  mesh->shiftIndex=0;

  mesh->lev_updates = (iint *) calloc(mesh->MRABNlevels,sizeof(iint));

  printf("cfl = %g %g\n", cfl_small,cfl_large);
  printf("dt = %g\n", dt);
  printf("max wave speed = %g\n", sqrt(3.)*mesh->sqrtRT);
  printf("global h in [%g,%g]\n", glmin, glmax);
  
  // errorStep
  mesh->errorStep = 100*mesh->Nq;

  printf("dt = %g\n", mesh->dt);

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", (rank+1)%2);
  //sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 0");
  //sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");

  boltzmannPlotVTUQuad3D(mesh, "bah.vtu", 0);
  
  occa::kernelInfo kernelInfo;

  // fixed to set up quad info on device too
  boltzmannOccaSetupQuad3D(mesh, deviceConfig,  kernelInfo);

  // quad stuff

  kernelInfo.addDefine("p_Nq", mesh->Nq);
  kernelInfo.addDefine("p_nrhs",mesh->Nrhs);

  mesh->o_D  = mesh->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D);
  mesh->o_vgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Np*mesh->Nvgeo*sizeof(dfloat),
			mesh->vgeo);
  
  //initialization isn't strictly necessary here.
  mesh->o_qFilter =
    mesh->device.malloc(mesh->Nrhs*mesh->Nelements*mesh->Nfields*mesh->Np*sizeof(dfloat),mesh->rhsq);

  mesh->o_qFiltered =
    mesh->device.malloc(mesh->Nrhs*mesh->Nelements*mesh->Nfields*mesh->Np*sizeof(dfloat),mesh->rhsq);

  mesh->o_qCorr =
    mesh->device.malloc(mesh->Nrhs*mesh->Nelements*mesh->Nfields*mesh->Np*sizeof(dfloat),mesh->rhsq);

  mesh->o_qPreCorr =
    mesh->device.malloc(mesh->Nrhs*mesh->Nelements*mesh->Nfields*mesh->Np*sizeof(dfloat),mesh->fQM);

  mesh->o_prerhsq =
    mesh->device.malloc(mesh->Nelements*mesh->Nfields*mesh->Np*sizeof(dfloat),mesh->fQM);

  mesh->o_qPreFilter =
    mesh->device.malloc(mesh->Nelements*mesh->Nfields*mesh->Np*sizeof(dfloat),mesh->fQM);

  mesh->o_qPreFiltered =
    mesh->device.malloc(mesh->Nelements*mesh->Nfields*mesh->Np*sizeof(dfloat),mesh->fQM);
  
  mesh->o_sgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*mesh->Nsgeo*sizeof(dfloat),
			mesh->sgeo);

  mesh->o_MRABlevels = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*sizeof(iint),mesh->MRABlevel);
  mesh->o_MRABelementIds = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
  mesh->o_MRABhaloIds = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
  mesh->o_fQM  = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*mesh->Nfields*sizeof(dfloat),mesh->fQM);
  mesh->o_lev_updates = mesh->device.malloc(mesh->MRABNlevels*sizeof(iint),mesh->lev_updates);
  mesh->o_shift = mesh->device.malloc(mesh->MRABNlevels*sizeof(iint),mesh->MRABshiftIndex);
  
  printf("start\n");
  mesh->o_dualProjMatrix =
    mesh->device.malloc(mesh->Nq*mesh->Nq*3*sizeof(dfloat),mesh->dualProjMatrix);

  mesh->o_cubeFaceNumber =
    mesh->device.malloc(mesh->Nelements*sizeof(iint),mesh->cubeFaceNumber);

  mesh->o_EToE =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),mesh->EToE);
  printf("end\n");
  for (iint lev = 0; lev < mesh->MRABNlevels;lev++) {
    if (mesh->MRABNelements[lev])
      mesh->o_MRABelementIds[lev]
	= mesh->device.malloc(mesh->MRABNelements[lev]*sizeof(iint),
			      mesh->MRABelementIds[lev]);
    if (mesh->MRABNhaloElements[lev])
      mesh->o_MRABhaloIds[lev]
	= mesh->device.malloc(mesh->MRABNhaloElements[lev]*sizeof(iint),
			      mesh->MRABhaloIds[lev]);
  }

  //break out computation of rk coefficients
  //needs to be last so we have mrab levels
  rk_coeffs(mesh);

      
  // specialization for Boltzmann

  kernelInfo.addDefine("p_maxNodesVolume", mymax(mesh->cubNp,mesh->Np));

  int StoVmaxNodes = mymax(mesh->Np,mesh->Nfp*mesh->Nfaces);
  kernelInfo.addDefine("p_StoVmaxNodes",StoVmaxNodes);
  
  int maxNodes = mesh->Nfp;
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 256/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 256/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  // physics
  kernelInfo.addDefine("p_sqrtRT", mesh->sqrtRT);
  kernelInfo.addDefine("p_invsqrtRT", (dfloat)(1./mesh->sqrtRT));
  kernelInfo.addDefine("p_sqrt2", (dfloat)sqrt(2.));
  kernelInfo.addDefine("p_invsqrt2", (dfloat)sqrt(1./2.));
  kernelInfo.addDefine("p_isq12", (dfloat)sqrt(1./12.));
  kernelInfo.addDefine("p_isq6", (dfloat)sqrt(1./6.));
  kernelInfo.addDefine("p_tauInv", mesh->tauInv);

  kernelInfo.addDefine("p_invRadiusSq", 1./(mesh->sphereRadius*mesh->sphereRadius));

  kernelInfo.addDefine("p_fainv", (dfloat) 0.0); // turn off rotation

  mesh->volumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolumeQuad3D.okl",
				       "boltzmannVolumeSAQuad3D",
				       kernelInfo);

  mesh->volumePreKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolumeQuad3D.okl",
				       "boltzmannVolumeQuad3D",
				       kernelInfo);
  mesh->volumeCorrectionKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolumeCorrectionQuad3D.okl",
				       "boltzmannVolumeCorrectionSAQuad3D",
				       kernelInfo);
  mesh->volumeCorrPreKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolumeCorrectionQuad3D.okl",
				       "boltzmannVolumeCorrectionQuad3D",
				       kernelInfo);
  mesh->surfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurfaceQuad3D.okl",
				       "boltzmannSurfaceMRQuad3D",
				       kernelInfo);

  mesh->surfacePreKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurfaceQuad3D.okl",
				       "boltzmannSurfaceQuad3D",
				       kernelInfo);
  
  mesh->traceUpdateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdateQuad3D.okl",
				       "boltzmannMRSAABTraceUpdateQuad3D",
				       kernelInfo);
  
  mesh->updateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdateQuad3D.okl",
				       "boltzmannMRSAABUpdateQuad3D",
				       kernelInfo);

  mesh->updatePreKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdateQuad3D.okl",
				       "boltzmannLSERKUpdateQuad3D",
				       kernelInfo);

  mesh->traceUpdatePreKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdateQuad3D.okl",
				       "boltzmannLSERKTraceUpdateQuad3D",
				       kernelInfo);
  
  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
				       "meshHaloExtract2D",
				       kernelInfo);
  mesh->filterKernelH =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterHQuad3D.okl",
				       "boltzmannFilterHQuad3D",
				       kernelInfo);
  mesh->filterKernelV =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterVQuad3D.okl",
				       "boltzmannFilterVQuad3D",
				       kernelInfo);

  mesh->filterKernelHLSERK =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterHQuad3D.okl",
				       "boltzmannFilterHLSERKQuad3D",
				       kernelInfo);
  mesh->filterKernelVLSERK =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterVQuad3D.okl",
				       "boltzmannFilterVLSERKQuad3D",
				       kernelInfo);

  
  mesh->filterKernelq0H =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterHQuad3D.okl",
				       "boltzmannFilterHq0Quad3D",
				       kernelInfo);
  mesh->filterKernelq0V =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterVQuad3D.okl",
				       "boltzmannFilterVq0Quad3D",
				       kernelInfo);
  return solver;
}
