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
solver_t *boltzmannSetupQuad3D(mesh_t *mesh){

  solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));

  solver->mesh = mesh;

  mesh->Nfields = 10;

  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  mesh->fQM  = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*mesh->Nfields, sizeof(dfloat));
  mesh->fQP  = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*mesh->Nfields, sizeof(dfloat));

  // set temperature, gas constant, wave speeds
  mesh->RT = 9.;
  mesh->sqrtRT = sqrt(mesh->RT);
  dfloat nu = 2.e-4;  // was 6.e-3

  // initial conditions
  dfloat sR = mesh->sphereRadius;

  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){

      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];
      dfloat z = mesh->z[n + mesh->Np*e];

      // Brown Minion shear layer roll up
      dfloat bmRho = 40;
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

  mesh->tauInv = mesh->RT/nu; // TW

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

  dfloat cfl = 2; // depends on the stability region size (was .4)

  // dt ~ cfl (h/(N+1)^2)/(Lambda^2*fastest wave speed)
  dfloat dt = cfl*hmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*mesh->sqrtRT);

  dt = mymin(dt, cfl/mesh->tauInv);

  printf("hmin = %g\n", hmin);
  printf("hmax = %g\n", hmax);
  printf("cfl = %g\n", cfl);
  printf("dt = %g\n", dt);
  printf("max wave speed = %g\n", sqrt(3.)*mesh->sqrtRT);

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  //
  mesh->finalTime = 10;
  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  mesh->dt = mesh->finalTime/mesh->NtimeSteps;

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
  //  sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 1");
  //  sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");

  occa::kernelInfo kernelInfo;

  // fixed to set up quad info on device too
  boltzmannOccaSetupQuad3D(mesh, deviceConfig,  kernelInfo);

  // quad stuff

  kernelInfo.addDefine("p_Nq", mesh->Nq);

  printf("mesh->Nq = %d\n", mesh->Nq);
  mesh->o_D  = mesh->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D);

    mesh->o_vgeo =
      mesh->device.malloc(mesh->Nelements*mesh->Np*mesh->Nvgeo*sizeof(dfloat),
			  mesh->vgeo);

      mesh->o_sgeo =
	mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*mesh->Nsgeo*sizeof(dfloat),
			    mesh->sgeo);

      // specialization for Boltzmann

      kernelInfo.addDefine("p_maxNodesVolume", mymax(mesh->cubNp,mesh->Np));

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
					     "boltzmannVolumeQuad3D",
					     kernelInfo);
	printf("starting surface\n");
	  mesh->surfaceKernel =
	    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSurfaceQuad3D.okl",
					       "boltzmannSurfaceQuad3D",
					       kernelInfo);
	  printf("ending surface\n");

	    mesh->updateKernel =
	      mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdateQuad3D.okl",
						 "boltzmannUpdateQuad3D",
						 kernelInfo);

	      mesh->haloExtractKernel =
		mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
						   "meshHaloExtract2D",
						   kernelInfo);


	      return solver;
}
