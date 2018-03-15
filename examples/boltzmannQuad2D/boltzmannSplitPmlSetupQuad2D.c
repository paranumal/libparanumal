#include "boltzmannQuad2D.h"

// NBN: toggle use of 2nd stream
#define USE_2_STREAMS
// #undef USE_2_STREAMS

void boltzmannSplitPmlSetupQuad2D(mesh2D *mesh){

  mesh->Nfields = 8;
  
  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  mesh->pmlqx    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->rhspmlqx = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				 sizeof(dfloat));
  mesh->respmlqx = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  mesh->pmlqy    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->rhspmlqy = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				 sizeof(dfloat));
  mesh->respmlqy = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  mesh->pmlNT    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				    sizeof(dfloat));
  mesh->rhspmlNT = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				    sizeof(dfloat));
  mesh->respmlNT = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				    sizeof(dfloat));
  
  mesh->sigmax= (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
  mesh->sigmay= (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
  
  // set temperature, gas constant, wave speeds
  mesh->RT = 9.;
  mesh->sqrtRT = sqrt(mesh->RT);  
  
  // initial conditions
  // uniform flow
  dfloat rho = 1, u = 1, v = 0; //u = 1.f/sqrt(2.f), v = 1.f/sqrt(2.f); 
  dfloat sigma11 = 0, sigma12 = 0, sigma22 = 0;
  //  dfloat ramp = 0.5*(1.f+tanh(10.f*(0-.5f)));
  dfloat ramp = boltzmannRampFunction2D(0);
  dfloat q1bar = rho;
  dfloat q2bar = rho*u/mesh->sqrtRT;
  dfloat q3bar = rho*v/mesh->sqrtRT;
  dfloat q4bar = (rho*u*v - sigma12)/mesh->RT;
  dfloat q5bar = (rho*u*u - sigma11)/(sqrt(2.)*mesh->RT);
  dfloat q6bar = (rho*v*v - sigma22)/(sqrt(2.)*mesh->RT);

  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      dfloat t = 0;
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];

#if 0
      boltzmannCavitySolution2D(x, y, t,
				mesh->q+cnt, mesh->q+cnt+1, mesh->q+cnt+2);
#endif

#if 0
      boltzmannGaussianPulse2D(x, y, t,
			       mesh->q+cnt,
			       mesh->q+cnt+1,
			       mesh->q+cnt+2,
			       mesh->q+cnt+3,
			       mesh->q+cnt+4,
			       mesh->q+cnt+5);
#endif
      mesh->q[cnt+0] = q1bar; // uniform density, zero flow
      mesh->q[cnt+1] = ramp*q2bar;
      mesh->q[cnt+2] = ramp*q3bar;
      mesh->q[cnt+3] = ramp*ramp*q4bar;
      mesh->q[cnt+4] = ramp*ramp*q5bar;
      mesh->q[cnt+5] = ramp*ramp*q6bar;
    
      cnt += mesh->Nfields;

      iint id = mesh->Np*mesh->Nfields*e + n;
      mesh->pmlqx[id+0*mesh->Np] = 0.f*q1bar;
      mesh->pmlqx[id+1*mesh->Np] = 0.f*q2bar;
      mesh->pmlqx[id+2*mesh->Np] = 0.f*q3bar;
      mesh->pmlqx[id+3*mesh->Np] = 0.f*q4bar;
      mesh->pmlqx[id+4*mesh->Np] = 0.f*q5bar;
      mesh->pmlqx[id+5*mesh->Np] = 0.f*q6bar;

      mesh->pmlqy[id+0*mesh->Np] = 0.f*q1bar;
      mesh->pmlqy[id+1*mesh->Np] = 0.f*q2bar;
      mesh->pmlqy[id+2*mesh->Np] = 0.f*q3bar;
      mesh->pmlqy[id+3*mesh->Np] = 0.f*q4bar;
      mesh->pmlqy[id+4*mesh->Np] = 0.f*q5bar;
      mesh->pmlqy[id+5*mesh->Np] = 0.f*q6bar;
    }
  }

  // set BGK collision relaxation rate
  // nu = R*T*tau
  // 1/tau = RT/nu
  //  dfloat nu = 1.e-2/.5;
  //  dfloat nu = 1.e-3/.5;
  //  dfloat nu = 5.e-4;
  //    dfloat nu = 1.e-2; TW works for start up fence
  dfloat nu = 2.e-3;  // was 6.e-3
  mesh->tauInv = mesh->RT/nu; // TW
  
  // set penalty parameter
  //  mesh->Lambda2 = 0.5/(sqrt(3.)*mesh->sqrtRT);
  mesh->Lambda2 = 0.5/(mesh->sqrtRT);

  // find elements with center inside PML zone
  dfloat xmin = -4, xmax = 8, ymin = -4, ymax = 4;
  dfloat xsigma = 80, ysigma = 80;
  //    dfloat xsigma = 0, ysigma = 0;
  
  for(iint e=0;e<mesh->Nelements;++e){
    dfloat cx = 0, cy = 0;
    for(iint n=0;n<mesh->Nverts;++n){
      cx += mesh->EX[e*mesh->Nverts+n];
      cy += mesh->EY[e*mesh->Nverts+n];
    }
    cx /= mesh->Nverts;
    cy /= mesh->Nverts;
    
    // add element outside [xmin,xmax]x[ymin,ymax] to pml
#if 0
    if(cx<xmin || cx>xmax)
      mesh->sigmax[e] = xsigma;
    if(cy<ymin || cy>ymax)
      mesh->sigmay[e] = ysigma;
#endif
    
    for(iint n=0;n<mesh->Np;++n){
      dfloat x = mesh->x[n + e*mesh->Np];
      dfloat y = mesh->y[n + e*mesh->Np];
      //      if(cx<xmax+1 && cx>xmin-1 && cy<ymax+1 && cy>ymin-1){
      {
	if(cx>xmax)
	  //  mesh->sigmax[mesh->Np*e + n] = xsigma;
	  mesh->sigmax[mesh->Np*e + n] = xsigma*pow(x-xmax,2);
	if(cx<xmin)
	  //  mesh->sigmax[mesh->Np*e + n] = xsigma;
	  mesh->sigmax[mesh->Np*e + n] = xsigma*pow(x-xmin,2);
	if(cy>ymax)
	  //	  mesh->sigmay[mesh->Np*e + n] = ysigma;
	  mesh->sigmay[mesh->Np*e + n] = ysigma*pow(y-ymax,2);
	if(cy<ymin)
	  //  mesh->sigmay[mesh->Np*e + n] = ysigma;
	  mesh->sigmay[mesh->Np*e + n] = ysigma*pow(y-ymin,2);
      }
    }
  }

  // set time step
  dfloat hmin = 1e9, hmax = 0;
  for(iint e=0;e<mesh->Nelements;++e){  

    for(iint f=0;f<mesh->Nfaces;++f){
      for(iint n=0;n<mesh->Nfp;++n){
	iint sid = mesh->Nsgeo*(mesh->Nfp*mesh->Nfaces*e + mesh->Nfp*f+n);
	dfloat sJ   = mesh->sgeo[sid + SJID];
	dfloat invJ = mesh->sgeo[sid + IJID];
	
	// A = 0.5*h*L
	// => J*2 = 0.5*h*sJ*2
	// => h = 2*J/sJ
	
	dfloat hest = 2./(sJ*invJ);
	
	hmin = mymin(hmin, hest);
	hmax = mymax(hmax, hest);
      }
    }
  }
    
  dfloat cfl = .8; // depends on the stability region size (was .4)

  // dt ~ cfl (h/(N+1)^2)/(Lambda^2*fastest wave speed)
  dfloat dt = cfl*hmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(3.)*mesh->sqrtRT);
  
  printf("hmin = %g\n", hmin);
  printf("hmax = %g\n", hmax);
  printf("cfl = %g\n", cfl);
  printf("dt = %g\n", dt);
  printf("max wave speed = %g\n", sqrt(3.)*mesh->sqrtRT);
  
  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  //
  mesh->finalTime = 50;
  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  mesh->dt = mesh->finalTime/mesh->NtimeSteps;

  // errorStep
  mesh->errorStep = 10000;

  printf("dt = %g\n", mesh->dt);

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", (rank+1)%3);
  //  sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 1");
  //  sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");	  

  occa::kernelInfo kernelInfo;

  // fixed to set up quad info on device too
  meshOccaSetup2D(mesh, deviceConfig,  kernelInfo);

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

  
  //pml variables
  mesh->o_pmlqx =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);
  mesh->o_rhspmlqx =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqx);
  mesh->o_respmlqx =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->respmlqx);

  mesh->o_pmlqy =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
  mesh->o_rhspmlqy =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlqy);
  mesh->o_respmlqy =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->respmlqy);

  mesh->o_pmlNT =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->pmlNT);
  mesh->o_rhspmlNT =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhspmlNT);
  mesh->o_respmlNT =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->respmlNT);

  
  mesh->o_sigmax =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), mesh->sigmax);

  mesh->o_sigmay =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), mesh->sigmay);

  // specialization for Boltzmann

  kernelInfo.addDefine("p_maxNodesVolume", mymax(mesh->cubNp,mesh->Np));
  
  kernelInfo.addDefine("p_pmlAlpha", (float).2);
  
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 512/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  // physics 
  kernelInfo.addDefine("p_Lambda2", 0.5f);
  kernelInfo.addDefine("p_sqrtRT", mesh->sqrtRT);
  kernelInfo.addDefine("p_sqrt2", (float)sqrt(2.));
  kernelInfo.addDefine("p_invsqrt2", (float)sqrt(1./2.));
  kernelInfo.addDefine("p_isq12", (float)sqrt(1./12.));
  kernelInfo.addDefine("p_isq6", (float)sqrt(1./6.));
  kernelInfo.addDefine("p_tauInv", mesh->tauInv);


  kernelInfo.addDefine("p_q1bar", q1bar);
  kernelInfo.addDefine("p_q2bar", q2bar);
  kernelInfo.addDefine("p_q3bar", q3bar);
  kernelInfo.addDefine("p_q4bar", q4bar);
  kernelInfo.addDefine("p_q5bar", q5bar);
  kernelInfo.addDefine("p_q6bar", q6bar);
  kernelInfo.addDefine("p_alpha0", (float).01f);


  mesh->volumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSplitPmlVolumeQuad2D.okl",
				       "boltzmannSplitPmlVolumeQuad2D",
				       kernelInfo);
  printf("starting surface\n");
  mesh->surfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSplitPmlSurfaceQuad2D.okl",
				       "boltzmannSplitPmlSurfaceQuad2D",
				       kernelInfo);
  printf("ending surface\n");

#if 0
  mesh->relaxationKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSplitPmlRelaxationQuad2D.okl",
				       "boltzmannSplitPmlRelaxationQuad2D",
				       kernelInfo);

  mesh->boltzmannPartialSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannPartialSurfaceQuad2D.okl",
				       "boltzmannPartialSurfaceQuad2D",
				       kernelInfo);
#endif

  mesh->updateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/boltzmannSplitPmlUpdateQuad2D.okl",
				       "boltzmannSplitPmlUpdateQuad2D",
				       kernelInfo);

  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
				       "meshHaloExtract2D",
				       kernelInfo);
  
}
