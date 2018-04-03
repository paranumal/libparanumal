#include "mpi.h"
#include <math.h>
#include "cnsQuad2D.h"

cns_t *cnsSetupQuad2D(mesh2D *mesh){

  cns_t *cns = (cns_t*) calloc(1, sizeof(cns_t));

  mesh->Nfields = 4;

  cns->Nfields = mesh->Nfields;
  cns->Nstresses = 3;
  cns->mesh = mesh;
  
  // speed of sound (assuming isothermal unit bulk flow) = sqrt(RT)
  cns->RT = 10;

  // viscosity
  cns->mu = 1e-2;
  
  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  cns->viscousStresses = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*cns->Nstresses,
					   sizeof(dfloat));
  
  // fix this later (initial conditions)
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat t = 0;
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];

      dlong qbase = e*mesh->Np*mesh->Nfields + n;
      
      cnsGaussianPulse2D(x, y, t,
			 mesh->q+qbase,
			 mesh->q+qbase+mesh->Np,
			 mesh->q+qbase+2*mesh->Np);
    }
  }

  // set penalty parameter
  mesh->Lambda2 = 0.5;
  
  // set time step
  dfloat hmin = 1e9;
  for(int e=0;e<mesh->Nelements;++e){  

    for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      int sid = mesh->Nsgeo*(mesh->Nfaces*mesh->Nfp*e + n);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      if(invJ<0) printf("invJ = %g\n", invJ);
      
      // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
      // h = 0.5/(sJ/J)
      
      dfloat hest = .5/(sJ*invJ);

      hmin = mymin(hmin, hest);
    }
  }

  // need to change cfl and defn of dt
  dfloat cfl = .1; // depends on the stability region size

  dfloat dt = cfl*hmin/((mesh->N+1.)*(mesh->N+1.)*mesh->Lambda2);

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  //
  mesh->finalTime = 1;
  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  mesh->dt = mesh->finalTime/mesh->NtimeSteps;

  // errorStep
  mesh->errorStep = 100;

  printf("dt = %g\n", mesh->dt);

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", rank%3);
  mesh->device.setup(deviceConfig);

 
  cns->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);
  
  cns->o_viscousStresses =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*cns->Nstresses*sizeof(dfloat), mesh->q);
  
  cns->o_rhsq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);

  cns->o_resq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->resq);
  
  mesh->o_D  = mesh->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D);
  
  mesh->o_LIFT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
			mesh->LIFT);

  cns->LIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Nfp*mesh->Nfaces;++m){
      cns->LIFTT[n + m*mesh->Np] = mesh->LIFT[n*mesh->Nfaces*mesh->Nfp+m];
    }
  }
  
  cns->o_LIFTT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
			cns->LIFTT);
  
  mesh->o_vgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nvgeo*sizeof(dfloat),
			mesh->vgeo);
  
  mesh->o_sgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->Nsgeo*sizeof(dfloat),
			mesh->sgeo);

  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(int),
			mesh->vmapM);

  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(int),
			mesh->vmapP);

  if(mesh->totalHaloPairs>0){
    // copy halo element list to DEVICE
    mesh->o_haloElementList =
      mesh->device.malloc(mesh->totalHaloPairs*sizeof(int), mesh->haloElementList);
    
    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat));
  }

  occa::kernelInfo kernelInfo;

  // generic occa device set up
  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

  //  p_RT, p_rbar, p_ubar, p_vbar, p_RXID, p_SXID, p_RYID, p_SYID
  // p_half, p_two, p_third, p_Nstresses
  
  kernelInfo.addDefine("p_Nfields", mesh->Nfields);
  kernelInfo.addDefine("p_Nstresses", cns->Nstresses);

  kernelInfo.addDefine("p_RT", cns->RT);

  dfloat sqrtRT = sqrt(cns->RT);
  kernelInfo.addDefine("p_sqrtRT", sqrtRT);
  
  kernelInfo.addDefine("p_rbar", cns->rbar);
  kernelInfo.addDefine("p_ubar", cns->ubar);
  kernelInfo.addDefine("p_vbar", cns->vbar);
  
  kernelInfo.addDefine("p_RXID", RXID);
  kernelInfo.addDefine("p_RYID", RYID);
  kernelInfo.addDefine("p_SXID", SXID);
  kernelInfo.addDefine("p_SYID", SYID);
  kernelInfo.addDefine("p_JID", JID);

  const dfloat p_one = 1.0, p_two = 2.0, p_half = 1./2., p_third = 1./3.;

  kernelInfo.addDefine("p_two", p_two);
  kernelInfo.addDefine("p_one", p_one);
  kernelInfo.addDefine("p_half", p_half);
  kernelInfo.addDefine("p_third", p_third);
  
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 512/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  kernelInfo.addDefine("p_Lambda2", 0.5f);
  
  cns->volumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsVolumeQuad2D.okl",
				       "cnsVolumeQuad2D",
				       kernelInfo);

  cns->surfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsSurfaceQuad2D.okl",
				       "cnsSurfaceQuad2D",
				       kernelInfo);

  cns->stressesVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsVolumeQuad2D.okl",
				       "cnsStressesVolumeQuad2D",
				       kernelInfo);

  cns->stressesSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsSurfaceQuad2D.okl",
				       "cnsStressesSurfaceQuad2D",
				       kernelInfo);
  

  cns->updateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/cnsUpdateQuad2D.okl",
				       "cnsUpdateQuad2D",
				       kernelInfo);
  
  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
				       "meshHaloExtract2D",
				       kernelInfo);

  return cns;
}
