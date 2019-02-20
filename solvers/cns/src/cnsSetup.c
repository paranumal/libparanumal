/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "cns.h"

cns_t *cnsSetup(mesh_t *mesh, setupAide &options){

  cns_t *cns = (cns_t*) calloc(1, sizeof(cns_t));

  options.getArgs("MESH DIMENSION", cns->dim);
  options.getArgs("ELEMENT TYPE", cns->elementType);

  mesh->Nfields = (cns->dim==3) ? 4:3;
  cns->Nfields = mesh->Nfields;

  cns->Nstresses = (cns->dim==3) ? 6:3;
  cns->mesh = mesh;

  dlong Ntotal = mesh->Nelements*mesh->Np*mesh->Nfields;
  cns->Nblock = (Ntotal+blockSize-1)/blockSize;

  hlong localElements = (hlong) mesh->Nelements;
  MPI_Allreduce(&localElements, &(cns->totalElements), 1, MPI_HLONG, MPI_SUM, mesh->comm);

  // mean flow
  cns->rbar = 1;
  cns->ubar = 0.2;
  cns->vbar = 0;
  cns->wbar = 0;

  // viscosity
  dfloat Re = 5000;
  dfloat mu = 1;

  int check;

  check = options.getArgs("RBAR", cns->rbar);
  if(!check) printf("WARNING setup file does not include RBAR\n");

  check = options.getArgs("UBAR", cns->ubar);
  if(!check) printf("WARNING setup file does not include UBAR\n");

  check = options.getArgs("VBAR", cns->vbar);
  if(!check) printf("WARNING setup file does not include VBAR\n");

  check = options.getArgs("WBAR", cns->wbar);
  if(!check) printf("WARNING setup file does not include WBAR\n");

  check = options.getArgs("VISCOSITY", cns->mu);
  if(!check) printf("WARNING setup file does not include VISCOSITY\n");

  dfloat soundSpeed = 5;
  check = options.getArgs("SPEED OF SOUND", soundSpeed);
  if(!check) printf("WARNING setup file does not include MACH\n");

  // speed of sound (assuming isothermal unit bulk flow) = sqrt(RT)
  cns->RT = soundSpeed*soundSpeed;

  cns->outputForceStep = 0;

  options.getArgs("TSTEPS FOR FORCE OUTPUT",   cns->outputForceStep);

  // compute samples of q at interpolation nodes
  //  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
  //                                sizeof(dfloat));

  cns->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
			       sizeof(dfloat));

  cns->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
                                sizeof(dfloat));

  if (options.compareArgs("TIME INTEGRATOR","LSERK4")){
    cns->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
                                sizeof(dfloat));
  }

  if (options.compareArgs("TIME INTEGRATOR","DOPRI5")){
    int NrkStages = 7;
    cns->rkq  = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
          sizeof(dfloat));
    cns->rkrhsq = (dfloat*) calloc(NrkStages*mesh->Nelements*mesh->Np*mesh->Nfields,
          sizeof(dfloat));
    cns->rkerr  = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
          sizeof(dfloat));

    cns->errtmp = (dfloat*) calloc(cns->Nblock, sizeof(dfloat));

    // Dormand Prince -order (4) 5 with PID timestep control
    int Nrk = 7;
    dfloat rkC[7] = {0.0, 0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0};
    dfloat rkA[7*7] ={             0.0,             0.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                                   0.2,             0.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                              3.0/40.0,        9.0/40.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                             44.0/45.0,      -56.0/15.0,       32.0/9.0,          0.0,             0.0,       0.0, 0.0,
                        19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0,             0.0,       0.0, 0.0,
                         9017.0/3168.0,     -355.0/33.0, 46732.0/5247.0,   49.0/176.0, -5103.0/18656.0,       0.0, 0.0,
                            35.0/384.0,             0.0,   500.0/1113.0,  125.0/192.0,  -2187.0/6784.0, 11.0/84.0, 0.0 };
    dfloat rkE[7] = {71.0/57600.0,  0.0, -71.0/16695.0, 71.0/1920.0, -17253.0/339200.0, 22.0/525.0, -1.0/40.0 };

    cns->Nrk = Nrk;
    cns->rkC = (dfloat*) calloc(cns->Nrk, sizeof(dfloat));
    cns->rkE = (dfloat*) calloc(cns->Nrk, sizeof(dfloat));
    cns->rkA = (dfloat*) calloc(cns->Nrk*cns->Nrk, sizeof(dfloat));

    cns->rkoutB = (dfloat*) calloc(cns->Nrk, sizeof(dfloat));

    memcpy(cns->rkC, rkC, cns->Nrk*sizeof(dfloat));
    memcpy(cns->rkE, rkE, cns->Nrk*sizeof(dfloat));
    memcpy(cns->rkA, rkA, cns->Nrk*cns->Nrk*sizeof(dfloat));

    cns->ATOL    = 1.0; options.getArgs("ABSOLUTE TOLERANCE",   cns->ATOL);
    cns->RTOL    = 1.0; options.getArgs("RELATIVE TOLERANCE",   cns->RTOL);
    cns->dtMIN   = 1.0; options.getArgs("MINUMUM TIME STEP SIZE",   cns->dtMIN);
    cns->dtMAX   = 1.0; options.getArgs("MAXIMUM TIME STEP SIZE",   cns->dtMAX);

    cns->safe = 0.9;   //safety factor

    //error control parameters
    cns->beta = 0.05;
    cns->factor1 = 0.2;
    cns->factor2 = 10.0;


    cns->exp1 = 0.2 - 0.75*cns->beta;
    cns->invfactor1 = 1.0/cns->factor1;
    cns->invfactor2 = 1.0/cns->factor2;
    cns->facold = 1E-4;

  }

  cns->viscousStresses = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*cns->Nstresses,
                                           sizeof(dfloat));

  cns->Vort = (dfloat*) calloc(3*mesh->Nelements*mesh->Np,sizeof(dfloat)); // 3 components (hard coded)

  dfloat fx, fy, fz, intfx, intfy, intfz;
  cnsBodyForce(0.0, &fx, &fy, &fz, &intfx, &intfy, &intfz);

  printf("setting up initial condition\n");

  if(options.compareArgs("INITIAL CONDITION", "BROWN-MINION")){
    cnsBrownMinionQuad3D(cns);
  }
  else{

    // fix this later (initial conditions)
    for(dlong e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
	dfloat t = 0;
	dfloat x = mesh->x[n + mesh->Np*e];
	dfloat y = mesh->y[n + mesh->Np*e];
	dfloat z = mesh->z[n + mesh->Np*e];

	dlong qbase = e*mesh->Np*mesh->Nfields + n;

#if 0
	cnsGaussianPulse(x, y, z, t,
                       cns->q+qbase,
			 cns->q+qbase+mesh->Np,
			 cns->q+qbase+2*mesh->Np,
			 cns->q+qbase+3*mesh->Np);
#else
	cns->q[qbase+0*mesh->Np] = cns->rbar;
	cns->q[qbase+1*mesh->Np] = cns->rbar*intfx;
	cns->q[qbase+2*mesh->Np] = cns->rbar*intfy;
	if(cns->dim==3)
	  cns->q[qbase+3*mesh->Np] = cns->rbar*intfz;
#endif
      }
    }
  }

  // set penalty parameter
  mesh->Lambda2 = 0.5;

  // set time step
  dfloat hmin = 1e9;
  for(dlong e=0;e<mesh->Nelements;++e){

    for(int f=0;f<mesh->Nfaces;++f){
      dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
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
  dfloat cfl = 0.5; // depends on the stability region size

  dfloat dtAdv  = hmin/((mesh->N+1.)*(mesh->N+1.)*sqrt(cns->RT));
  dfloat dtVisc = pow(hmin, 2)/(pow(mesh->N+1,4)*cns->mu);

  dfloat dt = cfl*mymin(dtAdv, dtVisc);
  dt = cfl*dtAdv;

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

  //
  options.getArgs("FINAL TIME", mesh->finalTime);

  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  if (options.compareArgs("TIME INTEGRATOR","LSERK4")){
    mesh->dt = mesh->finalTime/mesh->NtimeSteps;
  }

  if (mesh->rank ==0) printf("dtAdv = %lg (before cfl), dtVisc = %lg (before cfl), dt = %lg\n",
   dtAdv, dtVisc, dt);

  cns->frame = 0;
  // errorStep
  mesh->errorStep = 1000;

  if (mesh->rank ==0) printf("dt = %g\n", mesh->dt);

  // OCCA build stuff

  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  printf("occa setup\n");

  if(cns->dim==3){
    if(cns->elementType != QUADRILATERALS)
      meshOccaSetup3D(mesh, options, kernelInfo);
    else
      meshOccaSetupQuad3D(mesh, options, kernelInfo);
  }
  else
    meshOccaSetup2D(mesh, options, kernelInfo);

  printf("occa array setup\n");

  //add boundary data to kernel info
  string boundaryHeaderFileName;
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo["includes"] += (char*)boundaryHeaderFileName.c_str();

  cns->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), cns->q);

  cns->o_saveq =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), cns->q);


  cns->o_viscousStresses =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*cns->Nstresses*sizeof(dfloat),
                        cns->viscousStresses);

  cns->o_rhsq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), cns->rhsq);

  if (mesh->rank==0)
    cout << "TIME INTEGRATOR (" << options.getArgs("TIME INTEGRATOR") << ")" << endl;

  if (options.compareArgs("TIME INTEGRATOR","LSERK4")){
    cns->o_resq =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), cns->resq);
  }

  if (options.compareArgs("TIME INTEGRATOR","DOPRI5")){
    int NrkStages = 7;
    cns->o_rkq =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), cns->rkq);
    cns->o_rkrhsq =
      mesh->device.malloc(NrkStages*mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), cns->rkrhsq);
    cns->o_rkerr =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), cns->rkerr);

    cns->o_errtmp = mesh->device.malloc(cns->Nblock*sizeof(dfloat), cns->errtmp);

    cns->o_rkA = mesh->device.malloc(cns->Nrk*cns->Nrk*sizeof(dfloat), cns->rkA);
    cns->o_rkE = mesh->device.malloc(  cns->Nrk*sizeof(dfloat), cns->rkE);

    cns->o_rkoutB = mesh->device.malloc(cns->Nrk*sizeof(dfloat), cns->rkoutB);
  }


  cns->o_Vort = mesh->device.malloc(3*mesh->Np*mesh->Nelements*sizeof(dfloat), cns->Vort); // 3 components


  if(mesh->totalHaloPairs>0){
    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat));

    cns->o_haloStressesBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*cns->Nstresses*sizeof(dfloat));

    // MPI send buffer
    cns->haloBytes = mesh->totalHaloPairs*mesh->Np*cns->Nfields*sizeof(dfloat);
    cns->haloStressesBytes = mesh->totalHaloPairs*mesh->Np*cns->Nstresses*sizeof(dfloat);

    cns->o_haloBuffer = mesh->device.malloc(cns->haloBytes);
    cns->o_haloStressesBuffer = mesh->device.malloc(cns->haloStressesBytes);

    cns->sendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, cns->haloBytes, NULL, cns->o_sendBuffer, cns->h_sendBuffer);
    cns->recvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, cns->haloBytes, NULL, cns->o_recvBuffer, cns->h_recvBuffer);
    cns->sendStressesBuffer = (dfloat*) occaHostMallocPinned(mesh->device, cns->haloStressesBytes, NULL, cns->o_sendStressesBuffer, cns->h_sendStressesBuffer);
    cns->recvStressesBuffer = (dfloat*) occaHostMallocPinned(mesh->device, cns->haloStressesBytes, NULL, cns->o_recvStressesBuffer, cns->h_recvStressesBuffer);
  }

  kernelInfo["defines/" "p_Nfields"]= mesh->Nfields;
  kernelInfo["defines/" "p_Nstresses"]= cns->Nstresses;

  kernelInfo["defines/" "p_RT"]= cns->RT;

  dfloat sqrtRT = sqrt(cns->RT);
  kernelInfo["defines/" "p_sqrtRT"]= sqrtRT;

  kernelInfo["defines/" "p_rbar"]= cns->rbar;
  kernelInfo["defines/" "p_ubar"]= cns->ubar;
  kernelInfo["defines/" "p_vbar"]= cns->vbar;
  kernelInfo["defines/" "p_wbar"]= cns->wbar;


  const dfloat p_one = 1.0, p_two = 2.0, p_half = 1./2., p_third = 1./3., p_zero = 0;

  kernelInfo["defines/" "p_two"]= p_two;
  kernelInfo["defines/" "p_one"]= p_one;
  kernelInfo["defines/" "p_half"]= p_half;
  kernelInfo["defines/" "p_third"]= p_third;
  kernelInfo["defines/" "p_zero"]= p_zero;

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = 512/mesh->Np; // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  int cubMaxNodes = mymax(mesh->Np, (mesh->intNfp*mesh->Nfaces));
  kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
  int cubMaxNodes1 = mymax(mesh->Np, (mesh->intNfp));
  kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

  kernelInfo["defines/" "p_Lambda2"]= 0.5f;
  if(cns->elementType==QUADRILATERALS && mesh->dim==3){
    kernelInfo["defines/" "p_fainv"] = (dfloat) 0.0;
    kernelInfo["defines/" "p_invRadiusSq"] = (dfloat) 1./(mesh->sphereRadius*mesh->sphereRadius);
  }

  kernelInfo["defines/" "p_blockSize"]= blockSize;


  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  // set kernel name suffix
  char *suffix;

  if(cns->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(cns->elementType==QUADRILATERALS){
    if(cns->dim==2)
      suffix = strdup("Quad2D");
    else
      suffix = strdup("Quad3D");
  }
  if(cns->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(cns->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  printf("Building kernels\n");

  for (int r=0;r<mesh->size;r++) {
    if (r==mesh->rank) {

      // kernels from volume file
      sprintf(fileName, DCNS "/okl/cnsVolume%s.okl", suffix);
      sprintf(kernelName, "cnsVolume%s", suffix);

      cns->volumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "cnsStressesVolume%s", suffix);
      cns->stressesVolumeKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // kernels from surface file
      sprintf(fileName, DCNS "/okl/cnsSurface%s.okl", suffix);
      sprintf(kernelName, "cnsSurface%s", suffix);

      cns->surfaceKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "cnsStressesSurface%s", suffix);
      cns->stressesSurfaceKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      if(cns->elementType != HEXAHEDRA){ //remove later
	// kernels from cubature volume file
	sprintf(fileName, DCNS "/okl/cnsCubatureVolume%s.okl", suffix);
	sprintf(kernelName, "cnsCubatureVolume%s", suffix);

	cns->cubatureVolumeKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);

	// kernels from cubature surface file
	sprintf(fileName, DCNS "/okl/cnsCubatureSurface%s.okl", suffix);
	sprintf(kernelName, "cnsCubatureSurface%s", suffix);

	cns->cubatureSurfaceKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }

      // kernels from vorticity file
      sprintf(fileName, DCNS "/okl/cnsVorticity%s.okl", suffix);
      sprintf(kernelName, "cnsVorticity%s", suffix);

      cns->vorticityKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);


      // kernels from update file
      cns->updateKernel =
        mesh->device.buildKernel(DCNS "/okl/cnsUpdate.okl",
                                           "cnsUpdate",
                                           kernelInfo);

      cns->rkUpdateKernel =
        mesh->device.buildKernel(DCNS "/okl/cnsUpdate.okl",
                                           "cnsRkUpdate",
                                           kernelInfo);
      cns->rkStageKernel =
        mesh->device.buildKernel(DCNS "/okl/cnsUpdate.okl",
                                           "cnsRkStage",
                                           kernelInfo);

      cns->rkOutputKernel =
        mesh->device.buildKernel(DCNS "/okl/cnsUpdate.okl",
                                           "cnsRkOutput",
                                           kernelInfo);

      cns->rkErrorEstimateKernel =
        mesh->device.buildKernel(DCNS "/okl/cnsUpdate.okl",
                                           "cnsErrorEstimate",
                                           kernelInfo);

      // fix this later
      mesh->haloExtractKernel =
        mesh->device.buildKernel(DHOLMES "/okl/meshHaloExtract3D.okl",
                                           "meshHaloExtract3D",
				 kernelInfo);

      if(cns->elementType==QUADRILATERALS && mesh->dim==3){
	sprintf(kernelName, "cnsConstrain%s", suffix);
	sprintf(fileName, DCNS "/okl/cnsConstrain%s.okl", suffix);
	cns->constrainKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }
    }
    MPI_Barrier(mesh->comm);
  }

  printf("done building kernels\n");

  return cns;
}
