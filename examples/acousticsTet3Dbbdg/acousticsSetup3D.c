#include "acoustics3D.h"

iint factorial(iint n) {
  iint retval = 1;
  for (iint i = n; i > 1; --i) retval *= i;
  return retval;
}


void acousticsSetup3D(mesh3D *mesh){

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // set time step
  mesh->finalTime = .5;
  dfloat cfl = 0.5; // depends on the stability region size

  // set penalty parameter
  mesh->Lambda2 = 0.5;

  dfloat *EtoDT = (dfloat *) calloc(mesh->Nelements,sizeof(dfloat));
  dfloat hmin = 1e9;
  for(iint e=0;e<mesh->Nelements;++e){
    EtoDT[e] = 1e9;

    for(iint f=0;f<mesh->Nfaces;++f){
      iint sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
      // h = 0.5/(sJ/J)

      dfloat hest = .5/(sJ*invJ);

      // dt ~ cfl (h/(N+1)^2)/(Lambda^2*fastest wave speed)
      dfloat dtEst = cfl*hest/((mesh->N+1.)*(mesh->N+1.)*mesh->Lambda2);

      hmin = mymin(hmin, hest);
      EtoDT[e] = mymin(EtoDT[e], dtEst);
    }
  }

  //use dt on each element to setup MRAB
  int maxLevels = 100;
  meshMRABSetup3D(mesh,EtoDT,maxLevels);


  mesh->Nfields = 4;

  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
        sizeof(dfloat));
  mesh->fQM   = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*mesh->Nfields,
        sizeof(dfloat));
  mesh->fQP   = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*mesh->Nfields,
        sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(3*mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  // fix this later (initial conditions)
  iint cnt = 0;
  dfloat time = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];
      dfloat z = mesh->z[n + mesh->Np*e];

      cnt = e*mesh->Np*mesh->Nfields + n*mesh->Nfields;
      //acousticsGaussianPulse3D(x, y, z, time,
			//	mesh->q+cnt,
			//	mesh->q+cnt+1,
			//	mesh->q+cnt+2,
			//	mesh->q+cnt+3);

      mesh->q[cnt+0] = 0.;
      mesh->q[cnt+1] = 0.;
      mesh->q[cnt+2] = 0.;
      mesh->q[cnt+3] = 0.;
    }
  }

  //Transform to BB modal space
  #if USE_BERN
    dfloat qtmp[mesh->Nfields*mesh->Np];
    for (iint e =0;e<mesh->Nelements;e++){
      cnt = e*mesh->Np*mesh->Nfields;

      for (iint n=0; n<mesh->Np; n++){
        qtmp[n*mesh->Nfields + 0] = mesh->q[cnt+n*mesh->Nfields+0];
        qtmp[n*mesh->Nfields + 1] = mesh->q[cnt+n*mesh->Nfields+1];
        qtmp[n*mesh->Nfields + 2] = mesh->q[cnt+n*mesh->Nfields+2];
        qtmp[n*mesh->Nfields + 3] = mesh->q[cnt+n*mesh->Nfields+3];
        mesh->q[cnt+n*mesh->Nfields+0] = 0.0;
        mesh->q[cnt+n*mesh->Nfields+1] = 0.0;
        mesh->q[cnt+n*mesh->Nfields+2] = 0.0;
        mesh->q[cnt+n*mesh->Nfields+3] = 0.0;
      }
      for (iint n=0;n<mesh->Np;n++){
        for (iint m=0; m<mesh->Np; m++){
          mesh->q[cnt+n*mesh->Nfields + 0] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+0];
          mesh->q[cnt+n*mesh->Nfields + 1] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+1];
          mesh->q[cnt+n*mesh->Nfields + 2] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+2];
          mesh->q[cnt+n*mesh->Nfields + 3] += mesh->invVB[n*mesh->Np+m]*qtmp[m*mesh->Nfields+3];
        }
      }
    }

    //Change cubature Interp and Project matrices
    for (iint n=0;n<mesh->cubNp;n++) {
      dfloat r = mesh->cubr[n];
      dfloat s = mesh->cubs[n];
      dfloat t = mesh->cubt[n];

      dfloat l0 = -0.5*(1.+r+s+t); dfloat l1 = 0.5*(1.+r); dfloat l2 = 0.5*(1.+s); dfloat l3 = 0.5*(1.+t);

      iint cnt = 0;
      for (iint i=0;i<=mesh->N;i++){
        for (iint j=0;j<=mesh->N-i;j++){
          for (iint k=0;k<=mesh->N-i-j;k++){
            mesh->cubInterp[n*mesh->Np+cnt] = ((dfloat) factorial(mesh->N)/(factorial(i)*factorial(j)
                                            *factorial(k)*factorial(mesh->N-i-j-k)))
                                            *pow(l0,mesh->N-i-j-k)*pow(l1,k)*pow(l2,j)*pow(l3,i);
            cnt++;
          }
        }
      }
    }

    dfloat S[mesh->Np*mesh->cubNp];
    for (iint n=0;n<mesh->Np;n++) {
      for (iint m =0;m<mesh->cubNp;m++) {
        S[n*mesh->cubNp + m] = mesh->cubProject[n*mesh->cubNp + m];
      }
    }
    for (iint n=0;n<mesh->Np;n++) {
      for (iint m =0;m<mesh->cubNp;m++) {
        mesh->cubProject[n*mesh->cubNp + m] = 0.;
        for (iint i =0;i<mesh->Np;i++)
          mesh->cubProject[n*mesh->cubNp+m]
            += mesh->invVB[n*mesh->Np + i]*S[i*mesh->cubNp+m];
      }
    }
  #endif

  #if WADG
    // set heterogeneous c^2 for WADG
    mesh->c2 = (dfloat*) calloc(mesh->Nelements*mesh->cubNp,sizeof(dfloat));

    for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

      iint id = e*mesh->Nverts+0;

      dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
      dfloat xe2 = mesh->EX[id+1];
      dfloat xe3 = mesh->EX[id+2];
      dfloat xe4 = mesh->EX[id+3];

      dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
      dfloat ye2 = mesh->EY[id+1];
      dfloat ye3 = mesh->EY[id+2];
      dfloat ye4 = mesh->EY[id+3];

      dfloat ze1 = mesh->EZ[id+0]; /* y-coordinates of vertices */
      dfloat ze2 = mesh->EZ[id+1];
      dfloat ze3 = mesh->EZ[id+2];
      dfloat ze4 = mesh->EZ[id+3];

      for(iint n=0;n<mesh->cubNp;++n){ /* for each node */

        // cubature node coordinates
        dfloat rn = mesh->cubr[n];
        dfloat sn = mesh->cubs[n];
        dfloat tn = mesh->cubt[n];

        /* physical coordinate of interpolation node */
        dfloat x = -0.5*(rn+sn+tn+1.)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3 + 0.5*(1+tn)*xe4 ;
        dfloat y = -0.5*(rn+sn+tn+1.)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3 + 0.5*(1+tn)*ye4 ;
        dfloat z = -0.5*(rn+sn+tn+1.)*ze1 + 0.5*(1+rn)*ze2 + 0.5*(1+sn)*ze3 + 0.5*(1+tn)*ze4 ;

        // smoothly varying (sinusoidal) wavespeed
        //printf("M_PI = %f\n",M_PI);
        if (z<0.f) {
          mesh->c2[n + mesh->cubNp*e] = 1.0;//1.0 + 0.5*sin(M_PI*y);
        } else {
          mesh->c2[n + mesh->cubNp*e] = 0.2;
        }
      }
    }

  #endif

  // errorStep
  mesh->errorStep = 100;

  if (rank==0) {
    printf("hmin = %g\n", hmin);
    printf("cfl = %g\n", cfl);
    printf("dt = %g\n", mesh->dt);
    printf("max wave speed = %g\n", sqrt(3.)*mesh->sqrtRT);
  }

  // output mesh
  meshVTU3D(mesh, "foo.vtu");

  // OCCA build stuff
  char deviceConfig[BUFSIZ];

  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", 0);
  //sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 1");
  //sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 0);
  //sprintf(deviceConfig, "mode = Serial");

  mesh->device.setup(deviceConfig);

  occa::kernelInfo kernelInfo;

  // generic occa device set up
  meshOccaSetup3D(mesh, deviceConfig, kernelInfo);

  //reallocate rhsq for MRAB
  mesh->o_rhsq.free();
  mesh->o_rhsq =
    mesh->device.malloc(3*mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);

  if (mesh->totalHaloPairs) {
    //reallocate halo buffer for trace exchange
    mesh->o_haloBuffer.free();
    mesh->o_haloBuffer =
        mesh->device.malloc(mesh->totalHaloPairs*mesh->Nfp*mesh->Nfields*mesh->Nfaces*sizeof(dfloat));
  }

  mesh->o_c2 = mesh->device.malloc(mesh->Nelements*mesh->cubNp*sizeof(dfloat),
           mesh->c2);

  mesh->o_fQM = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*mesh->Nfields*sizeof(dfloat),
         mesh->fQM);
  mesh->o_fQP = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfp*mesh->Nfaces*mesh->Nfields*sizeof(dfloat),
         mesh->fQP);
  mesh->o_mapP = mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
         mesh->mapP);

  //set up pml
  acousticsPmlSetup3D(mesh);

  //set up source injection
  acousticsSourceSetup3D(mesh,kernelInfo);

  mesh->o_MRABelementIds = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
  mesh->o_MRABhaloIds = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
  for (iint lev=0;lev<mesh->MRABNlevels;lev++) {
    if (mesh->MRABNelements[lev])
      mesh->o_MRABelementIds[lev] = mesh->device.malloc(mesh->MRABNelements[lev]*sizeof(iint),
         mesh->MRABelementIds[lev]);
    if (mesh->MRABNhaloElements[lev])
      mesh->o_MRABhaloIds[lev] = mesh->device.malloc(mesh->MRABNhaloElements[lev]*sizeof(iint),
         mesh->MRABhaloIds[lev]);
  }


  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  int maxCubNodes = mymax(maxNodes,mesh->cubNp);

  kernelInfo.addDefine("p_maxNodes", maxNodes);
  kernelInfo.addDefine("p_maxCubNodes", maxCubNodes);

  int NblockV = 512/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  int NblockCub = 512/mesh->cubNp; // works for CUDA
  kernelInfo.addDefine("p_NblockCub", NblockCub);


  kernelInfo.addDefine("p_Lambda2", 0.5f);

  mesh->volumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsbbdgMRABVolume3D.okl",
				       "acousticsbbdgMRABVolume3D",
				       kernelInfo);

  mesh->surfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsbbdgMRABSurface3D.okl",
				       "acousticsbbdgMRABSurface3D",
				       kernelInfo);

  #if WADG
      mesh->updateKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABUpdate3D.okl",
              "acousticsMRABUpdate3D_wadg",
              kernelInfo);
      mesh->traceUpdateKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABUpdate3D.okl",
              "acousticsMRABTraceUpdate3D_wadg",
               kernelInfo);
      mesh->pmlUpdateKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABPmlUpdate3D.okl",
               "acousticsMRABPmlUpdate3D_wadg",
                 kernelInfo);
      mesh->pmlTraceUpdateKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABPmlUpdate3D.okl",
               "acousticsMRABPmlTraceUpdate3D_wadg",
                 kernelInfo);
  #else
    mesh->updateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABUpdate3D.okl",
	 			       "acousticsMRABUpdate3D",
				       kernelInfo);
    mesh->traceUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABUpdate3D.okl",
               "acousticsMRABTraceUpdate3D",
               kernelInfo);
    mesh->pmlUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABPmlUpdate3D.okl",
               "acousticsMRABPmlUpdate3D",
                 kernelInfo);
    mesh->pmlTraceUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABPmlUpdate3D.okl",
               "acousticsMRABPmlTraceUpdate3D",
                 kernelInfo);
  #endif

  mesh->pmlVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsbbdgMRABPmlVolume3D.okl",
             "acousticsbbdgMRABPmlVolume3D",
             kernelInfo);
  mesh->pmlSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsbbdgMRABPmlSurface3D.okl",
             "acousticsbbdgMRABPmlSurface3D",
             kernelInfo);

  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract3D.okl",
				       "meshHaloExtract3D",
				       kernelInfo);

}
