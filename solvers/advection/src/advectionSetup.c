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

#include "advection.h"

void interpolateHex3D(dfloat *I, dfloat *x, int N, dfloat *Ix, int M);

advection_t *advectionSetup(mesh_t *mesh, setupAide &newOptions, char* boundaryHeaderFileName){

  advection_t *advection = (advection_t*) calloc(1, sizeof(advection_t));

  newOptions.getArgs("MESH DIMENSION", advection->dim);
  newOptions.getArgs("ELEMENT TYPE", advection->elementType);

  mesh->Nfields = 1;
  advection->Nfields = mesh->Nfields;

  advection->mesh = mesh;

  dlong Ntotal = mesh->Nelements*mesh->Np*mesh->Nfields;
  advection->Nblock = (Ntotal+blockSize-1)/blockSize;

  hlong localElements = (hlong) mesh->Nelements;
  MPI_Allreduce(&localElements, &(advection->totalElements), 1, MPI_HLONG, MPI_SUM, mesh->comm);

  // viscosity
  int check;

  // compute samples of q at interpolation nodes
  advection->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				     sizeof(dfloat));
  advection->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  if (newOptions.compareArgs("TIME INTEGRATOR","LSERK4")){
    advection->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
		  		sizeof(dfloat));
  }

  if (newOptions.compareArgs("TIME INTEGRATOR","DOPRI5")){
    int NrkStages = 7;
    advection->rkq  = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
          sizeof(dfloat));
    advection->rkrhsq = (dfloat*) calloc(NrkStages*mesh->Nelements*mesh->Np*mesh->Nfields,
          sizeof(dfloat));
    advection->rkerr  = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
          sizeof(dfloat));

    advection->errtmp = (dfloat*) calloc(advection->Nblock, sizeof(dfloat));

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

    advection->Nrk = Nrk;
    advection->rkC = (dfloat*) calloc(advection->Nrk, sizeof(dfloat));
    advection->rkE = (dfloat*) calloc(advection->Nrk, sizeof(dfloat));
    advection->rkA = (dfloat*) calloc(advection->Nrk*advection->Nrk, sizeof(dfloat));

    memcpy(advection->rkC, rkC, advection->Nrk*sizeof(dfloat));
    memcpy(advection->rkE, rkE, advection->Nrk*sizeof(dfloat));
    memcpy(advection->rkA, rkA, advection->Nrk*advection->Nrk*sizeof(dfloat));

    advection->dtMIN = 1E-9; //minumum allowed timestep
    advection->ATOL = 1E-6;  //absolute error tolerance
    advection->RTOL = 1E-6;  //relative error tolerance
    advection->safe = 0.8;   //safety factor

    //error control parameters
    advection->beta = 0.05;
    advection->factor1 = 0.2;
    advection->factor2 = 10.0;


    advection->exp1 = 0.2 - 0.75*advection->beta;
    advection->invfactor1 = 1.0/advection->factor1;
    advection->invfactor2 = 1.0/advection->factor2;
    advection->facold = 1E-4;

  }

  // fix this later (initial conditions)
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat t = 0;
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];
      dfloat z = mesh->z[n + mesh->Np*e];

      dlong qbase = e*mesh->Np*mesh->Nfields + n;

      dfloat qn = 0;

      advectionGaussianPulse(x, y, z, t, &qn);
      advection->q[qbase+0*mesh->Np] = qn;
    }
  }

  // set time step
  dfloat hmin = 1e9;

  if(advection->elementType == TRIANGLES || advection->elementType == TETRAHEDRA){

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
  }
  else{
    for(dlong e=0;e<mesh->Nelements;++e){
      for(int f=0;f<mesh->Nfaces;++f){
	for(int n=0;n<mesh->Nfp;++n){
	  dlong sid = mesh->Nsgeo*(mesh->Nfaces*mesh->Nfp*e + mesh->Nfp*f + n);
	  dfloat sJ   = mesh->sgeo[sid + SJID];
	  dfloat invJ = mesh->sgeo[sid + IJID];

	  if(invJ<0) printf("invJ = %g\n", invJ);
	  dfloat hest = .5/(sJ*invJ);

	  hmin = mymin(hmin, hest);
	}
      }
    }
  }

  // need to change cfl and defn of dt
  dfloat cfl = 0.5; // depends on the stability region size

  dfloat dtAdv  = hmin/((mesh->N+1.)*(mesh->N+1.));
  dfloat dt = cfl*dtAdv;

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

  //
  if(newOptions.getArgs("TIME STEPS", mesh->NtimeSteps)){
    mesh->finalTime = mesh->NtimeSteps*mesh->dt;
  }
  else{
    newOptions.getArgs("FINAL TIME", mesh->finalTime);

    mesh->NtimeSteps = mesh->finalTime/mesh->dt;

    if (newOptions.compareArgs("TIME INTEGRATOR","LSERK4")){
      mesh->dt = mesh->finalTime/mesh->NtimeSteps;
    }
  }

  if (mesh->rank ==0) printf("dtAdv = %lg (before cfl), dt = %lg\n",
   dtAdv, dt);

  advection->frame = 0;
  // errorStep
  mesh->errorStep = 1000;
  if (newOptions.compareArgs("TIME INTEGRATOR","LSERK4")){
    newOptions.getArgs("ERROR STEP", mesh->errorStep);
  }

  if (mesh->rank ==0) printf("dt = %g\n", mesh->dt);

  // OCCA build stuff

  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["compiler_flags"] = " ";

  //  if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
  {
    kernelInfo["compiler_flags"] += " --ftz=true";
    kernelInfo["compiler_flags"] += " --prec-div=false";
    kernelInfo["compiler_flags"] += " --prec-sqrt=false";
    kernelInfo["compiler_flags"] += " --use_fast_math";
    kernelInfo["compiler_flags"] += " --fmad=true"; // compiler option for cuda
    kernelInfo["compiler_flags"] += " -Xptxas -dlcm=ca";
  }

  if(advection->dim==3)
    meshOccaSetup3D(mesh, newOptions, kernelInfo);
  else
    meshOccaSetup2D(mesh, newOptions, kernelInfo);

  //add boundary data to kernel info
  kernelInfo["includes"] += boundaryHeaderFileName;

  advection->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), advection->q);

  advection->o_qtmp0 =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), advection->q);

  advection->o_qtmp1 =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), advection->q);

  advection->o_qtmp2 =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), advection->q);

  advection->o_saveq =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), advection->q);

  advection->o_rhsq =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), advection->rhsq);

  int advectionForm = 0, advectionIntegration = 0, advectionMassType = 0, advectionCombined = 0;
  if (newOptions.compareArgs("ADVECTION FORMULATION","WEAK")){
    advectionForm = 0; // DEFAULT
  }
  if (newOptions.compareArgs("ADVECTION FORMULATION","SKEW")){
    advectionForm = 1;
  }

  if (newOptions.compareArgs("ADVECTION FORMULATION","NODAL")){
    advectionIntegration = 0; // DEFAULT
  }

  if (newOptions.compareArgs("ADVECTION FORMULATION","CUBATURE")){
    advectionIntegration = 1;
  }


  if (newOptions.compareArgs("ADVECTION FORMULATION","SEMDG")){
    advectionMassType = 0;  // DEFAULT
  }

  if (newOptions.compareArgs("ADVECTION FORMULATION","WADG")){
    advectionMassType = 1;
  }

  if (newOptions.compareArgs("ADVECTION FORMULATION","MASS")){
    advectionMassType = 2;
  }

  
  if (newOptions.compareArgs("ADVECTION FORMULATION","COMBINED")){
    advectionCombined = 1;
  }



  // non-constant advection velocity
  dfloat t = 0;
  hlong totalNelements = mesh->Nelements+mesh->totalHaloPairs;
  dfloat *advectionVelocityJW = (dfloat*) calloc(mesh->Np*totalNelements*mesh->dim,sizeof(dfloat));
  dfloat *advectionVelocityM  = (dfloat*) calloc(mesh->Nfaces*mesh->Nfp*totalNelements,sizeof(dfloat));
  dfloat *advectionVelocityP  = (dfloat*) calloc(mesh->Nfaces*mesh->Nfp*totalNelements,sizeof(dfloat));

  dfloat *cubAdvectionVelocityJW = (dfloat*) calloc(mesh->cubNp*totalNelements*mesh->dim,sizeof(dfloat));

  for(dlong e=0;e<mesh->Nelements;++e){

    for(int n=0;n<mesh->Np;++n){
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];
      dfloat z = mesh->z[n + mesh->Np*e];

      dlong gbase = e*mesh->Np*mesh->Nvgeo + n;

      dfloat JW = mesh->vgeo[gbase + mesh->Np*JWID];

      dfloat rx = mesh->vgeo[gbase+mesh->Np*RXID];
      dfloat sx = mesh->vgeo[gbase+mesh->Np*SXID];
      dfloat tx = mesh->vgeo[gbase+mesh->Np*TXID];
      dfloat ry = mesh->vgeo[gbase+mesh->Np*RYID];
      dfloat sy = mesh->vgeo[gbase+mesh->Np*SYID];
      dfloat ty = mesh->vgeo[gbase+mesh->Np*TYID];
      dfloat rz = mesh->vgeo[gbase+mesh->Np*RZID];
      dfloat sz = mesh->vgeo[gbase+mesh->Np*SZID];
      dfloat tz = mesh->vgeo[gbase+mesh->Np*TZID];

      dlong qbase = e*mesh->Np*mesh->dim + n;

      dfloat cx = -y;
      dfloat cy = +x;
      dfloat cz = 0;

      advectionVelocityJW[qbase + 0*mesh->Np] = JW*(rx*cx+ry*cy+rz*cz);
      advectionVelocityJW[qbase + 1*mesh->Np] = JW*(sx*cx+sy*cy+sz*cz);
      if(mesh->dim==3)
	advectionVelocityJW[qbase + 2*mesh->Np] = JW*(tx*cx+ty*cy+tz*cz);
    }


    for(int f=0;f<mesh->Nfaces;++f){
      for(int n=0;n<mesh->Nfp;++n){
	int id = e*mesh->Nfaces*mesh->Nfp + f*mesh->Nfp + n;
	int vid = mesh->vmapM[id];
	dfloat x = mesh->x[vid];
	dfloat y = mesh->y[vid];
	dfloat z = mesh->z[vid];

	dlong sbase = mesh->Nsgeo*(e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n);

	dfloat nx = mesh->sgeo[sbase + NXID];
	dfloat ny = mesh->sgeo[sbase + NYID];
	dfloat nz = mesh->sgeo[sbase + NZID];

	dlong gbase = e*mesh->Np*mesh->Nvgeo + (vid%mesh->Np);

	dfloat JW = mesh->vgeo[gbase + mesh->Np*JWID];

	dfloat cx = -y;
	dfloat cy = +x;
	dfloat cz = 0;

	dfloat cdotn = nx*cx + ny*cy + nz*cz;
	dfloat LIFT = (advectionMassType==2) ? mesh->sgeo[sbase+WSJID]:mesh->sgeo[sbase + WIJID]*mesh->sgeo[sbase+SJID];
	dfloat alpha = 1; // not allowed to use central here

	// fix later
	if(advectionForm==0){ // weak
	  advectionVelocityM[id] = LIFT*0.5*(cdotn - alpha*fabs(cdotn));
	  advectionVelocityP[id] = LIFT*0.5*(cdotn + alpha*fabs(cdotn));
	}
	if(advectionForm==1){ // skew
	  advectionVelocityM[id] = LIFT*0.5*(      - alpha*fabs(cdotn));
	  advectionVelocityP[id] = LIFT*0.5*(cdotn + alpha*fabs(cdotn));
	}

      }
    }
  }

  dfloat *cubx = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cuby = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cubz = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));

  for(dlong e=0;e<mesh->Nelements;++e){

    interpolateHex3D(mesh->cubInterp, mesh->x+e*mesh->Np, mesh->Nq, cubx, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, mesh->y+e*mesh->Np, mesh->Nq, cuby, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, mesh->z+e*mesh->Np, mesh->Nq, cubz, mesh->cubNq);

    for(int n=0;n<mesh->cubNp;++n){
      dfloat x = cubx[n];
      dfloat y = cuby[n];
      dfloat z = cubz[n];

      dlong gbase = e*mesh->cubNp*mesh->Nvgeo + n;

      dfloat JW = mesh->cubvgeo[gbase + mesh->cubNp*JWID];

      dfloat rx = mesh->cubvgeo[gbase+mesh->cubNp*RXID];
      dfloat sx = mesh->cubvgeo[gbase+mesh->cubNp*SXID];
      dfloat tx = mesh->cubvgeo[gbase+mesh->cubNp*TXID];
      dfloat ry = mesh->cubvgeo[gbase+mesh->cubNp*RYID];
      dfloat sy = mesh->cubvgeo[gbase+mesh->cubNp*SYID];
      dfloat ty = mesh->cubvgeo[gbase+mesh->cubNp*TYID];
      dfloat rz = mesh->cubvgeo[gbase+mesh->cubNp*RZID];
      dfloat sz = mesh->cubvgeo[gbase+mesh->cubNp*SZID];
      dfloat tz = mesh->cubvgeo[gbase+mesh->cubNp*TZID];

      dlong qbase = e*mesh->cubNp*mesh->dim + n;

      dfloat cx = -y;
      dfloat cy = +x;
      dfloat cz = 0;

      cubAdvectionVelocityJW[qbase + 0*mesh->cubNp] = JW*(rx*cx+ry*cy+rz*cz);
      cubAdvectionVelocityJW[qbase + 1*mesh->cubNp] = JW*(sx*cx+sy*cy+sz*cz);
      if(mesh->dim==3)
	cubAdvectionVelocityJW[qbase + 2*mesh->cubNp] = JW*(tx*cx+ty*cy+tz*cz);
    }
  }


  advection->o_advectionVelocityJW =
    mesh->device.malloc(mesh->Np*totalNelements*mesh->dim*sizeof(dfloat), advectionVelocityJW);

  advection->o_advectionVelocityM =
    mesh->device.malloc(mesh->Nfaces*mesh->Nfp*totalNelements*sizeof(dfloat), advectionVelocityM);

  advection->o_advectionVelocityP =
    mesh->device.malloc(mesh->Nfaces*mesh->Nfp*totalNelements*sizeof(dfloat), advectionVelocityP);

  advection->o_cubAdvectionVelocityJW =
    mesh->device.malloc(mesh->cubNp*totalNelements*mesh->dim*sizeof(dfloat), cubAdvectionVelocityJW);

  cout << "TIME INTEGRATOR (" << newOptions.getArgs("TIME INTEGRATOR") << ")" << endl;

  if (newOptions.compareArgs("TIME INTEGRATOR","LSERK4")){
    advection->o_resq =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), advection->resq);
  }

  if (newOptions.compareArgs("TIME INTEGRATOR","DOPRI5")){
    printf("setting up DOPRI5\n");
    int NrkStages = 7;
    advection->o_rkq =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), advection->rkq);
    advection->o_rkrhsq =
      mesh->device.malloc(NrkStages*mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), advection->rkrhsq);
    advection->o_rkerr =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), advection->rkerr);

    advection->o_errtmp = mesh->device.malloc(advection->Nblock*sizeof(dfloat), advection->errtmp);

    advection->o_rkA = mesh->device.malloc(advection->Nrk*advection->Nrk*sizeof(dfloat), advection->rkA);
    advection->o_rkE = mesh->device.malloc(  advection->Nrk*sizeof(dfloat), advection->rkE);
  }


  if(mesh->totalHaloPairs>0){
    // NOTE USE OF NFP NODES PER HALO FACE

    // MPI send buffer
    advection->haloBytes = mesh->totalHaloPairs*mesh->Nfp*advection->Nfields*sizeof(dfloat);
    
    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(advection->haloBytes*sizeof(dfloat));

    advection->o_haloBuffer =
      mesh->device.malloc(advection->haloBytes);

    advection->sendBuffer =
      (dfloat*) occaHostMallocPinned(mesh->device, advection->haloBytes, NULL, advection->o_sendBuffer);

    advection->recvBuffer =
      (dfloat*) occaHostMallocPinned(mesh->device, advection->haloBytes, NULL, advection->o_recvBuffer);
  }

  //  p_RT, p_rbar, p_ubar, p_vbar
  // p_half, p_two, p_third, p_Nstresses

  kernelInfo["defines/" "p_Nfields"]= mesh->Nfields;
  const dfloat p_one = 1.0, p_two = 2.0, p_half = 1./2., p_third = 1./3., p_zero = 0;

  kernelInfo["defines/" "p_two"]= p_two;
  kernelInfo["defines/" "p_one"]= p_one;
  kernelInfo["defines/" "p_half"]= p_half;
  kernelInfo["defines/" "p_third"]= p_third;
  kernelInfo["defines/" "p_zero"]= p_zero;

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = mymax(1,1024/mesh->Np); // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = mymax(1,1024/maxNodes); // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  int NblockAdvVol = 1; // works for CUDA
  kernelInfo["defines/" "p_NblockAdvVol"]= NblockAdvVol;

  int NblockAdvSur = 1; // works for CUDA
  kernelInfo["defines/" "p_NblockAdvSur"]= NblockAdvSur;

  int cubMaxNodes = mymax(mesh->Np, (mesh->intNfp*mesh->Nfaces));
  kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
  int cubMaxNodes1 = mymax(mesh->Np, (mesh->intNfp));
  kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

  kernelInfo["defines/" "p_Lambda2"]= 0.5f;

  kernelInfo["defines/" "p_blockSize"]= blockSize;


  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  //  std::cout << kernelInfo << std::endl;

  // set kernel name suffix
  char *suffix;

  if(advection->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(advection->elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(advection->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(advection->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  // kernels from volume file
  sprintf(fileName, DADVECTION "/okl/advectionVolume%s.okl", suffix);

  string kName = "advection";

  if(advectionCombined) kName += "Combined";

  kName += (advectionIntegration==0) ? "Nodal":"Cubature";
  kName += (advectionForm==0)        ? "Weak":"Skew";
  if(advectionMassType==0)
    kName += "SEMDG";
  if(advectionMassType==1)
    kName += "WADG";
  if(advectionMassType==2)
    kName += "SEMDG";
  kName += "Volume";

  sprintf(kernelName, "%s%s", kName.c_str(), suffix);

  printf("kernelName = %s\n", kernelName);

  advection->volumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

#if 0
  if(advectionForm==0) // weak
    sprintf(kernelName, "advectionCombined%s", suffix);

  if(advectionForm==1) // skew
    sprintf(kernelName, "advectionCombinedSkew%s", suffix);

  advection->combinedKernel =
    mesh->device.buildKernel(fileName, kernelName, kernelInfo);
#endif

  // kernels from surface file
  sprintf(fileName, DADVECTION "/okl/advectionSurface%s.okl", suffix);
  sprintf(kernelName, "advectionSurface%s", suffix);

  advection->surfaceKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);

  // kernels from update file
  advection->updateKernel =
    mesh->device.buildKernel(DADVECTION "/okl/advectionUpdate.okl",
				       "advectionUpdate",
				       kernelInfo);

  advection->rkUpdateKernel =
    mesh->device.buildKernel(DADVECTION "/okl/advectionUpdate.okl",
				       "advectionRkUpdate",
				       kernelInfo);
  advection->rkStageKernel =
    mesh->device.buildKernel(DADVECTION "/okl/advectionUpdate.okl",
				       "advectionRkStage",
				       kernelInfo);

  advection->rkErrorEstimateKernel =
    mesh->device.buildKernel(DADVECTION "/okl/advectionUpdate.okl",
				       "advectionErrorEstimate",
				       kernelInfo);

  // fix this later
  mesh->haloExtractKernel =
    mesh->device.buildKernel(DHOLMES "/okl/meshHaloExtract3D.okl",
				       "meshHaloExtract3D",
				       kernelInfo);
  mesh->haloGetKernel =
    mesh->device.buildKernel(DHOLMES "/okl/meshHaloGet.okl",
			     "meshHaloGet",
			     kernelInfo);
  mesh->haloPutKernel =
    mesh->device.buildKernel(DHOLMES "/okl/meshHaloPut.okl",
				       "meshHaloPut",
				       kernelInfo);

  sprintf(fileName, DADVECTION "/okl/advectionInvertMassMatrix%s.okl", suffix);
  sprintf(kernelName, "advectionInvertMassMatrix%s", suffix);

  advection->invertMassMatrixKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);

  sprintf(fileName, DADVECTION "/okl/advectionInvertMassMatrix%s.okl", suffix);
  sprintf(kernelName, "advectionCombinedNodalWeakMMDGVolume%s", suffix);

  advection->invertMassMatrixCombinedKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);
  
  
  return advection;
}
