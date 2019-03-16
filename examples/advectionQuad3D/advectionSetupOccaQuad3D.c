#include "advectionQuad3D.h"

void advectionSetupOccaQuad3D(solver_t *solver,occa::kernelInfo *kernelInfo) {

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  mesh_t *mesh = solver->mesh;
  
  // use rank to choose DEVICE
  //sprintf(solver->deviceConfig, "mode = CUDA, deviceID = %d", (rank+1)%2);
  //sprintf(solver->deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 0");
  sprintf(solver->deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(solver->deviceConfig, "mode = Serial");
  
  solver->device.setup(solver->deviceConfig);
  
  occa::initTimer(solver->device);
  
  dfloat *LIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Nfaces*mesh->Nfp;++m){
      LIFTT[n+m*mesh->Np] = mesh->LIFT[n*mesh->Nfp*mesh->Nfaces+m];
    }
  }

  dfloat *MLIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Nfaces*mesh->Nfp;++m){
      MLIFTT[n+m*mesh->Np] = mesh->MLIFT[n*mesh->Nfp*mesh->Nfaces+m];
    }
  }

  solver->o_D = solver->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D);
  solver->o_weakD = solver->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->MD);

  solver->o_mass = solver->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->weakD);

  solver->o_LIFT =
    solver->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
			  mesh->MLIFT);
  
  solver->o_LIFTT =
    solver->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
			  MLIFTT);

  solver->o_vgeo =
    solver->device.malloc(mesh->Nelements*mesh->Nvgeo*mesh->Np*sizeof(dfloat),
			mesh->vgeo);
  solver->o_sgeo =
    solver->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->Nfp*mesh->Nsgeo*sizeof(dfloat),
			mesh->sgeo);
  
  solver->o_vmapM =
    solver->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
			  mesh->vmapM);
  
  solver->o_vmapP =
    solver->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
			  mesh->vmapP);
  
  solver->o_mapP =
    solver->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
			  mesh->mapP);

  solver->o_x =
    solver->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),
			mesh->x);
  
  solver->o_y =
    solver->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),
			mesh->y);
  
  solver->o_z =
    solver->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),
			mesh->z);

  solver->o_dualProjMatrix =
    solver->device.malloc(mesh->Nq*mesh->Nq*3*sizeof(dfloat),mesh->dualProjMatrix);
  
  solver->o_cubeFaceNumber =
    solver->device.malloc(mesh->Nelements*sizeof(iint),mesh->cubeFaceNumber);

  solver->o_cubeDistance =
    solver->device.malloc(mesh->Nelements*sizeof(iint),mesh->cubeDistance);
  
  solver->o_EToE =
    solver->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),mesh->EToE);

   
  if(mesh->totalHaloPairs>0){
    // copy halo element list to DEVICE
    solver->o_haloElementList =
      solver->device.malloc(mesh->totalHaloPairs*sizeof(iint), mesh->haloElementList);

    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    solver->o_haloBuffer =
      solver->device.malloc(mesh->totalHaloPairs*mesh->Np*solver->Nfields*sizeof(dfloat));
  }

  //-------------------------------------
  // NBN: 2 streams for async MPI updates
  // {Vol, Surf, update}  run on q[0]
  // {halo-get, copy} run on q[1]
  //-------------------------------------
  mesh->stream0 = solver->device.getStream();
#ifdef USE_2_STREAMS
  mesh->stream1 = solver->device.createStream();  // NBN: second stream
#else
  mesh->stream1 = mesh->stream0;                // NBN: stream1 == stream0
#endif
  solver->device.setStream(mesh->stream0);
  //-------------------------------------

  kernelInfo->addDefine("p_Nfields", solver->Nfields);
  kernelInfo->addDefine("p_N", mesh->N);
  kernelInfo->addDefine("p_Nq", mesh->Nq);
  kernelInfo->addDefine("p_Np", mesh->Np);
  kernelInfo->addDefine("p_Nfp", mesh->Nfp);
  kernelInfo->addDefine("p_Nfaces", mesh->Nfaces);
  kernelInfo->addDefine("p_Nvgeo", mesh->Nvgeo);
  kernelInfo->addDefine("p_Nsgeo", mesh->Nsgeo);
  
  int maxNodes = mymax(mesh->Nfp*mesh->Nfaces,mesh->Np);
  kernelInfo->addDefine("p_maxNodes", maxNodes);

  int NblockV = 256/mesh->Np; // works for CUDA
  kernelInfo->addDefine("p_NblockV", NblockV);

  int NblockS = 256/maxNodes; // works for CUDA
  kernelInfo->addDefine("p_NblockS", NblockS);

  // physics
  kernelInfo->addDefine("p_sqrtRT", solver->sqrtRT);
  kernelInfo->addDefine("p_invsqrtRT", (dfloat)(1./solver->sqrtRT));
  kernelInfo->addDefine("p_sqrt2", (dfloat)sqrt(2.));
  kernelInfo->addDefine("p_invsqrt2", (dfloat)sqrt(1./2.));
  kernelInfo->addDefine("p_isq12", (dfloat)sqrt(1./12.));
  kernelInfo->addDefine("p_isq6", (dfloat)sqrt(1./6.));
  kernelInfo->addDefine("p_tauInv", solver->tauInv);

  //only used by dopri, hack a default value for mrsaab
  if (solver->blockSize == 0) solver->blockSize = 1;
  kernelInfo->addDefine("p_blockSize", solver->blockSize);

  //only used by mrsaab, hack a default value for dopri
  if (solver->Nrhs == 0) solver->Nrhs =1;
  kernelInfo->addDefine("p_nrhs",solver->Nrhs);
  
  //used in dopri reduce kernels
  kernelInfo->addParserFlag("automate-add-barriers", "disabled");
  
  //set to 0 for advection
  kernelInfo->addDefine("p_invRadiusSq", 0);//1./(mesh->sphereRadius*mesh->sphereRadius));

  kernelInfo->addDefine("p_fainv", (dfloat) 0.0); // turn off rotation
  
  if(sizeof(dfloat)==4){
    kernelInfo->addDefine("dfloat","float");
    kernelInfo->addDefine("dfloat2","float2");
    kernelInfo->addDefine("dfloat4","float4");
    kernelInfo->addDefine("dfloat8","float8");
  }
  else if(sizeof(dfloat)==8){
    kernelInfo->addDefine("dfloat","double");
    kernelInfo->addDefine("dfloat2","double2");
    kernelInfo->addDefine("dfloat4","double4");
    kernelInfo->addDefine("dfloat8","double8");
  }

  if(sizeof(iint)==4){
    kernelInfo->addDefine("iint","int");
  }
  else if(sizeof(iint)==8){
    kernelInfo->addDefine("iint","long long int");
  }
  
  if(solver->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo->addCompilerFlag("--ftz=true");
    kernelInfo->addCompilerFlag("--prec-div=false");
    kernelInfo->addCompilerFlag("--prec-sqrt=false");
    kernelInfo->addCompilerFlag("--use_fast_math");
    kernelInfo->addCompilerFlag("--fmad=true"); // compiler option for cuda
  }
}
