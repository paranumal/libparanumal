
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include  "mpi.h"

#include "mesh3D.h"

void boltzmannOccaSetupQuad3D(mesh_t *mesh, char *deviceConfig, occa::kernelInfo &kernelInfo){


  mesh->device.setup(deviceConfig);

  occa::initTimer(mesh->device);

  dfloat *LIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Nfaces*mesh->Nfp;++m){
      LIFTT[n+m*mesh->Np] = mesh->LIFT[n*mesh->Nfp*mesh->Nfaces+m];
    }
  }
  // OCCA allocate device memory (remember to go back for halo)
  mesh->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);

  mesh->o_rhsq =
    mesh->device.malloc(mesh->Np*mesh->Nrhs*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);
  mesh->o_resq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->resq);

  mesh->o_D = mesh->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D);

  mesh->o_LIFT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
        mesh->LIFT);

  mesh->o_LIFTT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
        LIFTT);

  mesh->o_vgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nvgeo*mesh->Np*sizeof(dfloat),
			mesh->vgeo);
  mesh->o_sgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->Nfp*mesh->Nsgeo*sizeof(dfloat),
			mesh->sgeo);
  
  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
        mesh->vmapM);

  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
        mesh->vmapP);


  mesh->o_mapP =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
        mesh->mapP);

  mesh->o_EToB =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),
        mesh->EToB);

  mesh->o_x =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),
        mesh->x);

  mesh->o_y =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),
        mesh->y);


  mesh->o_z =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),
        mesh->z);

  if(mesh->totalHaloPairs>0){
    // copy halo element list to DEVICE
    mesh->o_haloElementList =
      mesh->device.malloc(mesh->totalHaloPairs*sizeof(iint), mesh->haloElementList);

    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat));
  }

  //-------------------------------------
  // NBN: 2 streams for async MPI updates
  // {Vol, Surf, update}  run on q[0]
  // {halo-get, copy} run on q[1]
  //-------------------------------------
  mesh->stream0 = mesh->device.getStream();
#ifdef USE_2_STREAMS
  mesh->stream1 = mesh->device.createStream();  // NBN: second stream
#else
  mesh->stream1 = mesh->stream0;                // NBN: stream1 == stream0
#endif
  mesh->device.setStream(mesh->stream0);
  //-------------------------------------

  kernelInfo.addDefine("p_Nfields", mesh->Nfields);
  kernelInfo.addDefine("p_N", mesh->N);
  kernelInfo.addDefine("p_Nq", mesh->N+1);
  kernelInfo.addDefine("p_Np", mesh->Np);
  kernelInfo.addDefine("p_Nfp", mesh->Nfp);
  kernelInfo.addDefine("p_Nfaces", mesh->Nfaces);
  kernelInfo.addDefine("p_NfacesNfp", mesh->Nfp*mesh->Nfaces);
  kernelInfo.addDefine("p_Nvgeo", mesh->Nvgeo);
  kernelInfo.addDefine("p_Nsgeo", mesh->Nsgeo);
  kernelInfo.addDefine("p_Nggeo", mesh->Nggeo);

  kernelInfo.addDefine("p_NXID", NXID);
  kernelInfo.addDefine("p_NYID", NYID);
  kernelInfo.addDefine("p_NZID", NZID);
  kernelInfo.addDefine("p_SJID", SJID);
  kernelInfo.addDefine("p_IJID", IJID);
  kernelInfo.addDefine("p_WSJID", WSJID);
  kernelInfo.addDefine("p_IHID", IHID);

  kernelInfo.addDefine("p_max_EL_nnz", mesh->max_EL_nnz); // for Bernstein Bezier lift

  printf("adding p_cubNp\n");
  kernelInfo.addDefine("p_cubNp", mesh->cubNp);
  kernelInfo.addDefine("p_intNfp", mesh->intNfp);
  kernelInfo.addDefine("p_intNfpNfaces", mesh->intNfp*mesh->Nfaces);

  if(sizeof(dfloat)==4){
    kernelInfo.addDefine("dfloat","float");
    kernelInfo.addDefine("dfloat2","float2");
    kernelInfo.addDefine("dfloat4","float4");
    kernelInfo.addDefine("dfloat8","float8");
  }
  if(sizeof(dfloat)==8){
    kernelInfo.addDefine("dfloat","double");
    kernelInfo.addDefine("dfloat2","double2");
    kernelInfo.addDefine("dfloat4","double4");
    kernelInfo.addDefine("dfloat8","double8");
  }

  if(sizeof(iint)==4){
    kernelInfo.addDefine("iint","int");
  }
  if(sizeof(iint)==8){
    kernelInfo.addDefine("iint","long long int");
  }

  if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo.addCompilerFlag("--ftz=true");
    kernelInfo.addCompilerFlag("--prec-div=false");
    kernelInfo.addCompilerFlag("--prec-sqrt=false");
    kernelInfo.addCompilerFlag("--use_fast_math");
    kernelInfo.addCompilerFlag("--fmad=true"); // compiler option for cuda
  }
  
  kernelInfo.addDefine("p_RXID", RXID);
  kernelInfo.addDefine("p_SXID", SXID);
  kernelInfo.addDefine("p_TXID", TXID);

  kernelInfo.addDefine("p_RYID", RYID);
  kernelInfo.addDefine("p_SYID", SYID);
  kernelInfo.addDefine("p_TYID", TYID);

  kernelInfo.addDefine("p_RZID", RZID);
  kernelInfo.addDefine("p_SZID", SZID);
  kernelInfo.addDefine("p_TZID", TZID);
  
  kernelInfo.addDefine("p_JID", JID);
  kernelInfo.addDefine("p_JWID", JWID);
}
