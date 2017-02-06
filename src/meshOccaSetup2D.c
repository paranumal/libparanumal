
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include  "mpi.h"

#include "mesh2D.h"

void meshOccaSetup2D(mesh2D *mesh, char *deviceConfig, occa::kernelInfo &kernelInfo){
  
  mesh->device.setup(deviceConfig);
  
  // build Dr, Ds, LIFT transposes
  dfloat *DrT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *DsT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Np;++m){
      DrT[n+m*mesh->Np] = mesh->Dr[n*mesh->Np+m];
      DsT[n+m*mesh->Np] = mesh->Ds[n*mesh->Np+m];
    }
  }

  dfloat *LIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Nfaces*mesh->Nfp;++m){
      LIFTT[n+m*mesh->Np] = mesh->LIFT[n*mesh->Nfp*mesh->Nfaces+m];
    }
  }

  // build volume cubature matrix transposes
  iint cubNpBlocked = mesh->Np*((mesh->cubNp+mesh->Np-1)/mesh->Np);
  dfloat *cubDrWT = (dfloat*) calloc(cubNpBlocked*mesh->Np, sizeof(dfloat));
  dfloat *cubDsWT = (dfloat*) calloc(cubNpBlocked*mesh->Np, sizeof(dfloat));
  dfloat *cubProjectT = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  dfloat *cubInterpT = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->cubNp;++m){
      cubDrWT[n+m*mesh->Np] = mesh->cubDrW[n*mesh->cubNp+m];
      cubDsWT[n+m*mesh->Np] = mesh->cubDsW[n*mesh->cubNp+m];


      cubProjectT[n+m*mesh->Np] = mesh->cubProject[n*mesh->cubNp+m];
      cubInterpT[m+n*mesh->cubNp] = mesh->cubInterp[m*mesh->Np+n];
      //      printf("%g @ ", cubInterpT[m+n*mesh->cubNp]);
    }
  }

  // build surface integration matrix transposes
  dfloat *intLIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  dfloat *intInterpT = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Nfaces*mesh->intNfp;++m){
      intLIFTT[n+m*mesh->Np] = mesh->intLIFT[n*mesh->intNfp*mesh->Nfaces+m];
    }
  }
  for(int n=0;n<mesh->intNfp*mesh->Nfaces;++n){
    for(int m=0;m<mesh->Nfp;++m){
      intInterpT[n+m*mesh->Nfaces*mesh->intNfp] = mesh->intInterp[n*mesh->Nfp + m];
    }
  }

  // =============== BB operators [added by JC] ===============
#if 1

  // deriv operators: transpose from row major to column major
  iint *D1ids = (iint*) calloc(mesh->Np*3,sizeof(iint));
  iint *D2ids = (iint*) calloc(mesh->Np*3,sizeof(iint));
  iint *D3ids = (iint*) calloc(mesh->Np*3,sizeof(iint));
  dfloat *Dvals = (dfloat*) calloc(mesh->Np*3,sizeof(dfloat));  

  dfloat *VBq = (dfloat*) calloc(mesh->Np*mesh->cubNp,sizeof(dfloat));
  dfloat *PBq = (dfloat*) calloc(mesh->Np*mesh->cubNp,sizeof(dfloat));

  dfloat *L0vals = (dfloat*) calloc(mesh->Nfp*3,sizeof(dfloat)); // tridiag
  iint *ELids = (iint*) calloc(1+mesh->Np*mesh->max_EL_nnz,sizeof(iint));
  dfloat *ELvals = (dfloat*) calloc(1+mesh->Np*mesh->max_EL_nnz,sizeof(dfloat));
  
  if(mesh->max_EL_nnz){
    
    for (iint i = 0; i < mesh->Np; ++i){
      for (iint j = 0; j < 3; ++j){
	D1ids[i+j*mesh->Np] = mesh->D1ids[j+i*3];
	D2ids[i+j*mesh->Np] = mesh->D2ids[j+i*3];
	D3ids[i+j*mesh->Np] = mesh->D3ids[j+i*3];      
	Dvals[i+j*mesh->Np] = mesh->Dvals[j+i*3];    
      }
    }
    
    for (iint i = 0; i < mesh->cubNp; ++i){
      for (iint j = 0; j < mesh->Np; ++j){
	VBq[i+j*mesh->cubNp] = mesh->VBq[j+i*mesh->Np];
	PBq[j+i*mesh->Np] = mesh->PBq[i+j*mesh->cubNp];
      }
    }
    
    
    for (iint i = 0; i < mesh->Nfp; ++i){
      for (iint j = 0; j < 3; ++j){
	L0vals[i+j*mesh->Nfp] = mesh->L0vals[j+i*3];
      }
    }
    
    for (iint i = 0; i < mesh->Np; ++i){
      for (iint j = 0; j < mesh->max_EL_nnz; ++j){
	ELids[i + j*mesh->Np] = mesh->ELids[j+i*mesh->max_EL_nnz];
	ELvals[i + j*mesh->Np] = mesh->ELvals[j+i*mesh->max_EL_nnz]; // ???
      }
    }
  }
#endif  
  // =============== end BB stuff =============================


  
  mesh->intx = (dfloat*) calloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  mesh->inty = (dfloat*) calloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint f=0;f<mesh->Nfaces;++f){
      for(iint n=0;n<mesh->intNfp;++n){
	dfloat ix = 0, iy = 0;
	for(iint m=0;m<mesh->Nfp;++m){
	  iint vid = mesh->vmapM[m+f*mesh->Nfp+e*mesh->Nfp*mesh->Nfaces];
	  dfloat xm = mesh->x[vid];
	  dfloat ym = mesh->y[vid];
	  dfloat Inm = mesh->intInterp[n+f*mesh->intNfp+m*mesh->intNfp*mesh->Nfaces];
	  ix += Inm*xm;
	  iy += Inm*ym;
	}
	iint id = n + f*mesh->intNfp + e*mesh->Nfaces*mesh->intNfp;
	mesh->intx[id] = ix;
	mesh->inty[id] = iy;
      }
    }
  }
  
  
  // find elements that have all neighbors on this process
  iint *internalElementIds = (iint*) calloc(mesh->Nelements, sizeof(iint));
  iint *notInternalElementIds = (iint*) calloc(mesh->Nelements, sizeof(iint));

  iint Ninterior = 0, NnotInterior = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    iint flag = 0;
    for(iint f=0;f<mesh->Nfaces;++f)
      if(mesh->EToP[e*mesh->Nfaces+f]!=-1)
	flag = 1;
    if(!flag)
      internalElementIds[Ninterior++] = e;
    else
      notInternalElementIds[NnotInterior++] = e;
  }

  printf("NinteriorElements = %d, NnotInternalElements = %d\n", Ninterior, NnotInterior);
  
  mesh->NinternalElements = Ninterior;
  mesh->NnotInternalElements = NnotInterior;
  if(Ninterior)
    mesh->o_internalElementIds    = mesh->device.malloc(Ninterior*sizeof(iint), internalElementIds);
  
  if(NnotInterior>0)
    mesh->o_notInternalElementIds = mesh->device.malloc(NnotInterior*sizeof(iint), notInternalElementIds);
  
  // OCCA allocate device memory (remember to go back for halo)
  mesh->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);
  mesh->o_rhsq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);
  mesh->o_resq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->resq);


  mesh->o_D = mesh->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D);
  

  mesh->o_Dr = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
				   mesh->Dr);

  mesh->o_Ds = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
				   mesh->Ds);

  mesh->o_DrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
				    DrT);

  mesh->o_DsT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
				    DsT);

  mesh->o_LIFT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
			mesh->LIFT);

  mesh->o_LIFTT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
			LIFTT);

  // HARD CODED FOR QUADS
  mesh->o_vgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nvgeo*mesh->Np*sizeof(dfloat),
			mesh->vgeo);
  
  // HARD CODED FOR QUADS
  mesh->o_sgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->Nfp*mesh->Nsgeo*sizeof(dfloat),
			mesh->sgeo);

  mesh->o_ggeo =
    mesh->device.malloc(mesh->Nelements*mesh->Np*mesh->Nggeo*sizeof(dfloat),
			mesh->ggeo);
  

  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
			mesh->vmapM);

  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
			mesh->vmapP);

  mesh->o_EToB =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),
			mesh->EToB);

  mesh->o_x =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),
			mesh->x);

  mesh->o_y =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),
			mesh->y);


  mesh->o_intx =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
			mesh->intx);

  mesh->o_inty =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
			mesh->inty);

  
  if(mesh->totalHaloPairs>0){
    // copy halo element list to DEVICE
    mesh->o_haloElementList =
      mesh->device.malloc(mesh->totalHaloPairs*sizeof(iint), mesh->haloElementList);
    
    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat));
  }

  mesh->o_cubInterpT =
    mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),
			cubInterpT);
  
  mesh->o_cubProjectT =
    mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),
			cubProjectT);

  mesh->o_cubDrWT =
    mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),
			cubDrWT);
  
  mesh->o_cubDsWT =
    mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),
			cubDsWT);

  mesh->o_intInterpT =
    mesh->device.malloc(mesh->Nfp*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
			intInterpT);

  mesh->o_intLIFTT =
    mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
			intLIFTT);
  
  // =============== Bernstein-Bezier allocations [added by JC] ============
#if 1
  if(mesh->max_EL_nnz){
  mesh->o_D1ids = mesh->device.malloc(mesh->Np*3*sizeof(iint),D1ids);
  mesh->o_D2ids = mesh->device.malloc(mesh->Np*3*sizeof(iint),D2ids);
  mesh->o_D3ids = mesh->device.malloc(mesh->Np*3*sizeof(iint),D3ids);
  mesh->o_Dvals = mesh->device.malloc(mesh->Np*3*sizeof(dfloat),Dvals);

  mesh->o_VBq = mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),VBq);
  mesh->o_PBq = mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),PBq);

  mesh->o_L0vals = mesh->device.malloc(mesh->Nfp*3*sizeof(dfloat),L0vals);
  mesh->o_ELids =
    mesh->device.malloc(mesh->Np*mesh->max_EL_nnz*sizeof(iint),ELids);
  mesh->o_ELvals =
    mesh->device.malloc(mesh->Np*mesh->max_EL_nnz*sizeof(dfloat),ELvals);
  }
#endif
  // =============== end Bernstein-Bezier section [added by JC] ============  


  
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
  kernelInfo.addDefine("p_SJID", SJID);
  kernelInfo.addDefine("p_IJID", IJID);
  kernelInfo.addDefine("p_WSJID", WSJID);
  kernelInfo.addDefine("p_IHID", IHID);
    
  kernelInfo.addDefine("p_max_EL_nnz", mesh->max_EL_nnz); // for Bernstein Bezier lift

  kernelInfo.addDefine("p_cubNp", mesh->cubNp);
  kernelInfo.addDefine("p_intNfp", mesh->intNfp);
  kernelInfo.addDefine("p_intNfpNfaces", mesh->intNfp*mesh->Nfaces);

  if(sizeof(dfloat)==4){
    kernelInfo.addDefine("dfloat","float");
    kernelInfo.addDefine("dfloat4","float4");
    kernelInfo.addDefine("dfloat8","float8");
  }
  if(sizeof(dfloat)==8){
    kernelInfo.addDefine("dfloat","double");
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


  kernelInfo.addDefine("p_G00ID", G00ID);
  kernelInfo.addDefine("p_G01ID", G01ID);
  kernelInfo.addDefine("p_G11ID", G11ID);
  kernelInfo.addDefine("p_GWJID", GWJID);


  kernelInfo.addDefine("p_RXID", RXID);
  kernelInfo.addDefine("p_SXID", SXID);

  kernelInfo.addDefine("p_RYID", RYID);
  kernelInfo.addDefine("p_SYID", SYID);

  kernelInfo.addDefine("p_JWID", JWID);


  
}
