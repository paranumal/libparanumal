
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include  "mpi.h"

#include "mesh3D.h"

void meshOccaSetup3D(mesh3D *mesh, char *deviceConfig, occa::kernelInfo &kernelInfo){

  mesh->device.setup(deviceConfig);

  occa::initTimer(mesh->device);

  // OCCA allocate device memory (remember to go back for halo)
  mesh->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);
  mesh->o_rhsq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);
  mesh->o_resq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->resq);

  if(mesh->Nfaces!=4){
    mesh->o_D = mesh->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D);
  }

  if(mesh->Nfaces==4){

    // build Dr, Ds, LIFT transposes
    dfloat *DrT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
    dfloat *DsT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
    dfloat *DtT = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
      	DrT[n+m*mesh->Np] = mesh->Dr[n*mesh->Np+m];
      	DsT[n+m*mesh->Np] = mesh->Ds[n*mesh->Np+m];
      	DtT[n+m*mesh->Np] = mesh->Dt[n*mesh->Np+m];
      }
    }

    dfloat *LIFTT = (dfloat*) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Nfaces*mesh->Nfp;++m){
	      LIFTT[n+m*mesh->Np] = mesh->LIFT[n*mesh->Nfp*mesh->Nfaces+m];
      }
    }

    // =============== BB operators [added by NC] ===============
    
    // deriv operators: transpose from row major to column major
    iint *D0ids = (iint*) calloc(mesh->Np*4,sizeof(iint));
    iint *D1ids = (iint*) calloc(mesh->Np*4,sizeof(iint));
    iint *D2ids = (iint*) calloc(mesh->Np*4,sizeof(iint));
    iint *D3ids = (iint*) calloc(mesh->Np*4,sizeof(iint));
    dfloat *Dvals = (dfloat*) calloc(mesh->Np*4,sizeof(dfloat));  

    iint *L0ids = (iint*) calloc(mesh->Nfp*7,sizeof(iint)); 
    dfloat *L0vals = (dfloat*) calloc(mesh->Nfp*7,sizeof(dfloat)); // tridiag
    iint *ELids = (iint*) calloc(mesh->Np*mesh->max_EL_nnz,sizeof(iint));
    dfloat *ELvals = (dfloat*) calloc(mesh->Np*mesh->max_EL_nnz,sizeof(dfloat));
    
    for (iint i = 0; i < mesh->Np; ++i){
      for (iint j = 0; j < 4; ++j){
        D0ids[i+j*mesh->Np] = mesh->D0ids[j+i*4];
        D1ids[i+j*mesh->Np] = mesh->D1ids[j+i*4];
        D2ids[i+j*mesh->Np] = mesh->D2ids[j+i*4];
        D3ids[i+j*mesh->Np] = mesh->D3ids[j+i*4];      
        Dvals[i+j*mesh->Np] = mesh->Dvals[j+i*4];    
      }
    }
    
    for (iint i = 0; i < mesh->Nfp; ++i){
      for (iint j = 0; j < 7; ++j){
        L0ids [i+j*mesh->Nfp] = mesh->L0ids [j+i*7];
        L0vals[i+j*mesh->Nfp] = mesh->L0vals[j+i*7];
      }
    }
    
    for (iint i = 0; i < mesh->Np; ++i){
      for (iint j = 0; j < mesh->max_EL_nnz; ++j){
        ELids [i + j*mesh->Np] = mesh->ELids [j+i*mesh->max_EL_nnz];
        ELvals[i + j*mesh->Np] = mesh->ELvals[j+i*mesh->max_EL_nnz]; 
      }
    } 
    // =============== end BB stuff =============================
    
    mesh->o_Dr = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), mesh->Dr);
    mesh->o_Ds = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), mesh->Ds);
    mesh->o_Dt = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), mesh->Dt);

    mesh->o_DrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), DrT);
    mesh->o_DsT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), DsT);
    mesh->o_DtT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), DtT);

    mesh->o_LIFT =
      mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
			  mesh->LIFT);

    mesh->o_LIFTT =
      mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat),
			  LIFTT);

    dfloat *cubDrWT = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
    dfloat *cubDsWT = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
    dfloat *cubDtWT = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
    dfloat *cubProjectT = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
    dfloat *cubInterpT = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->cubNp;++m){
      	cubDrWT[n+m*mesh->Np] = mesh->cubDrW[n*mesh->cubNp+m];
      	cubDsWT[n+m*mesh->Np] = mesh->cubDsW[n*mesh->cubNp+m];
      	cubDtWT[n+m*mesh->Np] = mesh->cubDtW[n*mesh->cubNp+m];
      	
      	cubProjectT[n+m*mesh->Np] = mesh->cubProject[n*mesh->cubNp+m];
      	cubInterpT[m+n*mesh->cubNp] = mesh->cubInterp[m*mesh->Np+n];
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

    mesh->o_cubDtWT =
      mesh->device.malloc(mesh->Np*mesh->cubNp*sizeof(dfloat),
    	  cubDtWT);

    mesh->o_intInterpT =
      mesh->device.malloc(mesh->Nfp*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
        intInterpT);

    mesh->o_intLIFTT =
      mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
        intLIFTT);


    // =============== Bernstein-Bezier allocations [added by NC] ============
    mesh->o_D0ids = mesh->device.malloc(mesh->Np*4*sizeof(iint),D0ids);
    mesh->o_D1ids = mesh->device.malloc(mesh->Np*4*sizeof(iint),D1ids);
    mesh->o_D2ids = mesh->device.malloc(mesh->Np*4*sizeof(iint),D2ids);
    mesh->o_D3ids = mesh->device.malloc(mesh->Np*4*sizeof(iint),D3ids);
    mesh->o_Dvals = mesh->device.malloc(mesh->Np*4*sizeof(dfloat),Dvals);

    mesh->o_L0ids  = mesh->device.malloc(mesh->Nfp*7*sizeof(iint),L0ids);
    mesh->o_L0vals = mesh->device.malloc(mesh->Nfp*7*sizeof(dfloat),L0vals);
    mesh->o_ELids  = mesh->device.malloc(mesh->Np*mesh->max_EL_nnz*sizeof(iint),ELids);
    mesh->o_ELvals = mesh->device.malloc(mesh->Np*mesh->max_EL_nnz*sizeof(dfloat),ELvals);
    // =============== end Bernstein-Bezier section [added by NC] ============  
  }

  //build element stiffness matrices
  dfloat *SrrT, *SrsT, *SrtT;
  dfloat *SsrT, *SssT, *SstT;
  dfloat *StrT, *StsT, *SttT;
  if (mesh->Nverts==4) {
    mesh->Srr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Srs = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Srt = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Ssr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Sss = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Sst = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Str = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Sts = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    mesh->Stt = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    for (iint n=0;n<mesh->Np;n++) {
      for (iint m=0;m<mesh->Np;m++) {
        for (iint k=0;k<mesh->Np;k++) {
          for (iint l=0;l<mesh->Np;l++) {
            mesh->Srr[m+n*mesh->Np] += mesh->Dr[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dr[m+k*mesh->Np];
            mesh->Srs[m+n*mesh->Np] += mesh->Dr[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Ds[m+k*mesh->Np];
            mesh->Srt[m+n*mesh->Np] += mesh->Dr[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dt[m+k*mesh->Np];
            mesh->Ssr[m+n*mesh->Np] += mesh->Ds[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dr[m+k*mesh->Np];
            mesh->Sss[m+n*mesh->Np] += mesh->Ds[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Ds[m+k*mesh->Np];
            mesh->Sst[m+n*mesh->Np] += mesh->Ds[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dt[m+k*mesh->Np];
            mesh->Str[m+n*mesh->Np] += mesh->Dt[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dr[m+k*mesh->Np];
            mesh->Sts[m+n*mesh->Np] += mesh->Dt[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Ds[m+k*mesh->Np];
            mesh->Stt[m+n*mesh->Np] += mesh->Dt[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dt[m+k*mesh->Np];
          }
        } 
      }
    }
    SrrT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    SrsT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    SrtT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    SsrT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    SssT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    SstT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    StrT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    StsT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    SttT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    for (iint n=0;n<mesh->Np;n++) {
      for (iint m=0;m<mesh->Np;m++) {  
        SrrT[m+n*mesh->Np] = mesh->Srr[n+m*mesh->Np];
        SrsT[m+n*mesh->Np] = mesh->Srs[n+m*mesh->Np];
        SrtT[m+n*mesh->Np] = mesh->Srt[n+m*mesh->Np];
        SsrT[m+n*mesh->Np] = mesh->Ssr[n+m*mesh->Np];
        SssT[m+n*mesh->Np] = mesh->Sss[n+m*mesh->Np];
        SstT[m+n*mesh->Np] = mesh->Sst[n+m*mesh->Np];
        StrT[m+n*mesh->Np] = mesh->Str[n+m*mesh->Np];
        StsT[m+n*mesh->Np] = mesh->Sts[n+m*mesh->Np];
        SttT[m+n*mesh->Np] = mesh->Stt[n+m*mesh->Np];
      }
    }
  }

 #if 1
 // printf("Integration number of points: %d \n",mesh->intNfp);
  mesh->intx = (dfloat*) calloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  mesh->inty = (dfloat*) calloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));
  mesh->intz = (dfloat*) calloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp, sizeof(dfloat));

  for(iint e=0;e<mesh->Nelements;++e){
    for(iint f=0;f<mesh->Nfaces;++f){
      for(iint n=0;n<mesh->intNfp;++n){
        dfloat ix = 0, iy = 0, iz=0;
        for(iint m=0;m<mesh->Nfp;++m){
          iint vid = mesh->vmapM[m+f*mesh->Nfp+e*mesh->Nfp*mesh->Nfaces];
          dfloat xm = mesh->x[vid];
          dfloat ym = mesh->y[vid];
          dfloat zm = mesh->z[vid];
          dfloat Inm = mesh->intInterp[m+n*mesh->Nfp+f*mesh->intNfp*mesh->Nfp]; // Fixed
          ix += Inm*xm;
          iy += Inm*ym;
          iz += Inm*zm;
        }
        iint id = n + f*mesh->intNfp + e*mesh->Nfaces*mesh->intNfp;
        mesh->intx[id] = ix;
        mesh->inty[id] = iy;
        mesh->intz[id] = iz;
      }
    }
  }




  #endif

  printf("Nverts = %d, Nfaces = %d\n",mesh->Nverts,mesh->Nfaces);
  if (mesh->Nverts==8){     // hardcoded for hexes

    mesh->o_vgeo =
      mesh->device.malloc(mesh->Nelements*mesh->Np*mesh->Nvgeo*sizeof(dfloat),
                          mesh->vgeo);

    mesh->o_sgeo =
      mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->Nfp*mesh->Nsgeo*sizeof(dfloat),
                          mesh->sgeo);

    mesh->o_ggeo =
      mesh->device.malloc(mesh->Nelements*mesh->Np*mesh->Nggeo*sizeof(dfloat),
        mesh->ggeo);

  }else if (mesh->Nverts==4){     // for tets
    mesh->o_MM =
      mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),
        mesh->MM);

    mesh->o_vgeo =
      mesh->device.malloc(mesh->Nelements*mesh->Nvgeo*sizeof(dfloat),
                          mesh->vgeo);

    mesh->o_sgeo =
      mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->Nsgeo*sizeof(dfloat),
                          mesh->sgeo);

    mesh->o_ggeo =
      mesh->device.malloc(mesh->Nelements*mesh->Nggeo*sizeof(dfloat),
        mesh->ggeo);

    mesh->o_SrrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SrrT);
    mesh->o_SrsT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SrsT);
    mesh->o_SrtT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SrtT);
    mesh->o_SsrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SsrT);
    mesh->o_SssT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SssT);
    mesh->o_SstT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SstT);
    mesh->o_StrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), StrT);
    mesh->o_StsT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), StsT);
    mesh->o_SttT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SttT);

  }else{
    printf("Nverts = %d: unknown element type!\n",mesh->Nverts);
  }

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
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), mesh->x);

  mesh->o_y =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), mesh->y);

  mesh->o_z =
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), mesh->z);


  mesh->o_intx =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
      mesh->intx);

  mesh->o_inty =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
      mesh->inty);

   mesh->o_intz =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->intNfp*sizeof(dfloat),
      mesh->intz);

  if(mesh->totalHaloPairs>0){
    // copy halo element list to DEVICE
    mesh->o_haloElementList =
      mesh->device.malloc(mesh->totalHaloPairs*sizeof(iint), mesh->haloElementList);

    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat));
  }

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

  kernelInfo.addDefine("p_max_EL_nnz", mesh->max_EL_nnz); // for Bernstein Bezier lift

  kernelInfo.addDefine("p_NXID", NXID);
  kernelInfo.addDefine("p_NYID", NYID);
  kernelInfo.addDefine("p_NZID", NZID);
  kernelInfo.addDefine("p_SJID", SJID);
  kernelInfo.addDefine("p_IJID", IJID);
  kernelInfo.addDefine("p_WSJID", WSJID);
  kernelInfo.addDefine("p_IHID", IHID);


  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 256/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 256/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  kernelInfo.addDefine("p_Lambda2", 0.5f);  

  kernelInfo.addDefine("p_cubNp", mesh->cubNp);
  kernelInfo.addDefine("p_intNfp", mesh->intNfp);
  kernelInfo.addDefine("p_intNfpNfaces", mesh->intNfp*mesh->Nfaces);

  if(sizeof(dfloat)==4){
    kernelInfo.addDefine("dfloat","float");
    kernelInfo.addDefine("dfloat4","float4");
  }
  if(sizeof(dfloat)==8){
    kernelInfo.addDefine("dfloat","double");
    kernelInfo.addDefine("dfloat4","double4");
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
  kernelInfo.addDefine("p_G02ID", G02ID);
  kernelInfo.addDefine("p_G11ID", G11ID);
  kernelInfo.addDefine("p_G12ID", G12ID);
  kernelInfo.addDefine("p_G22ID", G22ID);
  kernelInfo.addDefine("p_GWJID", GWJID);


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
