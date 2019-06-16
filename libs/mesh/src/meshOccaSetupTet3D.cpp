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

#include "mesh.hpp"
#include "mesh3D.hpp"

void meshTet3D::OccaSetup(occa::properties &kernelInfo){

  this->mesh3D::OccaSetup(kernelInfo);

  //build inverse of mass matrix
  invMM = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  for (int n=0;n<Np*Np;n++)
    invMM[n] = MM[n];
  matrixInverse(Np,invMM);

  //set surface mass matrix
  sMT = (dfloat *) calloc(Np*Nfaces*Nfp,sizeof(dfloat));
  for (int n=0;n<Np;n++) {
    for (int m=0;m<Nfp*Nfaces;m++) {
      dfloat MSnm = 0;
      for (int i=0;i<Np;i++){
        MSnm += MM[n+i*Np]*LIFT[m+i*Nfp*Nfaces];
      }
      sMT[n+m*Np]  = MSnm;
    }
  }

  // build Dr, Ds, LIFT transposes
  dfloat *DrT = (dfloat*) calloc(Np*Np, sizeof(dfloat));
  dfloat *DsT = (dfloat*) calloc(Np*Np, sizeof(dfloat));
  dfloat *DtT = (dfloat*) calloc(Np*Np, sizeof(dfloat));
  for(int n=0;n<Np;++n){
    for(int m=0;m<Np;++m){
      DrT[n+m*Np] = Dr[n*Np+m];
      DsT[n+m*Np] = Ds[n*Np+m];
      DtT[n+m*Np] = Dt[n*Np+m];
    }
  }

  // build Dr, Ds transposes
  dfloat *DrstT = (dfloat*) calloc(3*Np*Np, sizeof(dfloat));
  for(int n=0;n<Np;++n){
    for(int m=0;m<Np;++m){
      DrstT[n+m*Np] = Dr[n*Np+m];
      DrstT[n+m*Np+Np*Np] = Ds[n*Np+m];
      DrstT[n+m*Np+2*Np*Np] = Dt[n*Np+m];
    }
  }

  dfloat *LIFTT = (dfloat*) calloc(Np*Nfaces*Nfp, sizeof(dfloat));
  for(int n=0;n<Np;++n){
    for(int m=0;m<Nfaces*Nfp;++m){
            LIFTT[n+m*Np] = LIFT[n*Nfp*Nfaces+m];
    }
  }

  // =============== BB operators [added by NC] ===============

  // deriv operators: transpose from row major to column major
  int *D0idsT = (int*) calloc(Np*4,sizeof(int));
  int *D1idsT = (int*) calloc(Np*4,sizeof(int));
  int *D2idsT = (int*) calloc(Np*4,sizeof(int));
  int *D3idsT = (int*) calloc(Np*4,sizeof(int));
  dfloat *DvalsT = (dfloat*) calloc(Np*4,sizeof(dfloat));

  int *L0idsT = (int*) calloc(Nfp*7,sizeof(int));
  dfloat *L0valsT = (dfloat*) calloc(Nfp*7,sizeof(dfloat)); // tridiag
  int *ELidsT = (int*) calloc(Np*max_EL_nnz,sizeof(int));
  dfloat *ELvalsT = (dfloat*) calloc(Np*max_EL_nnz,sizeof(dfloat));

  for (int i = 0; i < Np; ++i){
    for (int j = 0; j < 4; ++j){
      D0idsT[i+j*Np] = D0ids[j+i*4];
      D1idsT[i+j*Np] = D1ids[j+i*4];
      D2idsT[i+j*Np] = D2ids[j+i*4];
      D3idsT[i+j*Np] = D3ids[j+i*4];
      DvalsT[i+j*Np] = Dvals[j+i*4];
    }
  }

  for (int i = 0; i < Nfp; ++i){
    for (int j = 0; j < 7; ++j){
      L0idsT [i+j*Nfp] = L0ids [j+i*7];
      L0valsT[i+j*Nfp] = L0vals[j+i*7];
    }
  }

  for (int i = 0; i < Np; ++i){
    for (int j = 0; j < max_EL_nnz; ++j){
      ELidsT [i + j*Np] = ELids [j+i*max_EL_nnz];
      ELvalsT[i + j*Np] = ELvals[j+i*max_EL_nnz];
    }
  }
    // =============== end BB stuff =============================

  if(cubNp){
    dfloat *cubDrWT = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
    dfloat *cubDsWT = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
    dfloat *cubDtWT = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
    dfloat *cubDrstWT = (dfloat*) calloc(3*cubNp*Np, sizeof(dfloat));
    dfloat *cubProjectT = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
    dfloat *cubInterpT = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
    for(int n=0;n<Np;++n){
      for(int m=0;m<cubNp;++m){
        cubDrWT[n+m*Np] = cubDrW[n*cubNp+m];
        cubDsWT[n+m*Np] = cubDsW[n*cubNp+m];
        cubDtWT[n+m*Np] = cubDtW[n*cubNp+m];

        cubDrstWT[n+m*Np] = cubDrW[n*cubNp+m];
        cubDrstWT[n+m*Np+cubNp*Np] = cubDsW[n*cubNp+m];
        cubDrstWT[n+m*Np+2*cubNp*Np] = cubDtW[n*cubNp+m];

        cubProjectT[n+m*Np] = cubProject[n*cubNp+m];
        cubInterpT[m+n*cubNp] = cubInterp[m*Np+n];
      }
    }

    o_cubInterpT =
      device.malloc(Np*cubNp*sizeof(dfloat),
                          cubInterpT);

    o_cubProjectT =
      device.malloc(Np*cubNp*sizeof(dfloat),
                          cubProjectT);

    o_cubDrWT =
      device.malloc(Np*cubNp*sizeof(dfloat),
                          cubDrWT);

    o_cubDsWT =
      device.malloc(Np*cubNp*sizeof(dfloat),
                          cubDsWT);

    o_cubDtWT =
      device.malloc(Np*cubNp*sizeof(dfloat),
                          cubDtWT);

    o_cubDWmatrices = device.malloc(3*cubNp*Np*sizeof(dfloat), cubDrstWT);
  }

  if(intNfp){
    // build surface integration matrix transposes
    dfloat *intLIFTT = (dfloat*) calloc(Np*Nfaces*intNfp, sizeof(dfloat));
    dfloat *intInterpT = (dfloat*) calloc(Nfp*Nfaces*intNfp, sizeof(dfloat));
    for(int n=0;n<Np;++n){
      for(int m=0;m<Nfaces*intNfp;++m){
        intLIFTT[n+m*Np] = intLIFT[n*intNfp*Nfaces+m];
      }
    }
    for(int n=0;n<intNfp*Nfaces;++n){
      for(int m=0;m<Nfp;++m){
        intInterpT[n+m*Nfaces*intNfp] = intInterp[n*Nfp + m];
      }
    }

    o_intInterpT =
      device.malloc(Nfp*Nfaces*intNfp*sizeof(dfloat),
                          intInterpT);

    o_intLIFTT =
      device.malloc(Np*Nfaces*intNfp*sizeof(dfloat),
                          intLIFTT);

    // printf("Integration number of points: %d \n",intNfp);
    intx = (dfloat*) calloc(Nelements*Nfaces*intNfp, sizeof(dfloat));
    inty = (dfloat*) calloc(Nelements*Nfaces*intNfp, sizeof(dfloat));
    intz = (dfloat*) calloc(Nelements*Nfaces*intNfp, sizeof(dfloat));

    for(dlong e=0;e<Nelements;++e){
      for(int f=0;f<Nfaces;++f){
        for(int n=0;n<intNfp;++n){
          dfloat ix = 0, iy = 0, iz=0;
          for(int m=0;m<Nfp;++m){
            dlong vid = vmapM[m+f*Nfp+e*Nfp*Nfaces];
            dfloat xm = x[vid];
            dfloat ym = y[vid];
            dfloat zm = z[vid];
            dfloat Inm = intInterp[m+n*Nfp+f*intNfp*Nfp]; // Fixed
            ix += Inm*xm;
            iy += Inm*ym;
            iz += Inm*zm;
          }
          dlong id = n + f*intNfp + e*Nfaces*intNfp;
          intx[id] = ix;
          inty[id] = iy;
          intz[id] = iz;
        }
      }
    }

    o_intx =
      device.malloc(Nelements*Nfaces*intNfp*sizeof(dfloat),
                          intx);

    o_inty =
      device.malloc(Nelements*Nfaces*intNfp*sizeof(dfloat),
                          inty);

    o_intz =
      device.malloc(Nelements*Nfaces*intNfp*sizeof(dfloat),
                          intz);

  }

  // =============== Bernstein-Bezier allocations [added by NC] ============

  // create packed BB indexes
  o_D0ids = device.malloc(Np*4*sizeof(int),D0idsT);
  o_D1ids = device.malloc(Np*4*sizeof(int),D1idsT);
  o_D2ids = device.malloc(Np*4*sizeof(int),D2idsT);
  o_D3ids = device.malloc(Np*4*sizeof(int),D3idsT);
  o_Dvals = device.malloc(Np*4*sizeof(dfloat),DvalsT);

  unsigned char *packedDids = (unsigned char*) malloc(Np*3*4*sizeof(unsigned char));

  // for(int n=0;n<4*Np;++n){
  //   if(D1ids[n]<D0ids[n]) printf("bugger: D0id > D1id\n");
  //   if(D2ids[n]<D0ids[n]) printf("bugger: D0id > D2id\n");
  //   if(D3ids[n]<D0ids[n]) printf("bugger: D0id > D3id\n");
  // }

  for(int n=0;n<Np;++n){

    packedDids[n*4+0] = D1idsT[n+0*Np]-D0idsT[n+0*Np];
    packedDids[n*4+1] = D1idsT[n+1*Np]-D0idsT[n+1*Np];
    packedDids[n*4+2] = D1idsT[n+2*Np]-D0idsT[n+2*Np];
    packedDids[n*4+3] = D1idsT[n+3*Np]-D0idsT[n+3*Np];

    packedDids[4*Np+n*4+0] = D2idsT[n+0*Np]-D0idsT[n+0*Np];
    packedDids[4*Np+n*4+1] = D2idsT[n+1*Np]-D0idsT[n+1*Np];
    packedDids[4*Np+n*4+2] = D2idsT[n+2*Np]-D0idsT[n+2*Np];
    packedDids[4*Np+n*4+3] = D2idsT[n+3*Np]-D0idsT[n+3*Np];

    packedDids[8*Np+n*4+0] = D3idsT[n+0*Np]-D0idsT[n+0*Np];
    packedDids[8*Np+n*4+1] = D3idsT[n+1*Np]-D0idsT[n+1*Np];
    packedDids[8*Np+n*4+2] = D3idsT[n+2*Np]-D0idsT[n+2*Np];
    packedDids[8*Np+n*4+3] = D3idsT[n+3*Np]-D0idsT[n+3*Np];
  }


  o_packedDids = device.malloc(Np*3*4*sizeof(unsigned char),packedDids);

  o_L0ids  = device.malloc(Nfp*7*sizeof(int),L0idsT);
  o_L0vals = device.malloc(Nfp*7*sizeof(dfloat),L0valsT);
  o_ELids  = device.malloc(Np*max_EL_nnz*sizeof(int),ELidsT);
  o_ELvals = device.malloc(Np*max_EL_nnz*sizeof(dfloat),ELvalsT);
  // =============== end Bernstein-Bezier section [added by NC] ============

  //build element stiffness matrices
  dfloat *SrrT, *SrsT, *SrtT;
  dfloat *SsrT, *SssT, *SstT;
  dfloat *StrT, *StsT, *SttT;

  Srr = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  Srs = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  Srt = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  Ssr = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  Sss = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  Sst = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  Str = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  Sts = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  Stt = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  for (int n=0;n<Np;n++) {
    for (int m=0;m<Np;m++) {
      for (int k=0;k<Np;k++) {
        for (int l=0;l<Np;l++) {
          Srr[m+n*Np] += Dr[n+l*Np]*MM[k+l*Np]*Dr[m+k*Np];
          Srs[m+n*Np] += Dr[n+l*Np]*MM[k+l*Np]*Ds[m+k*Np];
          Srt[m+n*Np] += Dr[n+l*Np]*MM[k+l*Np]*Dt[m+k*Np];
          Ssr[m+n*Np] += Ds[n+l*Np]*MM[k+l*Np]*Dr[m+k*Np];
          Sss[m+n*Np] += Ds[n+l*Np]*MM[k+l*Np]*Ds[m+k*Np];
          Sst[m+n*Np] += Ds[n+l*Np]*MM[k+l*Np]*Dt[m+k*Np];
          Str[m+n*Np] += Dt[n+l*Np]*MM[k+l*Np]*Dr[m+k*Np];
          Sts[m+n*Np] += Dt[n+l*Np]*MM[k+l*Np]*Ds[m+k*Np];
          Stt[m+n*Np] += Dt[n+l*Np]*MM[k+l*Np]*Dt[m+k*Np];
        }
      }
    }
  }
  SrrT = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  SrsT = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  SrtT = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  SsrT = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  SssT = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  SstT = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  StrT = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  StsT = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  SttT = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  for (int n=0;n<Np;n++) {
    for (int m=0;m<Np;m++) {
      #if 0
      SrrT[m+n*Np] = Srr[n+m*Np];
      SrsT[m+n*Np] = Srs[n+m*Np];
      SrtT[m+n*Np] = Srt[n+m*Np];
      SsrT[m+n*Np] = Ssr[n+m*Np];
      SssT[m+n*Np] = Sss[n+m*Np];
      SstT[m+n*Np] = Sst[n+m*Np];
      StrT[m+n*Np] = Str[n+m*Np];
      StsT[m+n*Np] = Sts[n+m*Np];
      SttT[m+n*Np] = Stt[n+m*Np];
      #else
      SrrT[m+n*Np] = Srr[n+m*Np];
      SrsT[m+n*Np] = Srs[n+m*Np]+Ssr[n+m*Np];
      SrtT[m+n*Np] = Srt[n+m*Np]+Str[n+m*Np];
      SssT[m+n*Np] = Sss[n+m*Np];
      SstT[m+n*Np] = Sst[n+m*Np]+Sts[n+m*Np];
      SttT[m+n*Np] = Stt[n+m*Np];
      #endif
    }
  }

  dfloat *ST = (dfloat*) calloc(6*Np*Np, sizeof(dfloat));
  for(int n=0;n<Np;++n){
    for(int m=0;m<Np;++m){
      ST[n+m*Np+0*Np*Np] = Srr[n*Np+m];
      ST[n+m*Np+1*Np*Np] = Srs[n*Np+m]+Ssr[n*Np+m];
      ST[n+m*Np+2*Np*Np] = Srt[n*Np+m]+Str[n*Np+m];
      ST[n+m*Np+3*Np*Np] = Sss[n*Np+m];
      ST[n+m*Np+4*Np*Np] = Sst[n*Np+m]+Sts[n*Np+m];
      ST[n+m*Np+5*Np*Np] = Stt[n*Np+m];
    }
  }


  o_Dr = device.malloc(Np*Np*sizeof(dfloat), Dr);
  o_Ds = device.malloc(Np*Np*sizeof(dfloat), Ds);
  o_Dt = device.malloc(Np*Np*sizeof(dfloat), Dt);

  o_DrT = device.malloc(Np*Np*sizeof(dfloat), DrT);
  o_DsT = device.malloc(Np*Np*sizeof(dfloat), DsT);
  o_DtT = device.malloc(Np*Np*sizeof(dfloat), DtT);

  o_Dmatrices = device.malloc(3*Np*Np*sizeof(dfloat), DrstT);

  o_MM = device.malloc(Np*Np*sizeof(dfloat), MM);

  o_sMT = device.malloc(Np*Nfaces*Nfp*sizeof(dfloat), sMT);

  o_LIFT =
    device.malloc(Np*Nfaces*Nfp*sizeof(dfloat),
                        LIFT);

  o_LIFTT =
    device.malloc(Np*Nfaces*Nfp*sizeof(dfloat),
                        LIFTT);

  o_MM =
    device.malloc(Np*Np*sizeof(dfloat),
      MM);

  o_vgeo =
    device.malloc(Nelements*Nvgeo*sizeof(dfloat),
                        vgeo);

  o_sgeo =
    device.malloc(Nelements*Nfaces*Nsgeo*sizeof(dfloat),
                        sgeo);

  o_cubsgeo = o_sgeo; //dummy cubature geo factors

  o_ggeo =
    device.malloc(Nelements*Nggeo*sizeof(dfloat),
      ggeo);

  o_cubvgeo =   device.malloc(sizeof(dfloat));// dummy

  o_SrrT = device.malloc(Np*Np*sizeof(dfloat), SrrT);
  o_SrsT = device.malloc(Np*Np*sizeof(dfloat), SrsT);
  o_SrtT = device.malloc(Np*Np*sizeof(dfloat), SrtT);
  o_SsrT = device.malloc(Np*Np*sizeof(dfloat), SsrT);
  o_SssT = device.malloc(Np*Np*sizeof(dfloat), SssT);
  o_SstT = device.malloc(Np*Np*sizeof(dfloat), SstT);
  o_StrT = device.malloc(Np*Np*sizeof(dfloat), StrT);
  o_StsT = device.malloc(Np*Np*sizeof(dfloat), StsT);
  o_SttT = device.malloc(Np*Np*sizeof(dfloat), SttT);

  o_Smatrices = device.malloc(6*Np*Np*sizeof(dfloat), ST);

  free(DrstT); free(ST);

  free(D0idsT);
  free(D1idsT);
  free(D2idsT);
  free(D3idsT);
  free(DvalsT);

  free(L0idsT);
  free(L0valsT);
  free(ELidsT);
  free(ELvalsT);
}
