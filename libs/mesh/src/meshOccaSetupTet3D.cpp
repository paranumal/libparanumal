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

void meshTet3D::OccaSetup(){

  this->mesh3D::OccaSetup();

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
    device.malloc((Nelements+totalHaloPairs)*Nvgeo*sizeof(dfloat),
                        vgeo);

  o_sgeo =
    device.malloc(Nelements*Nfaces*Nsgeo*sizeof(dfloat),
                        sgeo);

  o_ggeo =
    device.malloc(Nelements*Nggeo*sizeof(dfloat),
      ggeo);

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
}
