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

void meshTri3D::OccaSetup(){

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
  for(int n=0;n<Np;++n){
    for(int m=0;m<Np;++m){
      DrT[n+m*Np] = Dr[n*Np+m];
      DsT[n+m*Np] = Ds[n*Np+m];
    }
  }

  // build Dr, Ds transposes
  dfloat *DrsT = (dfloat*) calloc(2*Np*Np, sizeof(dfloat));
  for(int n=0;n<Np;++n){
    for(int m=0;m<Np;++m){
      DrsT[n+m*Np] = Dr[n*Np+m];
      DrsT[n+m*Np+Np*Np] = Ds[n*Np+m];
    }
  }

  dfloat *LIFTT = (dfloat*) calloc(Np*Nfaces*Nfp, sizeof(dfloat));
  for(int n=0;n<Np;++n){
    for(int m=0;m<Nfaces*Nfp;++m){
      LIFTT[n+m*Np] = LIFT[n*Nfp*Nfaces+m];
    }
  }

  // build volume cubature matrix transposes
  int cubNpBlocked = Np*((cubNp+Np-1)/Np);
  dfloat *cubDrWT = (dfloat*) calloc(cubNpBlocked*Np, sizeof(dfloat));
  dfloat *cubDsWT = (dfloat*) calloc(cubNpBlocked*Np, sizeof(dfloat));
  dfloat *cubDrsWT = (dfloat*) calloc(2*cubNp*Np, sizeof(dfloat));
  dfloat *cubProjectT = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
  dfloat *cubInterpT = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
  for(int n=0;n<Np;++n){
    for(int m=0;m<cubNp;++m){
      cubDrWT[n+m*Np] = cubDrW[n*cubNp+m];
      cubDsWT[n+m*Np] = cubDsW[n*cubNp+m];

      cubDrsWT[n+m*Np] = cubDrW[n*cubNp+m];
      cubDrsWT[n+m*Np+cubNp*Np] = cubDsW[n*cubNp+m];

      cubProjectT[n+m*Np] = cubProject[n*cubNp+m];
      cubInterpT[m+n*cubNp] = cubInterp[m*Np+n];
    }
  }

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

  //build element stiffness matrices
  dfloat *SrrT, *SrsT, *SsrT, *SssT;
  Srr = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  Srs = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  Ssr = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  Sss = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  for (int n=0;n<Np;n++) {
    for (int m=0;m<Np;m++) {
      for (int k=0;k<Np;k++) {
        for (int l=0;l<Np;l++) {
          Srr[m+n*Np] += Dr[n+l*Np]*MM[k+l*Np]*Dr[m+k*Np];
          Srs[m+n*Np] += Dr[n+l*Np]*MM[k+l*Np]*Ds[m+k*Np];
          Ssr[m+n*Np] += Ds[n+l*Np]*MM[k+l*Np]*Dr[m+k*Np];
          Sss[m+n*Np] += Ds[n+l*Np]*MM[k+l*Np]*Ds[m+k*Np];
        }
      }
    }
  }
  SrrT = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  SrsT = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  SsrT = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  SssT = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  for (int n=0;n<Np;n++) {
    for (int m=0;m<Np;m++) {
      SrrT[m+n*Np] = Srr[n+m*Np];
      SrsT[m+n*Np] = Srs[n+m*Np];
      SsrT[m+n*Np] = Ssr[n+m*Np];
      SssT[m+n*Np] = Sss[n+m*Np];
    }
  }

  dfloat *ST = (dfloat*) calloc(3*Np*Np, sizeof(dfloat));
  for(int n=0;n<Np;++n){
    for(int m=0;m<Np;++m){
      ST[n+m*Np+0*Np*Np] = Srr[n*Np+m];
      ST[n+m*Np+1*Np*Np] = Srs[n*Np+m]+Ssr[n*Np+m];
      ST[n+m*Np+2*Np*Np] = Sss[n*Np+m];
    }
  }

  intx = (dfloat*) calloc(Nelements*Nfaces*intNfp, sizeof(dfloat));
  inty = (dfloat*) calloc(Nelements*Nfaces*intNfp, sizeof(dfloat));
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      for(int n=0;n<intNfp;++n){
        dfloat ix = 0, iy = 0;
        for(int m=0;m<Nfp;++m){
          dlong vid = vmapM[m+f*Nfp+e*Nfp*Nfaces];
          dfloat xm = x[vid];
          dfloat ym = y[vid];

          dfloat Inm = intInterp[m+n*Nfp+f*intNfp*Nfp];
          ix += Inm*xm;
          iy += Inm*ym;
        }
        dlong id = n + f*intNfp + e*Nfaces*intNfp;
        intx[id] = ix;
        inty[id] = iy;
      }
    }
  }

  o_Dr = device.malloc(Np*Np*sizeof(dfloat),
                                   Dr);

  o_Ds = device.malloc(Np*Np*sizeof(dfloat),
                                   Ds);

  o_DrT = device.malloc(Np*Np*sizeof(dfloat),
                                    DrT);

  o_DsT = device.malloc(Np*Np*sizeof(dfloat),
                                    DsT);

  o_DtT = device.malloc(Np*Np*sizeof(dfloat),
                                    DsT); // note: dummy allocated with DsT

  o_Dmatrices = device.malloc(2*Np*Np*sizeof(dfloat), DrsT);

  o_MM = device.malloc(Np*Np*sizeof(dfloat), MM);

  o_sMT = device.malloc(Np*Nfaces*Nfp*sizeof(dfloat), sMT);

  o_LIFT =
    device.malloc(Np*Nfaces*Nfp*sizeof(dfloat),
                        LIFT);

  o_LIFTT =
    device.malloc(Np*Nfaces*Nfp*sizeof(dfloat),
                        LIFTT);

  o_vgeo =
    device.malloc((Nelements+totalHaloPairs)*Nvgeo*sizeof(dfloat),
                        vgeo);

  o_cubvgeo =   device.malloc(sizeof(dfloat));// dummy

  o_sgeo =
    device.malloc(Nelements*Nfaces*Nsgeo*sizeof(dfloat),
                        sgeo);

  o_ggeo =
    device.malloc(Nelements*Nggeo*sizeof(dfloat),
                        ggeo);

  o_cubsgeo = o_sgeo; //dummy cubature geo factors

  o_SrrT = device.malloc(Np*Np*sizeof(dfloat), SrrT);
  o_SrsT = device.malloc(Np*Np*sizeof(dfloat), SrsT);
  o_SsrT = device.malloc(Np*Np*sizeof(dfloat), SsrT);
  o_SssT = device.malloc(Np*Np*sizeof(dfloat), SssT);
  o_Smatrices = device.malloc(3*Np*Np*sizeof(dfloat), ST);

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
                        cubDsWT); // dummy to align with 3d

  o_cubDWmatrices = device.malloc(2*cubNp*Np*sizeof(dfloat), cubDrsWT);

  o_intInterpT =
    device.malloc(Nfp*Nfaces*intNfp*sizeof(dfloat),
                        intInterpT);

  o_intLIFTT =
    device.malloc(Np*Nfaces*intNfp*sizeof(dfloat),
                        intLIFTT);

  o_intx =
    device.malloc(Nelements*Nfaces*intNfp*sizeof(dfloat),
                        intx);

  o_inty =
    device.malloc(Nelements*Nfaces*intNfp*sizeof(dfloat),
                        inty);

  o_intz =
    device.malloc(Nelements*Nfaces*intNfp*sizeof(dfloat),
                        inty); // dummy to align with 3d

  free(DrsT); free(ST);
}
