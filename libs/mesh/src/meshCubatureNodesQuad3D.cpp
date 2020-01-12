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

void meshQuad3D::CubatureNodes(){

  cubx = (dfloat*) calloc(Nelements*cubNp,sizeof(dfloat));
  cuby = (dfloat*) calloc(Nelements*cubNp,sizeof(dfloat));
  cubz = (dfloat*) calloc(Nelements*cubNp,sizeof(dfloat));

  //temp arrays
  dfloat *Ix1 = (dfloat*) calloc(Nq*cubNq, sizeof(dfloat));
  dfloat *Iy1 = (dfloat*) calloc(Nq*cubNq, sizeof(dfloat));
  dfloat *Iz1 = (dfloat*) calloc(Nq*cubNq, sizeof(dfloat));

  for(dlong e=0;e<Nelements;++e){ /* for each element */

    dfloat *xe = x + e*Np;
    dfloat *ye = y + e*Np;
    dfloat *ze = z + e*Np;
    dfloat *cubxe = cubx + e*cubNp;
    dfloat *cubye = cuby + e*cubNp;
    dfloat *cubze = cubz + e*cubNp;

    //interpolate physical coordinates to cubature
    for(int j=0;j<Nq;++j){
      for(int i=0;i<cubNq;++i){
        Ix1[j*cubNq+i] = 0;
        Iy1[j*cubNq+i] = 0;
        Iz1[j*cubNq+i] = 0;
        for(int n=0;n<Nq;++n){
          Ix1[j*cubNq+i] += cubInterp[i*Nq + n]*xe[j*Nq+n];
          Iy1[j*cubNq+i] += cubInterp[i*Nq + n]*ye[j*Nq+n];
          Iz1[j*cubNq+i] += cubInterp[i*Nq + n]*ze[j*Nq+n];
        }
      }
    }

    for(int j=0;j<cubNq;++j){
      for(int i=0;i<cubNq;++i){
        dfloat xlin = 0.0;
        dfloat ylin = 0.0;
        dfloat zlin = 0.0;
        for(int n=0;n<Nq;++n){
          xlin += cubInterp[j*Nq + n]*Ix1[n*cubNq+i];
          ylin += cubInterp[j*Nq + n]*Iy1[n*cubNq+i];
          zlin += cubInterp[j*Nq + n]*Iz1[n*cubNq+i];
        }

        // project to sphere
        dfloat rlin = sqrt(xlin*xlin+ylin*ylin+zlin*zlin);
        cubxe[j*cubNq+i] = xlin/rlin;
        cubye[j*cubNq+i] = ylin/rlin;
        cubze[j*cubNq+i] = zlin/rlin;
      }
    }
  }

  free(Ix1);
  free(Iy1);

  o_cubx = device.malloc(Nelements*cubNp*sizeof(dfloat), cubx);
  o_cuby = device.malloc(Nelements*cubNp*sizeof(dfloat), cuby);
  o_cubz = device.malloc(Nelements*cubNp*sizeof(dfloat), cubz);

  //Face cubature
  intx = (dfloat*) calloc(Nelements*Nfaces*cubNq, sizeof(dfloat));
  inty = (dfloat*) calloc(Nelements*Nfaces*cubNq, sizeof(dfloat));
  intz = (dfloat*) calloc(Nelements*Nfaces*cubNq, sizeof(dfloat));
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      for(int n=0;n<cubNq;++n){
        dfloat ix = 0, iy = 0, iz = 0;
        for(int m=0;m<Nq;++m){
          dlong vid = vmapM[m+f*Nfp+e*Nfp*Nfaces];
          dfloat xm = x[vid];
          dfloat ym = y[vid];
          dfloat zm = z[vid];

          dfloat Inm = cubInterp[m+n*Nq];
          ix += Inm*xm;
          iy += Inm*ym;
          iz += Inm*zm;
        }
        dlong id = n + f*cubNq + e*Nfaces*cubNq;
        // project to sphere
        dfloat rlin = sqrt(ix*ix+iy*iy+iz*iz);
        intx[id] = ix/rlin;
        inty[id] = iy/rlin;
        intz[id] = iz/rlin;
      }
    }
  }

  o_intx = device.malloc(Nelements*Nfaces*cubNq*sizeof(dfloat), intx);
  o_inty = device.malloc(Nelements*Nfaces*cubNq*sizeof(dfloat), inty);
  o_intz = device.malloc(Nelements*Nfaces*cubNq*sizeof(dfloat), intz);
}
