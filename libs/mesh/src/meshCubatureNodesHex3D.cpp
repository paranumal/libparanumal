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

void meshHex3D::CubatureNodes(){

  cubx = (dfloat*) calloc(Nelements*cubNp,sizeof(dfloat));
  cuby = (dfloat*) calloc(Nelements*cubNp,sizeof(dfloat));
  cubz = (dfloat*) calloc(Nelements*cubNp,sizeof(dfloat));

  //temp arrays
  dfloat *Ix1 = (dfloat*) calloc(Nq*Nq*cubNq, sizeof(dfloat));
  dfloat *Iy1 = (dfloat*) calloc(Nq*Nq*cubNq, sizeof(dfloat));
  dfloat *Iz1 = (dfloat*) calloc(Nq*Nq*cubNq, sizeof(dfloat));

  dfloat *Ix2 = (dfloat*) calloc(Nq*cubNq*cubNq, sizeof(dfloat));
  dfloat *Iy2 = (dfloat*) calloc(Nq*cubNq*cubNq, sizeof(dfloat));
  dfloat *Iz2 = (dfloat*) calloc(Nq*cubNq*cubNq, sizeof(dfloat));

  for(dlong e=0;e<Nelements;++e){ /* for each element */

    dfloat *xe = x + e*Np;
    dfloat *ye = y + e*Np;
    dfloat *ze = z + e*Np;
    dfloat *cubxe = cubx + e*cubNp;
    dfloat *cubye = cuby + e*cubNp;
    dfloat *cubze = cubz + e*cubNp;

    //interpolate physical coordinates to cubature
    for(int k=0;k<Nq;++k){
      for(int j=0;j<Nq;++j){
        for(int i=0;i<cubNq;++i){
          Ix1[k*Nq*cubNq+j*cubNq+i] = 0;
          Iy1[k*Nq*cubNq+j*cubNq+i] = 0;
          Iz1[k*Nq*cubNq+j*cubNq+i] = 0;
          for(int n=0;n<Nq;++n){
            Ix1[k*Nq*cubNq+j*cubNq+i] += cubInterp[i*Nq + n]*xe[k*Nq*Nq+j*Nq+n];
            Iy1[k*Nq*cubNq+j*cubNq+i] += cubInterp[i*Nq + n]*ye[k*Nq*Nq+j*Nq+n];
            Iz1[k*Nq*cubNq+j*cubNq+i] += cubInterp[i*Nq + n]*ze[k*Nq*Nq+j*Nq+n];
          }
        }
      }
    }

    for(int k=0;k<Nq;++k){
      for(int j=0;j<cubNq;++j){
        for(int i=0;i<cubNq;++i){
          Ix2[k*cubNq*cubNq+j*cubNq+i] = 0;
          Iy2[k*cubNq*cubNq+j*cubNq+i] = 0;
          Iz2[k*cubNq*cubNq+j*cubNq+i] = 0;
          for(int n=0;n<Nq;++n){
            Ix2[k*cubNq*cubNq+j*cubNq+i] += cubInterp[j*Nq + n]*Ix1[k*Nq*cubNq+n*cubNq+i];
            Iy2[k*cubNq*cubNq+j*cubNq+i] += cubInterp[j*Nq + n]*Iy1[k*Nq*cubNq+n*cubNq+i];
            Iz2[k*cubNq*cubNq+j*cubNq+i] += cubInterp[j*Nq + n]*Iz1[k*Nq*cubNq+n*cubNq+i];
          }
        }
      }
    }

    for(int k=0;k<cubNq;++k){
      for(int j=0;j<cubNq;++j){
        for(int i=0;i<cubNq;++i){
          cubxe[k*cubNq*cubNq+j*cubNq+i] = 0;
          cubye[k*cubNq*cubNq+j*cubNq+i] = 0;
          cubze[k*cubNq*cubNq+j*cubNq+i] = 0;
          for(int n=0;n<Nq;++n){
            cubxe[k*cubNq*cubNq+j*cubNq+i] += cubInterp[k*Nq + n]*Ix2[n*cubNq*cubNq+j*cubNq+i];
            cubye[k*cubNq*cubNq+j*cubNq+i] += cubInterp[k*Nq + n]*Iy2[n*cubNq*cubNq+j*cubNq+i];
            cubze[k*cubNq*cubNq+j*cubNq+i] += cubInterp[k*Nq + n]*Iz2[n*cubNq*cubNq+j*cubNq+i];
          }
        }
      }
    }
  }

  free(Ix1); free(Iy1); free(Iz1);
  free(Ix2); free(Iy2); free(Iz2);

  o_cubx = device.malloc(Nelements*cubNp*sizeof(dfloat), cubx);
  o_cuby = device.malloc(Nelements*cubNp*sizeof(dfloat), cuby);
  o_cubz = device.malloc(Nelements*cubNp*sizeof(dfloat), cubz);

  //Face cubature
  intx = (dfloat*) calloc(Nelements*Nfaces*cubNfp, sizeof(dfloat));
  inty = (dfloat*) calloc(Nelements*Nfaces*cubNfp, sizeof(dfloat));
  intz = (dfloat*) calloc(Nelements*Nfaces*cubNfp, sizeof(dfloat));

  dfloat *ix = (dfloat *) calloc(cubNq*Nq,sizeof(dfloat));
  dfloat *iy = (dfloat *) calloc(cubNq*Nq,sizeof(dfloat));
  dfloat *iz = (dfloat *) calloc(cubNq*Nq,sizeof(dfloat));
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      //interpolate in i
      for(int ny=0;ny<Nq;++ny){
        for(int nx=0;nx<cubNq;++nx){
          ix[nx+cubNq*ny] = 0;
          iy[nx+cubNq*ny] = 0;
          iz[nx+cubNq*ny] = 0;

          for(int m=0;m<Nq;++m){
            dlong vid = m+ny*Nq+f*Nfp+e*Nfp*Nfaces;
            dlong idM = vmapM[vid];

            dfloat xm = x[idM];
            dfloat ym = y[idM];
            dfloat zm = z[idM];

            dfloat Inm = cubInterp[m+nx*Nq];
            ix[nx+cubNq*ny] += Inm*xm;
            iy[nx+cubNq*ny] += Inm*ym;
            iz[nx+cubNq*ny] += Inm*zm;
          }
        }
      }

      //interpolate in j and store
      for(int ny=0;ny<cubNq;++ny){
        for(int nx=0;nx<cubNq;++nx){
          dfloat xn=0.0, yn=0.0, zn=0.0;

          for(int m=0;m<Nq;++m){
            dfloat xm = ix[nx + m*cubNq];
            dfloat ym = iy[nx + m*cubNq];
            dfloat zm = iz[nx + m*cubNq];

            dfloat Inm = cubInterp[m+ny*Nq];
            xn += Inm*xm;
            yn += Inm*ym;
            zn += Inm*zm;
          }

          dlong id = nx + ny*cubNq + f*cubNfp + e*Nfaces*cubNfp;
          intx[id] = xn;
          inty[id] = yn;
          intz[id] = zn;
        }
      }
    }
  }
  free(ix); free(iy); free(iz);

  o_intx = device.malloc(Nelements*Nfaces*cubNfp*sizeof(dfloat), intx);
  o_inty = device.malloc(Nelements*Nfaces*cubNfp*sizeof(dfloat), inty);
  o_intz = device.malloc(Nelements*Nfaces*cubNfp*sizeof(dfloat), intz);
}
