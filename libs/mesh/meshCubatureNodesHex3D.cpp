/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

namespace libp {

void mesh_t::CubaturePhysicalNodesHex3D(){

  cubx.malloc(Nelements*cubNp);
  cuby.malloc(Nelements*cubNp);
  cubz.malloc(Nelements*cubNp);

  //temp arrays
  memory<dfloat> Ix1(Nq*Nq*cubNq);
  memory<dfloat> Iy1(Nq*Nq*cubNq);
  memory<dfloat> Iz1(Nq*Nq*cubNq);

  memory<dfloat> Ix2(Nq*cubNq*cubNq);
  memory<dfloat> Iy2(Nq*cubNq*cubNq);
  memory<dfloat> Iz2(Nq*cubNq*cubNq);

  for(dlong e=0;e<Nelements;++e){ /* for each element */

    dfloat *xe = x.ptr() + e*Np;
    dfloat *ye = y.ptr() + e*Np;
    dfloat *ze = z.ptr() + e*Np;
    dfloat *cubxe = cubx.ptr() + e*cubNp;
    dfloat *cubye = cuby.ptr() + e*cubNp;
    dfloat *cubze = cubz.ptr() + e*cubNp;

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

  o_cubx = platform.malloc<dfloat>(Nelements*cubNp, cubx);
  o_cuby = platform.malloc<dfloat>(Nelements*cubNp, cuby);
  o_cubz = platform.malloc<dfloat>(Nelements*cubNp, cubz);

  //Face cubature
  intx.malloc(Nelements*Nfaces*cubNfp);
  inty.malloc(Nelements*Nfaces*cubNfp);
  intz.malloc(Nelements*Nfaces*cubNfp);

  memory<dfloat> ix(cubNq*Nq);
  memory<dfloat> iy(cubNq*Nq);
  memory<dfloat> iz(cubNq*Nq);
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

  o_intx = platform.malloc<dfloat>(Nelements*Nfaces*cubNfp, intx);
  o_inty = platform.malloc<dfloat>(Nelements*Nfaces*cubNfp, inty);
  o_intz = platform.malloc<dfloat>(Nelements*Nfaces*cubNfp, intz);
}

} //namespace libp
