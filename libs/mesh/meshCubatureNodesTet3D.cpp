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

void mesh_t::CubaturePhysicalNodesTet3D(){

  if(cubNp){
    cubx.malloc(Nelements*cubNp);
    cuby.malloc(Nelements*cubNp);
    cubz.malloc(Nelements*cubNp);

    dlong cnt = 0;
    for(dlong e=0;e<Nelements;++e){ /* for each element */

      dlong id = e*Nverts+0;

      dfloat xe1 = EX[id+0]; /* x-coordinates of vertices */
      dfloat xe2 = EX[id+1];
      dfloat xe3 = EX[id+2];
      dfloat xe4 = EX[id+3];

      dfloat ye1 = EY[id+0]; /* y-coordinates of vertices */
      dfloat ye2 = EY[id+1];
      dfloat ye3 = EY[id+2];
      dfloat ye4 = EY[id+3];

      dfloat ze1 = EZ[id+0]; /* z-coordinates of vertices */
      dfloat ze2 = EZ[id+1];
      dfloat ze3 = EZ[id+2];
      dfloat ze4 = EZ[id+3];

      for(int n=0;n<cubNp;++n){ /* for each node */

        /* (r,s,t) coordinates of interpolation nodes*/
        dfloat rn = cubr[n];
        dfloat sn = cubs[n];
        dfloat tn = cubt[n];

        /* physical coordinate of interpolation node */
        cubx[cnt] = -0.5*(1+rn+sn+tn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3 + 0.5*(1+tn)*xe4;
        cuby[cnt] = -0.5*(1+rn+sn+tn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3 + 0.5*(1+tn)*ye4;
        cubz[cnt] = -0.5*(1+rn+sn+tn)*ze1 + 0.5*(1+rn)*ze2 + 0.5*(1+sn)*ze3 + 0.5*(1+tn)*ze4;
        ++cnt;
      }
    }

    o_cubx = platform.malloc<dfloat>(Nelements*cubNp, cubx);
    o_cuby = platform.malloc<dfloat>(Nelements*cubNp, cuby);
    o_cubz = platform.malloc<dfloat>(Nelements*cubNp, cubz);
  }

  //Face cubature
  if(intNfp){
    // printf("Integration number of points: %d \n",intNfp);
    intx.malloc(Nelements*Nfaces*intNfp);
    inty.malloc(Nelements*Nfaces*intNfp);
    intz.malloc(Nelements*Nfaces*intNfp);

    for(dlong e=0;e<Nelements;++e){
      for(int f=0;f<Nfaces;++f){
        for(int n=0;n<intNfp;++n){
          dfloat ix = 0, iy = 0, iz=0;
          for(int m=0;m<Nfp;++m){
            dlong vid = vmapM[m+f*Nfp+e*Nfp*Nfaces];
            dfloat xm = x[vid];
            dfloat ym = y[vid];
            dfloat zm = z[vid];
            dfloat Inm = intInterp[m+n*Nfp+f*intNfp*Nfp];
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

    o_intx = platform.malloc<dfloat>(Nelements*Nfaces*intNfp, intx);
    o_inty = platform.malloc<dfloat>(Nelements*Nfaces*intNfp, inty);
    o_intz = platform.malloc<dfloat>(Nelements*Nfaces*intNfp, intz);
  }
}

} //namespace libp
