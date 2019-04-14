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
#include "mesh2D.hpp"

/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void meshQuad2D::SurfaceGeometricFactors(){

  /* unified storage array for geometric factors */
  Nsgeo = 7;
  sgeo = (dfloat*) calloc((Nelements+totalHaloPairs)*
                                Nsgeo*Nfp*Nfaces,
                                sizeof(dfloat));

  cubsgeo = (dfloat*) calloc((Nelements+totalHaloPairs)*
                                Nsgeo*cubNq*Nfaces,
                                sizeof(dfloat));

  for(dlong e=0;e<Nelements+totalHaloPairs;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*Nverts;

    dfloat *xe = EX + id;
    dfloat *ye = EY + id;

    for(int f=0;f<Nfaces;++f){ // for each face

      for(int i=0;i<Nfp;++i){  // for each node on face

        /* volume index of face node */
        int n = faceNodes[f*Nfp+i];

        /* local node coordinates */
        dfloat rn = r[n];
        dfloat sn = s[n];

        /* Jacobian matrix */
        dfloat xr = 0.25*( (1-sn)*(xe[1]-xe[0]) + (1+sn)*(xe[2]-xe[3]) );
        dfloat xs = 0.25*( (1-rn)*(xe[3]-xe[0]) + (1+rn)*(xe[2]-xe[1]) );
        dfloat yr = 0.25*( (1-sn)*(ye[1]-ye[0]) + (1+sn)*(ye[2]-ye[3]) );
        dfloat ys = 0.25*( (1-rn)*(ye[3]-ye[0]) + (1+rn)*(ye[2]-ye[1]) );

        /* compute geometric factors for affine coordinate transform*/
        dfloat J = xr*ys - xs*yr;

        /* face f normal and length */
        dfloat nx =   ye[(f+1)%Nverts]-ye[f];
        dfloat ny = -(xe[(f+1)%Nverts]-xe[f]);
        dfloat  d = norm2(nx,ny);

        /* output index */
        dlong base = Nsgeo*(Nfaces*Nfp*e + Nfp*f + i);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        sgeo[base+NXID] = nx/d;
        sgeo[base+NYID] = ny/d;
        sgeo[base+SJID] = d/2.;
        sgeo[base+IJID] = 1./J;

        sgeo[base+WIJID] = 1./(J*gllw[0]);
        sgeo[base+WSJID] = (d/2.)*gllw[i];
      }

      //geometric data for quadrature
      for(int i=0;i<cubNq;++i){  // for each quadrature node on face

        dfloat rn = 0., sn = 0.;

        /* interpolate local node coordinates */
        for (int j=0;j<Nfp;j++) {
          /* volume index of face node */
          int n = faceNodes[f*Nfp+j];

          rn += cubInterp[i*Nfp+j]*r[n];
          sn += cubInterp[i*Nfp+j]*s[n];
        }

        /* Jacobian matrix */
        dfloat xr = 0.25*( (1-sn)*(xe[1]-xe[0]) + (1+sn)*(xe[2]-xe[3]) );
        dfloat xs = 0.25*( (1-rn)*(xe[3]-xe[0]) + (1+rn)*(xe[2]-xe[1]) );
        dfloat yr = 0.25*( (1-sn)*(ye[1]-ye[0]) + (1+sn)*(ye[2]-ye[3]) );
        dfloat ys = 0.25*( (1-rn)*(ye[3]-ye[0]) + (1+rn)*(ye[2]-ye[1]) );

        /* compute geometric factors for affine coordinate transform*/
        dfloat J = xr*ys - xs*yr;

        /* face f normal and length */
        dfloat nx =   ye[(f+1)%Nverts]-ye[f];
        dfloat ny = -(xe[(f+1)%Nverts]-xe[f]);
        dfloat  d = norm2(nx,ny);

        /* output index */
        dlong base = Nsgeo*(Nfaces*cubNq*e + cubNq*f + i);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        cubsgeo[base+NXID] = nx/d;
        cubsgeo[base+NYID] = ny/d;
        cubsgeo[base+SJID] = d/2.;
        cubsgeo[base+IJID] = 1./J;

        cubsgeo[base+WIJID] = 1./(J*cubw[0]);
        cubsgeo[base+WSJID] = (d/2.)*cubw[i];
      }
    }
  }

  for(dlong e=0;e<Nelements;++e){ /* for each non-halo element */
    for(int n=0;n<Nfp*Nfaces;++n){
      dlong baseM = e*Nfp*Nfaces + n;
      dlong baseP = mapP[baseM];
      if(baseP<0) baseP = baseM;

      // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)
      dfloat hinvM = sgeo[baseM*Nsgeo + SJID]*sgeo[baseM*Nsgeo + IJID];
      dfloat hinvP = sgeo[baseP*Nsgeo + SJID]*sgeo[baseP*Nsgeo + IJID];

      sgeo[baseM*Nsgeo+IHID] = mymax(hinvM,hinvP);
      sgeo[baseP*Nsgeo+IHID] = mymax(hinvM,hinvP);
    }
  }
}
