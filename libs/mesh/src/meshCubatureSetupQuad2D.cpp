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
#include "mesh3D.hpp"

void meshQuad3D::CubatureSetup(){
  mesh_t *mesh_p = (mesh_t*) this;
  meshQuad2D* trimesh = (meshQuad2D*) mesh_p;
  trimesh->meshQuad2D::CubatureSetup();
}

void meshQuad2D::CubatureSetup(){

  /* Quadrature data */
  cubN = N+1;
  cubNq = cubN+1;
  cubNp = cubNq*cubNq;
  cubNfp = cubNq;
  intNfp = cubNq;

  // cubN+1 point Gauss-Legendre quadrature
  cubr = (dfloat *) malloc(cubNq*sizeof(dfloat));
  cubw = (dfloat *) malloc(cubNq*sizeof(dfloat));
  JacobiGQ(0, 0, cubN, cubr, cubw);

  // GLL to GL interpolation matrix
  cubInterp = (dfloat *) malloc(Nq*cubNq*sizeof(dfloat));
  InterpolationMatrix1D(N, Nq, gllz, cubNq, cubr, cubInterp);

  cubDW = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
  cubProject = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
  CubatureWeakDmatrices1D(N, Nq, gllz, cubNq, cubr, cubw,
                             cubDW, cubProject);

  // GLL to GL differentiation matrix
  cubD = (dfloat *) malloc(Nq*cubNq*sizeof(dfloat));
  Dmatrix1D(N, Nq, gllz, cubNq, cubr, cubD);

  // GL to GL differentiation matrix
  gjD = (dfloat *) malloc(cubNq*cubNq*sizeof(dfloat));
  Dmatrix1D(cubN, cubNq, cubr, cubNq, cubr, gjD);

  // add compile time constants to kernels
  props["defines/" "p_cubNq"]= cubNq;
  props["defines/" "p_cubNp"]= cubNp;
  props["defines/" "p_intNfp"]= intNfp;
  props["defines/" "p_intNfpNfaces"]= intNfp*Nfaces;
  props["defines/" "p_cubNfp"]= cubNfp;

  dfloat *cubDWT = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *cubProjectT = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *cubInterpT = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  for(int n=0;n<Nq;++n){
    for(int m=0;m<cubNq;++m){
      cubDWT[n+m*Nq] = cubDW[n*cubNq+m];
      cubProjectT[n+m*Nq] = cubProject[n*cubNq+m];
      cubInterpT[m+n*cubNq] = cubInterp[m*Nq+n];
    }
  }

  o_cubInterpT =
    device.malloc(Nq*cubNq*sizeof(dfloat),
        cubInterpT);

  o_cubProjectT =
    device.malloc(Nq*cubNq*sizeof(dfloat),
        cubProjectT);

  o_cubDWT =
    device.malloc(Nq*cubNq*sizeof(dfloat),
        cubDWT);

  o_cubDWmatrices = device.malloc(cubNq*Nq*sizeof(dfloat), cubDWT);

  o_intInterpT = device.malloc(cubNq*Nq*sizeof(dfloat));
  o_intInterpT.copyFrom(o_cubInterpT);

  o_intLIFTT = device.malloc(cubNq*Nq*sizeof(dfloat));
  o_intLIFTT.copyFrom(o_cubProjectT);

  free(cubDWT);
  free(cubProjectT);
  free(cubInterpT);

  cubvgeo = (dfloat*) calloc(Nelements*Nvgeo*cubNp, sizeof(dfloat));
  cubggeo = (dfloat*) calloc(Nelements*Nggeo*cubNp, sizeof(dfloat));

  cubsgeo = (dfloat*) calloc((Nelements+totalHaloPairs)*
                                Nsgeo*cubNq*Nfaces,
                                sizeof(dfloat));

  //temp arrays
  dfloat *xre = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *xse = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *yre = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *yse = (dfloat*) calloc(Np, sizeof(dfloat));

  dfloat *xre1 = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *xse1 = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *yre1 = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *yse1 = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));

  //geometric data for quadrature
  for(dlong e=0;e<Nelements;++e){ /* for each element */
    for(int j=0;j<Nq;++j){
      for(int i=0;i<Nq;++i){
        int n = i + j*Nq;

        //differentiate physical coordinates
        xre[n] = 0.0; xse[n] = 0.0;
        yre[n] = 0.0; yse[n] = 0.0;

        for(int m=0;m<Nq;++m){
          int idr = e*Np + j*Nq + m;
          int ids = e*Np + m*Nq + i;
          xre[n] += D[i*Nq+m]*x[idr];
          xse[n] += D[j*Nq+m]*x[ids];
          yre[n] += D[i*Nq+m]*y[idr];
          yse[n] += D[j*Nq+m]*y[ids];
        }
      }
    }

    //interpolate derivaties to cubature
    for(int j=0;j<Nq;++j){
      for(int i=0;i<cubNq;++i){
        xre1[j*cubNq+i] = 0.0; xse1[j*cubNq+i] = 0.0;
        yre1[j*cubNq+i] = 0.0; yse1[j*cubNq+i] = 0.0;
        for(int n=0;n<Nq;++n){
          xre1[j*cubNq+i] += cubInterp[i*Nq + n]*xre[j*Nq+n];
          xse1[j*cubNq+i] += cubInterp[i*Nq + n]*xse[j*Nq+n];
          yre1[j*cubNq+i] += cubInterp[i*Nq + n]*yre[j*Nq+n];
          yse1[j*cubNq+i] += cubInterp[i*Nq + n]*yse[j*Nq+n];
        }
      }
    }

    for(int j=0;j<cubNq;++j){
      for(int i=0;i<cubNq;++i){
        dfloat xr = 0.0, xs = 0.0;
        dfloat yr = 0.0, ys = 0.0;
        for(int n=0;n<Nq;++n){
          xr += cubInterp[j*Nq + n]*xre1[n*cubNq+i];
          xs += cubInterp[j*Nq + n]*xse1[n*cubNq+i];
          yr += cubInterp[j*Nq + n]*yre1[n*cubNq+i];
          ys += cubInterp[j*Nq + n]*yse1[n*cubNq+i];
        }

        /* compute geometric factors for affine coordinate transform*/
        dfloat J = xr*ys - xs*yr;

        if(J<1e-8) {
          stringstream ss;
          ss << "Negative J found at element " << e << "\n";
          LIBP_ABORT(ss.str())
        }
        dfloat rx =  ys/J;
        dfloat ry = -xs/J;
        dfloat sx = -yr/J;
        dfloat sy =  xr/J;

        dfloat JW = J*cubw[i]*cubw[j];

        /* store geometric factors */
        dlong base = Nvgeo*cubNp*e + i + j*cubNq;
        cubvgeo[base + cubNp*RXID] = rx;
        cubvgeo[base + cubNp*RYID] = ry;
        cubvgeo[base + cubNp*SXID] = sx;
        cubvgeo[base + cubNp*SYID] = sy;
        cubvgeo[base + cubNp*JID]  = J;
        cubvgeo[base + cubNp*JWID] = JW;
        cubvgeo[base + cubNp*IJWID] = 1./JW;

        /* store second order geometric factors */
        base = Nggeo*cubNp*e + i + j*cubNq;
        cubggeo[base + cubNp*G00ID] = JW*(rx*rx + ry*ry);
        cubggeo[base + cubNp*G01ID] = JW*(rx*sx + ry*sy);
        cubggeo[base + cubNp*G11ID] = JW*(sx*sx + sy*sy);
        cubggeo[base + cubNp*GWJID] = JW;
      }
    }

    for(int f=0;f<Nfaces;++f){ // for each face
      for(int m=0;m<cubNq;++m){  // for each node on face

        //interpolate derivatives of physical coordinates
        dfloat xr = 0.0, xs = 0.0;
        dfloat yr = 0.0, ys = 0.0;
        for(int n=0;n<Nfp;++n){  // for each node on face
          /* volume index of face node */
          int idn = faceNodes[f*Nfp+n];
          xr += cubInterp[m*Nq + n]*xre[idn];
          xs += cubInterp[m*Nq + n]*xse[idn];
          yr += cubInterp[m*Nq + n]*yre[idn];
          ys += cubInterp[m*Nq + n]*yse[idn];
        }

        /* compute geometric factors for affine coordinate transform*/
        dfloat J = xr*ys - xs*yr;

        dfloat rx =  ys/J;
        dfloat ry = -xs/J;
        dfloat sx = -yr/J;
        dfloat sy =  xr/J;

        /* face f normal and length */
        dfloat nx=0.0, ny=0.0;
        switch(f){
        case 0: nx = -sx; ny = -sy; break;
        case 1: nx = +rx; ny = +ry; break;
        case 2: nx = +sx; ny = +sy; break;
        case 3: nx = -rx; ny = -ry; break;
        }
        dfloat  sJ = sqrt((nx)*(nx)+(ny)*(ny));
        nx /= sJ; ny /= sJ;
        sJ *= J;

        /* output index */
        dlong base = Nsgeo*(Nfaces*cubNq*e + cubNq*f + m);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        cubsgeo[base+NXID] = nx;
        cubsgeo[base+NYID] = ny;
        cubsgeo[base+SJID] = sJ;
        cubsgeo[base+IJID] = 1./J;

        cubsgeo[base+WIJID] = 1./(J*cubw[0]);
        cubsgeo[base+WSJID] = sJ*cubw[m];
      }
    }
  }

  o_cubvgeo =
    device.malloc(Nelements*Nvgeo*cubNp*sizeof(dfloat),
        cubvgeo);

  o_cubggeo =
    device.malloc(Nelements*Nggeo*cubNp*sizeof(dfloat),
        cubvgeo);

  o_cubsgeo =
    device.malloc(Nelements*Nfaces*cubNq*Nsgeo*sizeof(dfloat),
        cubsgeo);

  free(xre); free(xse); free(yre); free(yse);
  free(xre1); free(xse1); free(yre1); free(yse1);
}
