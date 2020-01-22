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

void meshHex3D::CubatureSetup(){

  /* Quadrature data */
  cubN = N+1;
  cubNq = cubN+1;
  cubNp = cubNq*cubNq*cubNq;
  cubNfp = cubNq*cubNq;
  intNfp = cubNq*cubNq;

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

  o_cubD = device.malloc(cubNq*cubNq*sizeof(dfloat), cubD);

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

  cubsgeo = (dfloat*) calloc(Nelements*Nsgeo*cubNq*cubNq*Nfaces,
                                sizeof(dfloat));

  //temp arrays
  dfloat *xre = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *xse = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *xte = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *yre = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *yse = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *yte = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *zre = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *zse = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *zte = (dfloat*) calloc(Np, sizeof(dfloat));

  dfloat *xre1 = (dfloat*) calloc(Nq*Nq*cubNq, sizeof(dfloat));
  dfloat *xse1 = (dfloat*) calloc(Nq*Nq*cubNq, sizeof(dfloat));
  dfloat *xte1 = (dfloat*) calloc(Nq*Nq*cubNq, sizeof(dfloat));
  dfloat *yre1 = (dfloat*) calloc(Nq*Nq*cubNq, sizeof(dfloat));
  dfloat *yse1 = (dfloat*) calloc(Nq*Nq*cubNq, sizeof(dfloat));
  dfloat *yte1 = (dfloat*) calloc(Nq*Nq*cubNq, sizeof(dfloat));
  dfloat *zre1 = (dfloat*) calloc(Nq*Nq*cubNq, sizeof(dfloat));
  dfloat *zse1 = (dfloat*) calloc(Nq*Nq*cubNq, sizeof(dfloat));
  dfloat *zte1 = (dfloat*) calloc(Nq*Nq*cubNq, sizeof(dfloat));

  dfloat *xre2 = (dfloat*) calloc(Nq*cubNq*cubNq, sizeof(dfloat));
  dfloat *xse2 = (dfloat*) calloc(Nq*cubNq*cubNq, sizeof(dfloat));
  dfloat *xte2 = (dfloat*) calloc(Nq*cubNq*cubNq, sizeof(dfloat));
  dfloat *yre2 = (dfloat*) calloc(Nq*cubNq*cubNq, sizeof(dfloat));
  dfloat *yse2 = (dfloat*) calloc(Nq*cubNq*cubNq, sizeof(dfloat));
  dfloat *yte2 = (dfloat*) calloc(Nq*cubNq*cubNq, sizeof(dfloat));
  dfloat *zre2 = (dfloat*) calloc(Nq*cubNq*cubNq, sizeof(dfloat));
  dfloat *zse2 = (dfloat*) calloc(Nq*cubNq*cubNq, sizeof(dfloat));
  dfloat *zte2 = (dfloat*) calloc(Nq*cubNq*cubNq, sizeof(dfloat));

  //surface temp arrays
  dfloat *xr1 = (dfloat*) calloc(Nq*cubNq, sizeof(dfloat));
  dfloat *xs1 = (dfloat*) calloc(Nq*cubNq, sizeof(dfloat));
  dfloat *xt1 = (dfloat*) calloc(Nq*cubNq, sizeof(dfloat));
  dfloat *yr1 = (dfloat*) calloc(Nq*cubNq, sizeof(dfloat));
  dfloat *ys1 = (dfloat*) calloc(Nq*cubNq, sizeof(dfloat));
  dfloat *yt1 = (dfloat*) calloc(Nq*cubNq, sizeof(dfloat));
  dfloat *zr1 = (dfloat*) calloc(Nq*cubNq, sizeof(dfloat));
  dfloat *zs1 = (dfloat*) calloc(Nq*cubNq, sizeof(dfloat));
  dfloat *zt1 = (dfloat*) calloc(Nq*cubNq, sizeof(dfloat));

  //geometric data for quadrature
  for(dlong e=0;e<Nelements;++e){ /* for each element */
    for(int k=0;k<Nq;++k){
      for(int j=0;j<Nq;++j){
        for(int i=0;i<Nq;++i){
          int n = i + j*Nq + k*Nq*Nq;

          //differentiate physical coordinates
          xre[n] = 0; xse[n] = 0; xte[n] = 0;
          yre[n] = 0; yse[n] = 0; yte[n] = 0;
          zre[n] = 0; zse[n] = 0; zte[n] = 0;

          for(int m=0;m<Nq;++m){
            int idr = e*Np + k*Nq*Nq + j*Nq + m;
            int ids = e*Np + k*Nq*Nq + m*Nq + i;
            int idt = e*Np + m*Nq*Nq + j*Nq + i;
            xre[n] += D[i*Nq+m]*x[idr];
            xse[n] += D[j*Nq+m]*x[ids];
            xte[n] += D[k*Nq+m]*x[idt];
            yre[n] += D[i*Nq+m]*y[idr];
            yse[n] += D[j*Nq+m]*y[ids];
            yte[n] += D[k*Nq+m]*y[idt];
            zre[n] += D[i*Nq+m]*z[idr];
            zse[n] += D[j*Nq+m]*z[ids];
            zte[n] += D[k*Nq+m]*z[idt];
          }
        }
      }
    }

    //interpolate derivaties to cubature
    for(int k=0;k<Nq;++k){
      for(int j=0;j<Nq;++j){
        for(int i=0;i<cubNq;++i){
          dlong id = k*Nq*cubNq+j*cubNq+i;
          xre1[id] = 0.0; xse1[id] = 0.0;  xte1[id] = 0.0;
          yre1[id] = 0.0; yse1[id] = 0.0;  yte1[id] = 0.0;
          zre1[id] = 0.0; zse1[id] = 0.0;  zte1[id] = 0.0;
          for(int n=0;n<Nq;++n){
            dlong idn = k*Nq*Nq+j*Nq+n;
            xre1[id] += cubInterp[i*Nq + n]*xre[idn];
            xse1[id] += cubInterp[i*Nq + n]*xse[idn];
            xte1[id] += cubInterp[i*Nq + n]*xte[idn];
            yre1[id] += cubInterp[i*Nq + n]*yre[idn];
            yse1[id] += cubInterp[i*Nq + n]*yse[idn];
            yte1[id] += cubInterp[i*Nq + n]*yte[idn];
            zre1[id] += cubInterp[i*Nq + n]*zre[idn];
            zse1[id] += cubInterp[i*Nq + n]*zse[idn];
            zte1[id] += cubInterp[i*Nq + n]*zte[idn];
          }
        }
      }
    }

    for(int k=0;k<Nq;++k){
      for(int j=0;j<cubNq;++j){
        for(int i=0;i<cubNq;++i){
          dlong id = k*cubNq*cubNq+j*cubNq+i;
          xre2[id] = 0.0; xse2[id] = 0.0;  xte2[id] = 0.0;
          yre2[id] = 0.0; yse2[id] = 0.0;  yte2[id] = 0.0;
          zre2[id] = 0.0; zse2[id] = 0.0;  zte2[id] = 0.0;
          for(int n=0;n<Nq;++n){
            dlong idn = k*Nq*cubNq+n*cubNq+i;
            xre2[id] += cubInterp[j*Nq + n]*xre1[idn];
            xse2[id] += cubInterp[j*Nq + n]*xse1[idn];
            xte2[id] += cubInterp[j*Nq + n]*xte1[idn];
            yre2[id] += cubInterp[j*Nq + n]*yre1[idn];
            yse2[id] += cubInterp[j*Nq + n]*yse1[idn];
            yte2[id] += cubInterp[j*Nq + n]*yte1[idn];
            zre2[id] += cubInterp[j*Nq + n]*zre1[idn];
            zse2[id] += cubInterp[j*Nq + n]*zse1[idn];
            zte2[id] += cubInterp[j*Nq + n]*zte1[idn];
          }
        }
      }
    }

    for(int k=0;k<cubNq;++k){
      for(int j=0;j<cubNq;++j){
        for(int i=0;i<cubNq;++i){
          dfloat xr = 0.0, xs = 0.0, xt = 0.0;
          dfloat yr = 0.0, ys = 0.0, yt = 0.0;
          dfloat zr = 0.0, zs = 0.0, zt = 0.0;
          for(int n=0;n<Nq;++n){
            dlong idn = n*cubNq*cubNq+j*cubNq+i;
            xr += cubInterp[k*Nq + n]*xre2[idn];
            xs += cubInterp[k*Nq + n]*xse2[idn];
            xt += cubInterp[k*Nq + n]*xte2[idn];
            yr += cubInterp[k*Nq + n]*yre2[idn];
            ys += cubInterp[k*Nq + n]*yse2[idn];
            yt += cubInterp[k*Nq + n]*yte2[idn];
            zr += cubInterp[k*Nq + n]*zre2[idn];
            zs += cubInterp[k*Nq + n]*zse2[idn];
            zt += cubInterp[k*Nq + n]*zte2[idn];
          }

          /* compute geometric factors for affine coordinate transform*/
          dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);

          if(J<1e-8) {
            stringstream ss;
            ss << "Negative J found at element " << e << "\n";
            LIBP_ABORT(ss.str())
          }

          dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
          dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
          dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;

          dfloat JW = J*cubw[i]*cubw[j]*cubw[k];

          /* store geometric factors */
          dlong base = Nvgeo*cubNp*e + i + j*cubNq + k*cubNq*cubNq;
          cubvgeo[base + cubNp*RXID] = rx;
          cubvgeo[base + cubNp*RYID] = ry;
          cubvgeo[base + cubNp*RZID] = rz;

          cubvgeo[base + cubNp*SXID] = sx;
          cubvgeo[base + cubNp*SYID] = sy;
          cubvgeo[base + cubNp*SZID] = sz;

          cubvgeo[base + cubNp*TXID] = tx;
          cubvgeo[base + cubNp*TYID] = ty;
          cubvgeo[base + cubNp*TZID] = tz;

          cubvgeo[base + cubNp*JID]  = J;
          cubvgeo[base + cubNp*JWID] = JW;
          cubvgeo[base + cubNp*IJWID] = 1./JW;

          /* store second order geometric factors */
          base = Nggeo*cubNp*e + i + j*cubNq + k*cubNq*cubNq;
          cubggeo[base + cubNp*G00ID] = JW*(rx*rx + ry*ry + rz*rz);
          cubggeo[base + cubNp*G01ID] = JW*(rx*sx + ry*sy + rz*sz);
          cubggeo[base + cubNp*G02ID] = JW*(rx*tx + ry*ty + rz*tz);
          cubggeo[base + cubNp*G11ID] = JW*(sx*sx + sy*sy + sz*sz);
          cubggeo[base + cubNp*G12ID] = JW*(sx*tx + sy*ty + sz*tz);
          cubggeo[base + cubNp*G22ID] = JW*(tx*tx + ty*ty + tz*tz);
          cubggeo[base + cubNp*GWJID] = JW;
        }
      }
    }

    //surface geometric data for quadrature
    for(int f=0;f<Nfaces;++f){ // for each face

      for(int j=0;j<Nq;++j){
        for(int i=0;i<cubNq;++i){
          dlong id = j*cubNq+i;
          //interpolate derivatives of physical coordinates
          xr1[id] = 0.0; xs1[id] = 0.0; xt1[id] = 0.0;
          yr1[id] = 0.0; ys1[id] = 0.0; yt1[id] = 0.0;
          zr1[id] = 0.0; zs1[id] = 0.0; zt1[id] = 0.0;
          for(int n=0;n<Nq;++n){  // for each node on face
            /* volume index of face node */
            int idn = faceNodes[f*Nfp+j*Nq+n];
            xr1[id] += cubInterp[i*Nq + n]*xre[idn];
            xs1[id] += cubInterp[i*Nq + n]*xse[idn];
            xt1[id] += cubInterp[i*Nq + n]*xte[idn];
            yr1[id] += cubInterp[i*Nq + n]*yre[idn];
            ys1[id] += cubInterp[i*Nq + n]*yse[idn];
            yt1[id] += cubInterp[i*Nq + n]*yte[idn];
            zr1[id] += cubInterp[i*Nq + n]*zre[idn];
            zs1[id] += cubInterp[i*Nq + n]*zse[idn];
            zt1[id] += cubInterp[i*Nq + n]*zte[idn];
          }
        }
      }

      for(int j=0;j<cubNq;++j){
        for(int i=0;i<cubNq;++i){
          //interpolate derivatives of physical coordinates
          dfloat xr = 0.0, xs = 0.0, xt = 0.0;
          dfloat yr = 0.0, ys = 0.0, yt = 0.0;
          dfloat zr = 0.0, zs = 0.0, zt = 0.0;
          for(int n=0;n<Nq;++n){  // for each node on face
            /* volume index of face node */
            xr += cubInterp[j*Nq + n]*xr1[n*cubNq+i];
            xs += cubInterp[j*Nq + n]*xs1[n*cubNq+i];
            xt += cubInterp[j*Nq + n]*xt1[n*cubNq+i];
            yr += cubInterp[j*Nq + n]*yr1[n*cubNq+i];
            ys += cubInterp[j*Nq + n]*ys1[n*cubNq+i];
            yt += cubInterp[j*Nq + n]*yt1[n*cubNq+i];
            zr += cubInterp[j*Nq + n]*zr1[n*cubNq+i];
            zs += cubInterp[j*Nq + n]*zs1[n*cubNq+i];
            zt += cubInterp[j*Nq + n]*zt1[n*cubNq+i];
          }


          /* compute geometric factors for affine coordinate transform*/
          dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);

          dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
          dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
          dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;

          /* face f normal and length */
          dfloat nx=0.0, ny=0.0, nz=0.0;
          switch(f){
          case 0: nx = -tx; ny = -ty; nz = -tz; break;
          case 1: nx = -sx; ny = -sy; nz = -sz; break;
          case 2: nx = +rx; ny = +ry; nz = +rz; break;
          case 3: nx = +sx; ny = +sy; nz = +sz; break;
          case 4: nx = -rx; ny = -ry; nz = -rz; break;
          case 5: nx = +tx; ny = +ty; nz = +tz; break;
          }
          dfloat sJ = sqrt(nx*nx+ny*ny+nz*nz);
          nx /= sJ; ny /= sJ; nz /= sJ;
          sJ *= J;

          /* output index */
          dlong base = Nsgeo*(Nfaces*cubNq*cubNq*e + cubNq*cubNq*f + j*cubNq + i);

          /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
          cubsgeo[base+NXID] = nx;
          cubsgeo[base+NYID] = ny;
          cubsgeo[base+NZID] = nz;
          cubsgeo[base+SJID] = sJ;
          cubsgeo[base+IJID] = 1./J;

          cubsgeo[base+WIJID] = 1./(J*cubw[0]);
          cubsgeo[base+WSJID] = sJ*cubw[i]*cubw[j];
        }
      }
    }
  }


  o_cubvgeo =
    device.malloc(Nelements*Nvgeo*cubNp*sizeof(dfloat),
        cubvgeo);

  o_cubsgeo =
    device.malloc(Nelements*Nfaces*cubNq*cubNq*Nsgeo*sizeof(dfloat),
        cubsgeo);

  free(xre); free(xse); free(xte);
  free(yre); free(yse); free(yte);
  free(zre); free(zse); free(zte);
  free(xre1); free(xse1); free(xte1);
  free(yre1); free(yse1); free(yte1);
  free(zre1); free(zse1); free(zte1);
  free(xre2); free(xse2); free(xte2);
  free(yre2); free(yse2); free(yte2);
  free(zre2); free(zse2); free(zte2);

  free(xr1); free(xs1); free(xt1);
  free(yr1); free(ys1); free(yt1);
  free(zr1); free(zs1); free(zt1);
}
