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
#include "mesh3D.hpp"

void meshQuad3D::CubatureSetup(){

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
  InterpolationMatrix1D(N, Nq, r, cubNq, cubr, cubInterp); //uses the fact that r = gllz for 1:Nq

  //cubature project cubProject = cubInterp^T
  cubProject = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  matrixTranspose(cubNq, Nq, cubInterp, Nq, cubProject, cubNq);

  //cubature derivates matrix, cubD: differentiate on cubature nodes
  cubD = (dfloat *) malloc(cubNq*cubNq*sizeof(dfloat));
  Dmatrix1D(cubN, cubNq, cubr, cubNq, cubr, cubD);

  // weak cubature derivative cubPDT = cubProject * cubD^T
  cubPDT  = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  CubatureWeakDmatrix1D(Nq, cubNq, cubProject, cubD, cubPDT);

  // add compile time constants to kernels
  props["defines/" "p_cubNq"]= cubNq;
  props["defines/" "p_cubNp"]= cubNp;
  props["defines/" "p_intNfp"]= intNfp;
  props["defines/" "p_intNfpNfaces"]= intNfp*Nfaces;
  props["defines/" "p_cubNfp"]= cubNfp;

  // build transposes (we hold matrices as column major on device)
  dfloat *cubProjectT = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *cubInterpT   = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  matrixTranspose(cubNq, Nq, cubInterp, Nq, cubInterpT, cubNq);
  matrixTranspose(Nq, cubNq, cubProject, cubNq, cubProjectT, Nq);

  dfloat *cubPDTT     = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  matrixTranspose(Nq, cubNq, cubPDT, cubNq, cubPDTT, Nq);

  o_cubInterp   = device.malloc(Nq*cubNq*sizeof(dfloat), cubInterpT);
  o_cubProject = device.malloc(Nq*cubNq*sizeof(dfloat), cubProjectT);

  o_cubPDT = device.malloc(Nq*cubNq*sizeof(dfloat), cubPDTT);
  o_cubD = device.malloc(cubNq*cubNq*sizeof(dfloat), cubD);

  o_intInterp = o_cubInterp;
  o_intLIFT = o_cubProject;

  free(cubPDTT);
  free(cubProjectT);
  free(cubInterpT);

  cubvgeo = (dfloat*) calloc(Nelements*Nvgeo*cubNp, sizeof(dfloat));
  cubggeo = (dfloat*) calloc(Nelements*Nggeo*cubNp, sizeof(dfloat));

  cubsgeo = (dfloat*) calloc(Nelements*Nsgeo*cubNq*Nfaces, sizeof(dfloat));

  //temp arrays
  dfloat *xe = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *ye = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *ze = (dfloat*) calloc(Np, sizeof(dfloat));
  
  dfloat *xre = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *xse = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *yre = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *yse = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *zre = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *zse = (dfloat*) calloc(Np, sizeof(dfloat));

  dfloat *xre1 = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *xse1 = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *yre1 = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *yse1 = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *zre1 = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *zse1 = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));

  dfloat *xe1 = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *ye1 = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *ze1 = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  
  //geometric data for quadrature
  for(dlong e=0;e<Nelements;++e){ /* for each element */
    for(int j=0;j<Nq;++j){
      for(int i=0;i<Nq;++i){
        int n = i + j*Nq;

	int id = e*Np + n;
	xe[n] = x[id];
	ye[n] = y[id];
	ze[n] = z[id];
	
        //differentiate physical coordinates
        xre[n] = 0.0; xse[n] = 0.0;
        yre[n] = 0.0; yse[n] = 0.0;
	zre[n] = 0.0; zse[n] = 0.0;

        for(int m=0;m<Nq;++m){
          int idr = e*Np + j*Nq + m;
          int ids = e*Np + m*Nq + i;
          xre[n] += D[i*Nq+m]*x[idr];
          xse[n] += D[j*Nq+m]*x[ids];
          yre[n] += D[i*Nq+m]*y[idr];
          yse[n] += D[j*Nq+m]*y[ids];
          zre[n] += D[i*Nq+m]*z[idr];
          zse[n] += D[j*Nq+m]*z[ids];
        }
      }
    }

    //interpolate derivaties to cubature
    for(int j=0;j<Nq;++j){
      for(int i=0;i<cubNq;++i){

        xe1[j*cubNq+i] = 0.0;
        ye1[j*cubNq+i] = 0.0;
	ze1[j*cubNq+i] = 0.0;
	
        xre1[j*cubNq+i] = 0.0; xse1[j*cubNq+i] = 0.0;
        yre1[j*cubNq+i] = 0.0; yse1[j*cubNq+i] = 0.0;
	zre1[j*cubNq+i] = 0.0; zse1[j*cubNq+i] = 0.0;
        for(int n=0;n<Nq;++n){
          xe1[j*cubNq+i] += cubInterp[i*Nq + n]*xe[j*Nq+n];
          ye1[j*cubNq+i] += cubInterp[i*Nq + n]*ye[j*Nq+n];
          ze1[j*cubNq+i] += cubInterp[i*Nq + n]*ze[j*Nq+n];
	  
          xre1[j*cubNq+i] += cubInterp[i*Nq + n]*xre[j*Nq+n];
          xse1[j*cubNq+i] += cubInterp[i*Nq + n]*xse[j*Nq+n];
          yre1[j*cubNq+i] += cubInterp[i*Nq + n]*yre[j*Nq+n];
          yse1[j*cubNq+i] += cubInterp[i*Nq + n]*yse[j*Nq+n];
	  zre1[j*cubNq+i] += cubInterp[i*Nq + n]*zre[j*Nq+n];
          zse1[j*cubNq+i] += cubInterp[i*Nq + n]*zse[j*Nq+n];
        }
      }
    }

    for(int j=0;j<cubNq;++j){
      for(int i=0;i<cubNq;++i){
        dfloat xr = 0.0, xs = 0.0;
        dfloat yr = 0.0, ys = 0.0;
	dfloat zr = 0.0, zs = 0.0;
	dfloat xij = 0.0;
        dfloat yij = 0.0;
	dfloat zij = 0.0;
        for(int n=0;n<Nq;++n){
	  xij += cubInterp[j*Nq + n]*xe1[n*cubNq+i];
          yij += cubInterp[j*Nq + n]*ye1[n*cubNq+i];
	  zij += cubInterp[j*Nq + n]*ze1[n*cubNq+i];

          xr += cubInterp[j*Nq + n]*xre1[n*cubNq+i];
          xs += cubInterp[j*Nq + n]*xse1[n*cubNq+i];
          yr += cubInterp[j*Nq + n]*yre1[n*cubNq+i];
          ys += cubInterp[j*Nq + n]*yse1[n*cubNq+i];
	  zr += cubInterp[j*Nq + n]*zre1[n*cubNq+i];
          zs += cubInterp[j*Nq + n]*zse1[n*cubNq+i];
        }

	dfloat rx = ys*zij - zs*yij; // dXds x X
	dfloat ry = zs*xij - xs*zij;
	dfloat rz = xs*yij - ys*xij;
	
	dfloat sx = zr*yij - yr*zij; // -dXdr x X
	dfloat sy = xr*zij - zr*xij;
	dfloat sz = yr*xij - xr*yij;
	
	dfloat tx = yr*zs - zr*ys; // dXdr x dXds ~ X*|dXdr x dXds|/|X|
	dfloat ty = zr*xs - xr*zs;
	dfloat tz = xr*ys - yr*xs;
	
	dfloat Gx = tx, Gy = ty, Gz = tz;
	
	dfloat J = xij*tx + yij*ty + zij*tz;
	
        /* compute geometric factors for affine coordinate transform*/
        if(J<1e-8) {
          stringstream ss;
          ss << "Negative J found at element " << e << "\n";
          LIBP_ABORT(ss.str())
        }

	rx /= J;      sx /= J;      tx /= J;
	ry /= J;      sy /= J;      ty /= J;
	rz /= J;      sz /= J;      tz /= J;

        dfloat JW = J*cubw[i]*cubw[j];

        /* store geometric factors */
        dlong base = Nvgeo*cubNp*e + i + j*cubNq;
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
	int gbase = Nggeo*Np*e + j*Nq + i;
	ggeo[gbase + cubNp*G00ID] = JW*(rx*rx + ry*ry + rz*rz);
	ggeo[gbase + cubNp*G01ID] = JW*(rx*sx + ry*sy + rz*sz);
	ggeo[gbase + cubNp*G02ID] = JW*(rx*tx + ry*ty + rz*tz);
	ggeo[gbase + cubNp*G11ID] = JW*(sx*sx + sy*sy + sz*sz);
	ggeo[gbase + cubNp*G12ID] = JW*(sx*tx + sy*ty + sz*tz);
	
	ggeo[gbase + cubNp*G22ID] = JW*(tx*tx + ty*ty + tz*tz);
	ggeo[gbase + cubNp*GWJID] = JW;
      }
    }

    for(int f=0;f<Nfaces;++f){ // for each face
      for(int m=0;m<cubNq;++m){  // for each node on face

        //interpolate derivatives of physical coordinates
        dfloat xm = 0, xrm = 0.0, xsm = 0.0;
        dfloat ym = 0, yrm = 0.0, ysm = 0.0;
	dfloat zm = 0, zrm = 0.0, zsm = 0.0;
        for(int n=0;n<Nfp;++n){  // for each node on face
          /* volume index of face node */
          int idn = faceNodes[f*Nfp+n];
	  xm += cubInterp[m*Nq + n]*xe[idn];
          ym += cubInterp[m*Nq + n]*ye[idn];
          zm += cubInterp[m*Nq + n]*ze[idn];
          xrm += cubInterp[m*Nq + n]*xre[idn];
          xsm += cubInterp[m*Nq + n]*xse[idn];
          yrm += cubInterp[m*Nq + n]*yre[idn];
          ysm += cubInterp[m*Nq + n]*yse[idn];
          zrm += cubInterp[m*Nq + n]*zre[idn];
          zsm += cubInterp[m*Nq + n]*zse[idn];
        }

        dfloat txm = yrm*zsm - zrm*ysm;
        dfloat tym = zrm*xsm - xrm*zsm;
        dfloat tzm = xrm*ysm - yrm*xsm;

        dfloat Gx = txm, Gy = tym, Gz = tzm;
	
        /* compute geometric factors for affine coordinate transform*/
        dfloat J = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);

        dfloat nx=0.0, ny=0.0, nz=0.0;

        if(f==0){
          nx = yrm*zm - zrm*ym;
          ny = zrm*xm - xrm*zm;
          nz = xrm*ym - yrm*xm;
        }

        if(f==1){
          nx = ysm*zm - zsm*ym;
          ny = zsm*xm - xsm*zm;
          nz = xsm*ym - ysm*xm;
        }

        if(f==2){
          nx = -yrm*zm + zrm*ym;
          ny = -zrm*xm + xrm*zm;
          nz = -xrm*ym + yrm*xm;
        }

        if(f==3){
          nx = -ysm*zm + zsm*ym;
          ny = -zsm*xm + xsm*zm;
          nz = -xsm*ym + ysm*xm;
        }

        dfloat R = sqrt(xm*xm+ym*ym+zm*zm);

        nx /= R;
        ny /= R;
        nz /= R;

        dfloat sJ = sqrt(nx*nx+ny*ny+nz*nz);

        nx /= sJ;
        ny /= sJ;
        nz /= sJ;
	
        if(sJ<1e-8) {
	  stringstream ss;
	  ss << "Negative sJ found at element " << e << "\n";
	  LIBP_ABORT(ss.str())
        }

        /* output index */
        dlong base = Nsgeo*(Nfaces*cubNq*e + cubNq*f + m);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        cubsgeo[base+NXID] = nx;
        cubsgeo[base+NYID] = ny;
	cubsgeo[base+NZID] = nz;
        cubsgeo[base+SJID] = sJ;
        cubsgeo[base+IJID] = 1./J;

        cubsgeo[base+WIJID] = 1./(J*cubw[0]);
        cubsgeo[base+WSJID] = sJ*cubw[m];
      }
    }
  }

  o_cubvgeo = device.malloc(Nelements*Nvgeo*cubNp*sizeof(dfloat), cubvgeo);
  o_cubggeo = device.malloc(Nelements*Nggeo*cubNp*sizeof(dfloat), cubvgeo);
  o_cubsgeo = device.malloc(Nelements*Nfaces*cubNq*Nsgeo*sizeof(dfloat), cubsgeo);

  free(xre); free(xse); free(yre); free(yse);
  free(xre1); free(xse1); free(yre1); free(yse1);
}
