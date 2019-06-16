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

static void interpolateHex3D(dfloat *I, dfloat *x, int N, dfloat *Ix, int M){

  dfloat *Ix1 = (dfloat*) calloc(N*N*M, sizeof(dfloat));
  dfloat *Ix2 = (dfloat*) calloc(N*M*M, sizeof(dfloat));

  for(int k=0;k<N;++k){
    for(int j=0;j<N;++j){
      for(int i=0;i<M;++i){
        dfloat tmp = 0;
        for(int n=0;n<N;++n){
          tmp += I[i*N + n]*x[k*N*N+j*N+n];
        }
        Ix1[k*N*M+j*M+i] = tmp;
      }
    }
  }

  for(int k=0;k<N;++k){
    for(int j=0;j<M;++j){
      for(int i=0;i<M;++i){
        dfloat tmp = 0;
        for(int n=0;n<N;++n){
          tmp += I[j*N + n]*Ix1[k*N*M+n*M+i];
        }
        Ix2[k*M*M+j*M+i] = tmp;
      }
    }
  }

  for(int k=0;k<M;++k){
    for(int j=0;j<M;++j){
      for(int i=0;i<M;++i){
        dfloat tmp = 0;
        for(int n=0;n<N;++n){
          tmp += I[k*N + n]*Ix2[n*M*M+j*M+i];
        }
        Ix[k*M*M+j*M+i] = tmp;
      }
    }
  }

  free(Ix1);
  free(Ix2);

}

void meshHex3D::GeometricFactors(){

  /* unified storage array for geometric factors */
  Nvgeo = 12;

  /* note that we have volume geometric factors for each node */
  vgeo = (dfloat*) calloc((Nelements+totalHaloPairs)*Nvgeo*Np, sizeof(dfloat));

  cubvgeo = (dfloat*) calloc(Nelements*Nvgeo*cubNp, sizeof(dfloat));

  /* number of second order geometric factors */
  Nggeo = 7;

  ggeo    = (dfloat*) calloc(Nelements*Nggeo*Np,    sizeof(dfloat));
  cubggeo = (dfloat*) calloc(Nelements*Nggeo*cubNp, sizeof(dfloat));

  dfloat minJ = 1e9, maxJ = -1e9, maxSkew = 0;

  dfloat *xre = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *xse = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *xte = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *yre = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *yse = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *yte = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *zre = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *zse = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *zte = (dfloat*) calloc(Np, sizeof(dfloat));

  dfloat *cubxre = (dfloat*) calloc(cubNp, sizeof(dfloat));
  dfloat *cubxse = (dfloat*) calloc(cubNp, sizeof(dfloat));
  dfloat *cubxte = (dfloat*) calloc(cubNp, sizeof(dfloat));
  dfloat *cubyre = (dfloat*) calloc(cubNp, sizeof(dfloat));
  dfloat *cubyse = (dfloat*) calloc(cubNp, sizeof(dfloat));
  dfloat *cubyte = (dfloat*) calloc(cubNp, sizeof(dfloat));
  dfloat *cubzre = (dfloat*) calloc(cubNp, sizeof(dfloat));
  dfloat *cubzse = (dfloat*) calloc(cubNp, sizeof(dfloat));
  dfloat *cubzte = (dfloat*) calloc(cubNp, sizeof(dfloat));

  for(dlong e=0;e<Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    // dlong id = e*Nverts;

    // dfloat *xe = EX + id;
    // dfloat *ye = EY + id;
    // dfloat *ze = EZ + id;

    for(int n=0;n<Np;++n){
      xre[n] = 0; xse[n] = 0; xte[n] = 0;
      yre[n] = 0; yse[n] = 0; yte[n] = 0;
      zre[n] = 0; zse[n] = 0; zte[n] = 0;
    }

    for(int k=0;k<Nq;++k){
      for(int j=0;j<Nq;++j){
        for(int i=0;i<Nq;++i){

          int n = i + j*Nq + k*Nq*Nq;

          /* local node coordinates */
          // dfloat rn = r[n];
          // dfloat sn = s[n];
          // dfloat tn = t[n];

#if 0
          /* Jacobian matrix */
          dfloat xr = 0.125*( (1-tn)*(1-sn)*(xe[1]-xe[0]) + (1-tn)*(1+sn)*(xe[2]-xe[3]) + (1+tn)*(1-sn)*(xe[5]-xe[4]) + (1+tn)*(1+sn)*(xe[6]-xe[7]) );
          dfloat xs = 0.125*( (1-tn)*(1-rn)*(xe[3]-xe[0]) + (1-tn)*(1+rn)*(xe[2]-xe[1]) + (1+tn)*(1-rn)*(xe[7]-xe[4]) + (1+tn)*(1+rn)*(xe[6]-xe[5]) );
          dfloat xt = 0.125*( (1-rn)*(1-sn)*(xe[4]-xe[0]) + (1+rn)*(1-sn)*(xe[5]-xe[1]) + (1+rn)*(1+sn)*(xe[6]-xe[2]) + (1-rn)*(1+sn)*(xe[7]-xe[3]) );

          dfloat yr = 0.125*( (1-tn)*(1-sn)*(ye[1]-ye[0]) + (1-tn)*(1+sn)*(ye[2]-ye[3]) + (1+tn)*(1-sn)*(ye[5]-ye[4]) + (1+tn)*(1+sn)*(ye[6]-ye[7]) );
          dfloat ys = 0.125*( (1-tn)*(1-rn)*(ye[3]-ye[0]) + (1-tn)*(1+rn)*(ye[2]-ye[1]) + (1+tn)*(1-rn)*(ye[7]-ye[4]) + (1+tn)*(1+rn)*(ye[6]-ye[5]) );
          dfloat yt = 0.125*( (1-rn)*(1-sn)*(ye[4]-ye[0]) + (1+rn)*(1-sn)*(ye[5]-ye[1]) + (1+rn)*(1+sn)*(ye[6]-ye[2]) + (1-rn)*(1+sn)*(ye[7]-ye[3]) );

          dfloat zr = 0.125*( (1-tn)*(1-sn)*(ze[1]-ze[0]) + (1-tn)*(1+sn)*(ze[2]-ze[3]) + (1+tn)*(1-sn)*(ze[5]-ze[4]) + (1+tn)*(1+sn)*(ze[6]-ze[7]) );
          dfloat zs = 0.125*( (1-tn)*(1-rn)*(ze[3]-ze[0]) + (1-tn)*(1+rn)*(ze[2]-ze[1]) + (1+tn)*(1-rn)*(ze[7]-ze[4]) + (1+tn)*(1+rn)*(ze[6]-ze[5]) );
          dfloat zt = 0.125*( (1-rn)*(1-sn)*(ze[4]-ze[0]) + (1+rn)*(1-sn)*(ze[5]-ze[1]) + (1+rn)*(1+sn)*(ze[6]-ze[2]) + (1-rn)*(1+sn)*(ze[7]-ze[3]) );
#else
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

          dfloat xr = xre[n], xs = xse[n], xt = xte[n];
          dfloat yr = yre[n], ys = yse[n], yt = yte[n];
          dfloat zr = zre[n], zs = zse[n], zt = zte[n];
#endif

          /* compute geometric factors for affine coordinate transform*/
          dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);

          dfloat hr = sqrt(xr*xr+yr*yr+zr*zr);
          dfloat hs = sqrt(xs*xs+ys*ys+zs*zs);
          dfloat ht = sqrt(xt*xt+yt*yt+zt*zt);
          minJ = mymin(J, minJ);
          maxJ = mymax(J, maxJ);
          maxSkew = mymax(maxSkew, hr/hs);
          maxSkew = mymax(maxSkew, hr/ht);
          maxSkew = mymax(maxSkew, hs/hr);
          maxSkew = mymax(maxSkew, hs/ht);
          maxSkew = mymax(maxSkew, ht/hr);
          maxSkew = mymax(maxSkew, ht/hs);

          if(J<1e-12) {
            stringstream ss;
            ss << "Negative J found at element " << e << "\n";
            LIBP_ABORT(ss.str())
          }

          dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
          dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
          dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;

          dfloat JW = J*gllw[i]*gllw[j]*gllw[k];

          /* store geometric factors */
          vgeo[Nvgeo*Np*e + n + Np*RXID] = rx;
          vgeo[Nvgeo*Np*e + n + Np*RYID] = ry;
          vgeo[Nvgeo*Np*e + n + Np*RZID] = rz;

          vgeo[Nvgeo*Np*e + n + Np*SXID] = sx;
          vgeo[Nvgeo*Np*e + n + Np*SYID] = sy;
          vgeo[Nvgeo*Np*e + n + Np*SZID] = sz;

          vgeo[Nvgeo*Np*e + n + Np*TXID] = tx;
          vgeo[Nvgeo*Np*e + n + Np*TYID] = ty;
          vgeo[Nvgeo*Np*e + n + Np*TZID] = tz;

          vgeo[Nvgeo*Np*e + n + Np*JID]  = J;
          vgeo[Nvgeo*Np*e + n + Np*JWID] = JW;
          vgeo[Nvgeo*Np*e + n + Np*IJWID] = 1./JW;

          /* store second order geometric factors */
          ggeo[Nggeo*Np*e + n + Np*G00ID] = JW*(rx*rx + ry*ry + rz*rz);
          ggeo[Nggeo*Np*e + n + Np*G01ID] = JW*(rx*sx + ry*sy + rz*sz);
          ggeo[Nggeo*Np*e + n + Np*G02ID] = JW*(rx*tx + ry*ty + rz*tz);
          ggeo[Nggeo*Np*e + n + Np*G11ID] = JW*(sx*sx + sy*sy + sz*sz);
          ggeo[Nggeo*Np*e + n + Np*G12ID] = JW*(sx*tx + sy*ty + sz*tz);
          ggeo[Nggeo*Np*e + n + Np*G22ID] = JW*(tx*tx + ty*ty + tz*tz);
          ggeo[Nggeo*Np*e + n + Np*GWJID] = JW;
        }
      }
    }

    interpolateHex3D(cubInterp, xre, Nq, cubxre, cubNq);
    interpolateHex3D(cubInterp, xse, Nq, cubxse, cubNq);
    interpolateHex3D(cubInterp, xte, Nq, cubxte, cubNq);

    interpolateHex3D(cubInterp, yre, Nq, cubyre, cubNq);
    interpolateHex3D(cubInterp, yse, Nq, cubyse, cubNq);
    interpolateHex3D(cubInterp, yte, Nq, cubyte, cubNq);

    interpolateHex3D(cubInterp, zre, Nq, cubzre, cubNq);
    interpolateHex3D(cubInterp, zse, Nq, cubzse, cubNq);
    interpolateHex3D(cubInterp, zte, Nq, cubzte, cubNq);

    //geometric data for quadrature
    for(int k=0;k<cubNq;++k){
      for(int j=0;j<cubNq;++j){
        for(int i=0;i<cubNq;++i){

          int n = k*cubNq*cubNq + j*cubNq + i;

          // dfloat rn = cubr[i];
          // dfloat sn = cubr[j];
          // dfloat tn = cubr[k];

          /* Jacobian matrix */
          dfloat xr = cubxre[n], xs = cubxse[n], xt = cubxte[n];
          dfloat yr = cubyre[n], ys = cubyse[n], yt = cubyte[n];
          dfloat zr = cubzre[n], zs = cubzse[n], zt = cubzte[n];

          /* compute geometric factors for affine coordinate transform*/
          dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);

          if(J<1e-12) {
            stringstream ss;
            ss << "Negative cubature J found at element " << e << "\n";
            LIBP_ABORT(ss.str())
          }

          dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
          dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
          dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;

          dfloat JW = J*cubw[i]*cubw[j]*cubw[k];

          /* store geometric factors */
          dlong base = Nvgeo*cubNp*e + n;
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
          base = Nggeo*cubNp*e + n;
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
  }

  #if 0
    dfloat globalMinJ, globalMaxJ, globalMaxSkew;

    MPI_Reduce(&minJ, &globalMinJ, 1, MPI_DFLOAT, MPI_MIN, 0, comm);
    MPI_Reduce(&maxJ, &globalMaxJ, 1, MPI_DFLOAT, MPI_MAX, 0, comm);
    MPI_Reduce(&maxSkew, &globalMaxSkew, 1, MPI_DFLOAT, MPI_MAX, 0, comm);

    if(rank==0)
      printf("J in range [%g,%g] and max Skew = %g\n", globalMinJ, globalMaxJ, globalMaxSkew);
  #endif

  HaloExchange(vgeo, Nvgeo*Np, ogsDfloat);

  free(xre); free(xse); free(xte);
  free(yre); free(yse); free(yte);
  free(zre); free(zse); free(zte);

  free(cubxre); free(cubxse); free(cubxte);
  free(cubyre); free(cubyse); free(cubyte);
  free(cubzre); free(cubzse); free(cubzte);

}
