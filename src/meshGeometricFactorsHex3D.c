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


#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

void meshGeometricFactorsHex3D(mesh3D *mesh){

  /* unified storage array for geometric factors */
  mesh->Nvgeo = 12;
  
  /* note that we have volume geometric factors for each node */
  mesh->vgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo*mesh->Np, sizeof(dfloat));

  mesh->cubvgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo*mesh->cubNp, sizeof(dfloat));

  /* number of second order geometric factors */
  mesh->Nggeo = 7;
  mesh->ggeo = (dfloat*) calloc(mesh->Nelements*mesh->Nggeo*mesh->Np, sizeof(dfloat));

  dfloat minJ = 1e9, maxJ = -1e9, maxSkew = 0;
  
  for(dlong e=0;e<mesh->Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*mesh->Nverts;
    
    dfloat *xe = mesh->EX + id;
    dfloat *ye = mesh->EY + id;
    dfloat *ze = mesh->EZ + id;

    for(int k=0;k<mesh->Nq;++k){
      for(int j=0;j<mesh->Nq;++j){
        for(int i=0;i<mesh->Nq;++i){
          
          int n = i + j*mesh->Nq + k*mesh->Nq*mesh->Nq;

          /* local node coordinates */
          dfloat rn = mesh->r[n]; 
          dfloat sn = mesh->s[n];
          dfloat tn = mesh->t[n];
          
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
          
          if(J<1e-12) printf("J = %g !!!!!!!!!!!!!\n", J);
          
          dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
          dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
          dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;
          
          dfloat JW = J*mesh->gllw[i]*mesh->gllw[j]*mesh->gllw[k];
          
          /* store geometric factors */
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*RXID] = rx;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*RYID] = ry;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*RZID] = rz;
          
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*SXID] = sx;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*SYID] = sy;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*SZID] = sz;
          
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*TXID] = tx;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*TYID] = ty;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*TZID] = tz;
          
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*JID]  = J;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*JWID] = JW;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*IJWID] = 1./JW;

          /* store second order geometric factors */
          mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G00ID] = JW*(rx*rx + ry*ry + rz*rz);
          mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G01ID] = JW*(rx*sx + ry*sy + rz*sz);
          mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G02ID] = JW*(rx*tx + ry*ty + rz*tz);
          mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G11ID] = JW*(sx*sx + sy*sy + sz*sz);
          mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G12ID] = JW*(sx*tx + sy*ty + sz*tz);
          mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G22ID] = JW*(tx*tx + ty*ty + tz*tz);
          mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*GWJID] = JW;
        }
      }
    }

    //geometric data for quadrature
    for(int k=0;k<mesh->cubNq;++k){
      for(int j=0;j<mesh->cubNq;++j){
        for(int i=0;i<mesh->cubNq;++i){
          
          dfloat rn = mesh->cubr[i];
          dfloat sn = mesh->cubr[j];
          dfloat tn = mesh->cubr[k];

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
          
          /* compute geometric factors for affine coordinate transform*/
          dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
          
          dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
          dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
          dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;
          
          dfloat JW = J*mesh->cubw[i]*mesh->cubw[j]*mesh->cubw[k];
          
          /* store geometric factors */
          dlong base = mesh->Nvgeo*mesh->cubNp*e + i + j*mesh->cubNq + k*mesh->cubNq*mesh->cubNq;
          mesh->cubvgeo[base + mesh->cubNp*RXID] = rx;
          mesh->cubvgeo[base + mesh->cubNp*RYID] = ry;
          mesh->cubvgeo[base + mesh->cubNp*RZID] = rz;
          
          mesh->cubvgeo[base + mesh->cubNp*SXID] = sx;
          mesh->cubvgeo[base + mesh->cubNp*SYID] = sy;
          mesh->cubvgeo[base + mesh->cubNp*SZID] = sz;
          
          mesh->cubvgeo[base + mesh->cubNp*TXID] = tx;
          mesh->cubvgeo[base + mesh->cubNp*TYID] = ty;
          mesh->cubvgeo[base + mesh->cubNp*TZID] = tz;
          
          mesh->cubvgeo[base + mesh->cubNp*JID]  = J;
          mesh->cubvgeo[base + mesh->cubNp*JWID] = JW;
          mesh->cubvgeo[base + mesh->cubNp*IJWID] = 1./JW;
        }
      }
    }
  }

  printf("J in range [%g,%g] and max Skew = %g\n", minJ, maxJ, maxSkew);
}
