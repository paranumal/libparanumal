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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

void computeFrame(dfloat nx, dfloat ny, dfloat nz,
		  dfloat &tanx, dfloat &tany, dfloat &tanz,
		  dfloat &binx, dfloat &biny, dfloat &binz){

  dfloat rdotn, ranx, rany, ranz;
  do{
    ranx = drand48();
    rany = drand48();
    ranz = drand48();
    
    dfloat magran = sqrt(ranx*ranx+rany*rany+ranz*ranz);
    
    ranx /= magran;
    rany /= magran;
    ranz /= magran;
    
    rdotn = nx*ranx+ny*rany+nz*ranz;
  }while(fabs(rdotn)<1e-4);
  
  tanx = ny*ranz - nz*rany;
  tany = nz*ranx - nx*ranz;
  tanz = nx*rany - ny*ranx;

  dfloat magtan = sqrt(tanx*tanx+tany*tany+tanz*tanz);

  tanx /= magtan;
  tany /= magtan;
  tanz /= magtan;

  binx = ny*tanz - nz*tany;
  biny = nz*tanx - nx*tanz;
  binz = nx*tany - ny*tanx;

  dfloat magbin = sqrt(binx*binx+biny*biny+binz*binz);

  binx /= magbin;
  biny /= magbin;
  binz /= magbin;

  //  printf("nor = %g,%g,%g; tan = %g,%g,%g; bin = %g,%g,%g\n", nx, ny, nz, tanx, tany, tanz, binx, biny, binz);
}


/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void meshSurfaceGeometricFactorsHex3D(mesh3D *mesh){

  /* unified storage array for geometric factors */
  mesh->Nsgeo = 14;
  mesh->sgeo = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
                                mesh->Nsgeo*mesh->Nfp*mesh->Nfaces, 
                                sizeof(dfloat));

  mesh->cubsgeo = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
                                mesh->Nsgeo*mesh->cubNfp*mesh->Nfaces, 
                                sizeof(dfloat));

  for(dlong e=0;e<mesh->Nelements+mesh->totalHaloPairs;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*mesh->Nverts;

    dfloat *xe = mesh->EX + id;
    dfloat *ye = mesh->EY + id;
    dfloat *ze = mesh->EZ + id;
    
    for(int f=0;f<mesh->Nfaces;++f){ // for each face
      
      for(int i=0;i<mesh->Nfp;++i){  // for each node on face

        /* volume index of face node */
        int n = mesh->faceNodes[f*mesh->Nfp+i];

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

        /* determinant of Jacobian matrix */
        dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
        
        dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
        dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
        dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;
        
        /* face f normal and length */
        dfloat nx, ny, nz, d;
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
        dlong base = mesh->Nsgeo*(mesh->Nfaces*mesh->Nfp*e + mesh->Nfp*f + i);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        mesh->sgeo[base+NXID] = nx;
        mesh->sgeo[base+NYID] = ny;
        mesh->sgeo[base+NZID] = nz;
        mesh->sgeo[base+SJID] = sJ;
        mesh->sgeo[base+IJID] = 1./J;

        mesh->sgeo[base+WIJID] = 1./(J*mesh->gllw[0]);
        mesh->sgeo[base+WSJID] = sJ*mesh->gllw[i%mesh->Nq]*mesh->gllw[i/mesh->Nq];

	computeFrame(nx, ny, nz,
		     mesh->sgeo[base+STXID], mesh->sgeo[base+STYID], mesh->sgeo[base+STZID],
		     mesh->sgeo[base+SBXID], mesh->sgeo[base+SBYID], mesh->sgeo[base+SBZID]);
      }

      //geometric data for quadrature
      for(int i=0;i<mesh->cubNfp;++i){  // for each quadrature node on face

        dfloat rn, sn, tn;
        switch(f){
        case 0: rn = mesh->cubr[i%mesh->cubNq]; sn = mesh->cubr[i/mesh->cubNq]; tn = -1.0;                      break;
        case 1: rn = mesh->cubr[i%mesh->cubNq]; sn = -1.0;                      tn = mesh->cubr[i/mesh->cubNq]; break;
        case 2: rn = 1.0;                       sn = mesh->cubr[i%mesh->cubNq]; tn = mesh->cubr[i/mesh->cubNq]; break;
        case 3: rn = mesh->cubr[i%mesh->cubNq]; sn = 1.0;                       tn = mesh->cubr[i/mesh->cubNq]; break;
        case 4: rn = -1.0;                      sn = mesh->cubr[i%mesh->cubNq]; tn = mesh->cubr[i/mesh->cubNq]; break;
        case 5: rn = mesh->cubr[i%mesh->cubNq]; sn = mesh->cubr[i/mesh->cubNq]; tn = 1.0;                       break;
        }

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

        /* determinant of Jacobian matrix */
        dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
        
        dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
        dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
        dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;
        
        /* face f normal and length */
        dfloat nx, ny, nz, d;
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
        dlong base = mesh->Nsgeo*(mesh->Nfaces*mesh->cubNfp*e + mesh->cubNfp*f + i);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        mesh->cubsgeo[base+NXID] = nx;
        mesh->cubsgeo[base+NYID] = ny;
        mesh->cubsgeo[base+NZID] = nz;
        mesh->cubsgeo[base+SJID] = sJ;
        mesh->cubsgeo[base+IJID] = 1./J;

        mesh->cubsgeo[base+WIJID] = 1./(J*mesh->cubw[0]);
        mesh->cubsgeo[base+WSJID] = sJ*mesh->cubw[i%mesh->cubNq]*mesh->cubw[i/mesh->cubNq];

	computeFrame(nx, ny, nz,
		     mesh->cubsgeo[base+STXID], mesh->cubsgeo[base+STYID], mesh->cubsgeo[base+STZID],
		     mesh->cubsgeo[base+SBXID], mesh->cubsgeo[base+SBYID], mesh->cubsgeo[base+SBZID]);
	
	
      }
    }
  }

  for(dlong e=0;e<mesh->Nelements;++e){ /* for each non-halo element */
    for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      dlong baseM = e*mesh->Nfp*mesh->Nfaces + n;
      dlong baseP = mesh->mapP[baseM];
      // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)
      dfloat hinvM = mesh->sgeo[baseM*mesh->Nsgeo + SJID]*mesh->sgeo[baseM*mesh->Nsgeo + IJID];
      dfloat hinvP = mesh->sgeo[baseP*mesh->Nsgeo + SJID]*mesh->sgeo[baseP*mesh->Nsgeo + IJID];
      mesh->sgeo[baseM*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
      mesh->sgeo[baseP*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
    }
  }
}
