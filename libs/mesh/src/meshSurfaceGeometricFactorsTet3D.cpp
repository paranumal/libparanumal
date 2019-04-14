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

static void computeFrameTet3D(dfloat nx, dfloat ny, dfloat nz,
		  dfloat &tanx, dfloat &tany, dfloat &tanz,
		  dfloat &binx, dfloat &biny, dfloat &binz){

  dfloat ranx = drand48();
  dfloat rany = drand48();
  dfloat ranz = drand48();

  dfloat magran = sqrt(ranx*ranx+rany*rany+ranz*ranz);

  ranx /= magran;
  rany /= magran;
  ranz /= magran;

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

void meshTet3D::SurfaceGeometricFactors(){

  /* unified storage array for geometric factors */
  Nsgeo = 14;
  sgeo = (dfloat*) calloc((Nelements+totalHaloPairs)*
                            Nsgeo*Nfaces, sizeof(dfloat));

  for(dlong e=0;e<Nelements+totalHaloPairs;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*Nverts;
    dfloat xe1 = EX[id+0], ye1 = EY[id+0], ze1 = EZ[id+0];
    dfloat xe2 = EX[id+1], ye2 = EY[id+1], ze2 = EZ[id+1];
    dfloat xe3 = EX[id+2], ye3 = EY[id+2], ze3 = EZ[id+2];
    dfloat xe4 = EX[id+3], ye4 = EY[id+3], ze4 = EZ[id+3];

    /* Jacobian matrix */
    dfloat xr = 0.5*(xe2-xe1), xs = 0.5*(xe3-xe1), xt = 0.5*(xe4-xe1);
    dfloat yr = 0.5*(ye2-ye1), ys = 0.5*(ye3-ye1), yt = 0.5*(ye4-ye1);
    dfloat zr = 0.5*(ze2-ze1), zs = 0.5*(ze3-ze1), zt = 0.5*(ze4-ze1);

    /* compute geometric factors for affine coordinate transform*/
    dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
    dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
    dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
    dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;

    if(J<0) {
      stringstream ss;
      ss << "Negative J found at element " << e << "\n";
      LIBP_ABORT(ss.str())
    }

    /* face 1 */
    dlong base = Nsgeo*Nfaces*e;
    dfloat nx1 = -tx;
    dfloat ny1 = -ty;
    dfloat nz1 = -tz;
    dfloat sJ1 = norm3(nx1,ny1,nz1);

    sgeo[base+NXID] = nx1/sJ1;
    sgeo[base+NYID] = ny1/sJ1;
    sgeo[base+NZID] = nz1/sJ1;
    sgeo[base+SJID] = sJ1*J;
    sgeo[base+IJID] = 1./J;

    // generate local tangent and binormal using random vector
    computeFrameTet3D(nx1/sJ1, ny1/sJ1, nz1/sJ1,
		 sgeo[base+STXID], sgeo[base+STYID], sgeo[base+STZID],
		 sgeo[base+SBXID], sgeo[base+SBYID], sgeo[base+SBZID]);

    /* face 2 */
    base += Nsgeo;
    dfloat nx2 = -sx;
    dfloat ny2 = -sy;
    dfloat nz2 = -sz;
    dfloat sJ2 = norm3(nx2,ny2,nz2);

    sgeo[base+NXID] = nx2/sJ2;
    sgeo[base+NYID] = ny2/sJ2;
    sgeo[base+NZID] = nz2/sJ2;
    sgeo[base+SJID] = sJ2*J;
    sgeo[base+IJID] = 1./J;

    // generate local tangent and binormal using random vector
    computeFrameTet3D(nx2/sJ2, ny2/sJ2, nz2/sJ2,
		 sgeo[base+STXID], sgeo[base+STYID], sgeo[base+STZID],
		 sgeo[base+SBXID], sgeo[base+SBYID], sgeo[base+SBZID]);

    /* face 3 */
    base += Nsgeo;
    dfloat nx3 = rx+sx+tx;
    dfloat ny3 = ry+sy+ty;
    dfloat nz3 = rz+sz+tz;
    dfloat sJ3 = norm3(nx3,ny3,nz3);

    sgeo[base+NXID] = nx3/sJ3;
    sgeo[base+NYID] = ny3/sJ3;
    sgeo[base+NZID] = nz3/sJ3;
    sgeo[base+SJID] = sJ3*J;
    sgeo[base+IJID] = 1./J;

    // generate local tangent and binormal using random vector
    computeFrameTet3D(nx3/sJ3, ny3/sJ3, nz3/sJ3,
		 sgeo[base+STXID], sgeo[base+STYID], sgeo[base+STZID],
		 sgeo[base+SBXID], sgeo[base+SBYID], sgeo[base+SBZID]);

    /* face 4 */
    base += Nsgeo;
    dfloat nx4 = -rx;
    dfloat ny4 = -ry;
    dfloat nz4 = -rz;
    dfloat sJ4 = norm3(nx4,ny4,nz4);

    sgeo[base+NXID] = nx4/sJ4;
    sgeo[base+NYID] = ny4/sJ4;
    sgeo[base+NZID] = nz4/sJ4;
    sgeo[base+SJID] = sJ4*J;
    sgeo[base+IJID] = 1./J;

    // generate local tangent and binormal using random vector
    computeFrameTet3D(nx4/sJ4, ny4/sJ4, nz4/sJ4,
		 sgeo[base+STXID], sgeo[base+STYID], sgeo[base+STZID],
		 sgeo[base+SBXID], sgeo[base+SBYID], sgeo[base+SBZID]);

#if 0
    printf("N1=(%g,%g,%g),sJ1=%g\n", nx1/sJ1,ny1/sJ1,nz1/sJ1,sJ1*J);
    printf("N2=(%g,%g,%g),sJ2=%g\n", nx2/sJ2,ny2/sJ2,nz2/sJ2,sJ2*J);
    printf("N3=(%g,%g,%g),sJ3=%g\n", nx3/sJ3,ny3/sJ3,nz3/sJ3,sJ3*J);
    printf("N4=(%g,%g,%g),sJ4=%g\n", nx4/sJ4,ny4/sJ4,nz4/sJ4,sJ4*J);
#endif
  }

  for(dlong e=0;e<Nelements;++e){ /* for each non-halo element */
    for(int f=0;f<Nfaces;++f){
      dlong baseM = e*Nfaces + f;

      // awkward: (need to find eP,fP relative to bulk+halo)
      dlong idP = vmapP[e*Nfp*Nfaces+f*Nfp+0];
      dlong eP = (idP>=0) ? (idP/Np):e;

      int fP = EToF[baseM];
      fP = (fP==-1) ? f:fP;

      dlong baseP = eP*Nfaces + fP;

      // rescaling,  V = A*h/3 => (J*4/3) = (sJ*2)*h/3 => h  = 0.5*J/sJ
      dfloat hinvM = 0.5*sgeo[baseM*Nsgeo + SJID]*sgeo[baseM*Nsgeo + IJID];
      dfloat hinvP = 0.5*sgeo[baseP*Nsgeo + SJID]*sgeo[baseP*Nsgeo + IJID];

      sgeo[baseM*Nsgeo+IHID] = mymax(hinvM,hinvP);
      sgeo[baseP*Nsgeo+IHID] = mymax(hinvM,hinvP);
    }
  }
}
