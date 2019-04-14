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

/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void meshTri3D::SurfaceGeometricFactors(){

  /* unified storage array for geometric factors */
  Nsgeo = 14;
  sgeo = (dfloat*) calloc((Nelements+totalHaloPairs)*
				Nsgeo*Nfp*Nfaces,
				sizeof(dfloat));

  dfloat *xr = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *yr = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *zr = (dfloat*) calloc(Np, sizeof(dfloat));

  dfloat *xs = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *ys = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *zs = (dfloat*) calloc(Np, sizeof(dfloat));

  dfloat *J  = (dfloat*) calloc(Np, sizeof(dfloat));

  for(int e=0;e<Nelements+totalHaloPairs;++e){ /* for each element */

    for(int n=0;n<Np;++n){

      // dfloat x1 = x[n+e*Np];
      // dfloat y1 = y[n+e*Np];
      // dfloat z1 = z[n+e*Np];

      dfloat xrn = 0, yrn = 0, zrn = 0;
      dfloat xsn = 0, ysn = 0, zsn = 0;

      for(int m=0;m<Np;++m){

	dfloat Drn = Dr[n*Np+m];
	dfloat Dsn = Ds[n*Np+m];

	xrn += Drn*x[m+e*Np];
	yrn += Drn*y[m+e*Np];
	zrn += Drn*z[m+e*Np];

	xsn += Dsn*x[m+e*Np];
	ysn += Dsn*y[m+e*Np];
	zsn += Dsn*z[m+e*Np];

      }

      dfloat txn = yrn*zsn - zrn*ysn;
      dfloat tyn = zrn*xsn - xrn*zsn;
      dfloat tzn = xrn*ysn - yrn*xsn;

      dfloat Gx = txn, Gy = tyn, Gz = tzn;

      // dfloat Jn = x1*txn + y1*tyn + z1*tzn;

      xr[n] = xrn;
      yr[n] = yrn;
      zr[n] = zrn;

      xs[n] = xsn;
      ys[n] = ysn;
      zs[n] = zsn;

      J[n] = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);
    }

    for(int f=0;f<Nfaces;++f){
      for(int n=0;n<Nfp;++n){
	int id = faceNodes[n+f*Nfp];

	dfloat xid = x[id+e*Np];
	dfloat yid = y[id+e*Np];
	dfloat zid = z[id+e*Np];
	dfloat Jid = J[id];

	dfloat nx=0.0, ny=0.0, nz=0.0;

	if(f==0){
	  nx = yr[id]*zid - zr[id]*yid;
	  ny = zr[id]*xid - xr[id]*zid;
	  nz = xr[id]*yid - yr[id]*xid;
	}

	if(f==1){
	  nx = (ys[id]-yr[id])*zid - (zs[id]-zr[id])*yid;
	  ny = (zs[id]-zr[id])*xid - (xs[id]-xr[id])*zid;
	  nz = (xs[id]-xr[id])*yid - (ys[id]-yr[id])*xid;
	}

	if(f==2){
	  nx = -ys[id]*zid + zs[id]*yid;
	  ny = -zs[id]*xid + xs[id]*zid;
	  nz = -xs[id]*yid + ys[id]*xid;
	}

	dfloat R = sqrt(xid*xid+yid*yid+zid*zid);

	nx /= R;
	ny /= R;
	nz /= R;

	dfloat sJ = sqrt(nx*nx+ny*ny+nz*nz);

	nx /= sJ;
	ny /= sJ;
	nz /= sJ;

	if(sJ<1e-8) {
		stringstream ss;
		ss << "Negative J found at element " << e << "\n";
		LIBP_ABORT(ss.str())
	}

	int base = e*Nfp*Nfaces*Nsgeo + n + f*Nfp;

	sgeo[base+Nfp*Nfaces*NXID] = nx;
	sgeo[base+Nfp*Nfaces*NYID] = ny;
	sgeo[base+Nfp*Nfaces*NZID] = nz;
	sgeo[base+Nfp*Nfaces*SJID] = sJ;

	sgeo[base+Nfp*Nfaces*IJID] = 1./Jid;
      }
    }
  }

#if 0
  for(int e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      for(int n=0;n<Nfp;++n){
	int idM = n+f*Nfp+e*Nfaces*Nfp;
	int idP = mapP[idM];
	int eP = idP/(Nfp*Nfaces);
	int fP = (idP%(Nfp*Nfaces))/Nfp;
	int nP = (idP%Nfp);
	int baseM = e*Nfp*Nfaces*Nsgeo + f*Nfp + n;
	int baseP = eP*Nfp*Nfaces*Nsgeo + fP*Nfp + nP;
	printf("e,f,n=(%d,%d,%d)-(%d,%d,%d): xP-xM=(%g,%g,%g) : norP+norM=%g,%g,%g\n",
	       e,f,n,eP,fP,nP,
	       x[vmapP[idM]]-x[vmapM[idM]],
	       y[vmapP[idM]]-y[vmapM[idM]],
	       z[vmapP[idM]]-z[vmapM[idM]],
	       sgeo[baseM+NXID*Nfp*Nfaces]+sgeo[baseP+NXID*Nfp*Nfaces],
	       sgeo[baseM+NYID*Nfp*Nfaces]+sgeo[baseP+NYID*Nfp*Nfaces],
	       sgeo[baseM+NZID*Nfp*Nfaces]+sgeo[baseP+NZID*Nfp*Nfaces]);

      }
    }
  }
#endif
  // TW: omit 1/min(h) calculation
}
