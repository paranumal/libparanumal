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

/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void meshQuad3D::SurfaceGeometricFactors(){

  /* unified storage array for geometric factors */
  Nsgeo = 14; // fix later
  sgeo = (dfloat*) calloc((Nelements+totalHaloPairs)*
				Nsgeo*Nfp*Nfaces,
				sizeof(dfloat));

  cubsgeo = (dfloat*) calloc((Nelements+totalHaloPairs)*
                                Nsgeo*cubNq*Nfaces,
                                sizeof(dfloat));

  dfloat *cubx = (dfloat*) calloc((Nelements+totalHaloPairs)*
				  cubNq*Nfaces, sizeof(dfloat));

  dfloat *cuby = (dfloat*) calloc((Nelements+totalHaloPairs)*
				  cubNq*Nfaces, sizeof(dfloat));

  dfloat *cubz = (dfloat*) calloc((Nelements+totalHaloPairs)*
				  cubNq*Nfaces, sizeof(dfloat));



  dfloat *xr = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *yr = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *zr = (dfloat*) calloc(Np, sizeof(dfloat));

  dfloat *xs = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *ys = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *zs = (dfloat*) calloc(Np, sizeof(dfloat));

  dfloat *J  = (dfloat*) calloc(Np, sizeof(dfloat));

  for(int e=0;e<Nelements+totalHaloPairs;++e){ /* for each element */

    for(int j=0;j<Nq;++j){
      for(int i=0;i<Nq;++i){

	dfloat xn = x[i+j*Nq+e*Np];
	dfloat yn = y[i+j*Nq+e*Np];
	dfloat zn = z[i+j*Nq+e*Np];

	dfloat xrij = 0, yrij = 0, zrij = 0;
	dfloat xsij = 0, ysij = 0, zsij = 0;

	for(int n=0;n<Nq;++n){

	  dfloat Din = D[i*Nq+n];
	  dfloat Djn = D[j*Nq+n];

	  xrij += Din*x[n+j*Nq+e*Np];
	  yrij += Din*y[n+j*Nq+e*Np];
	  zrij += Din*z[n+j*Nq+e*Np];

	  xsij += Djn*x[i+n*Nq+e*Np];
	  ysij += Djn*y[i+n*Nq+e*Np];
	  zsij += Djn*z[i+n*Nq+e*Np];

	}

	dfloat txij = yrij*zsij - zrij*ysij;
	dfloat tyij = zrij*xsij - xrij*zsij;
	dfloat tzij = xrij*ysij - yrij*xsij;

	dfloat Gx = txij, Gy = tyij, Gz = tzij;

	dfloat Jij = xn*txij + yn*tyij + zn*tzij;

	xr[i+j*Nq] = xrij;
	yr[i+j*Nq] = yrij;
	zr[i+j*Nq] = zrij;

	xs[i+j*Nq] = xsij;
	ys[i+j*Nq] = ysij;
	zs[i+j*Nq] = zsij;

	J[i+j*Nq] = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);
      }
    }

    for(int f=0;f<Nfaces;++f){
      for(int n=0;n<Nq;++n){
	int id = faceNodes[n+f*Nq];

	dfloat xid = x[id+e*Np];
	dfloat yid = y[id+e*Np];
	dfloat zid = z[id+e*Np];
	dfloat Jid = J[id];

	dfloat nx, ny, nz;

	if(f==0){
	  nx = yr[id]*zid - zr[id]*yid;
	  ny = zr[id]*xid - xr[id]*zid;
	  nz = xr[id]*yid - yr[id]*xid;
	}

	if(f==1){
	  nx = ys[id]*zid - zs[id]*yid;
	  ny = zs[id]*xid - xs[id]*zid;
	  nz = xs[id]*yid - ys[id]*xid;
	}

	if(f==2){
	  nx = -yr[id]*zid + zr[id]*yid;
	  ny = -zr[id]*xid + xr[id]*zid;
	  nz = -xr[id]*yid + yr[id]*xid;
	}

	if(f==3){
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

	int base = Nsgeo*(e*Nq*Nfaces + n + f*Nq);

	sgeo[base+NXID] = nx;
	sgeo[base+NYID] = ny;
	sgeo[base+NZID] = nz;
	sgeo[base+SJID] = sJ;

	sgeo[base+IJID] = 1./Jid;

	sgeo[base+WIJID] = 1./(Jid*gllw[0]);
	sgeo[base+WSJID] = sJ*gllw[n];
      }
    }

    // interpolate geofacs to surface quadrature
    for(int f=0;f<Nfaces;++f){

      for(int n=0;n<cubNq;++n){
	dfloat cxr = 0, cxs = 0, cx = 0;
	dfloat cyr = 0, cys = 0, cy = 0;
	dfloat czr = 0, czs = 0, cz = 0;

	for(int i=0;i<Nq;++i){
	  int id = faceNodes[i+f*Nq];
	  dfloat cIni = cubInterp[n*Nq+i];
	  cxr += cIni*xr[id];
	  cxs += cIni*xs[id];
	  cyr += cIni*yr[id];
	  cys += cIni*ys[id];
	  czr += cIni*zr[id];
	  czs += cIni*zs[id];
	  cx  += cIni*x[id+e*Np];
	  cy  += cIni*y[id+e*Np];
	  cz  += cIni*z[id+e*Np];
	}

	cubx[e*cubNq*Nfaces+f*cubNq + n] = cx;
	cuby[e*cubNq*Nfaces+f*cubNq + n] = cy;
	cubz[e*cubNq*Nfaces+f*cubNq + n] = cz;

	dfloat Gx = cyr*czs - czr*cys;
	dfloat Gy = czr*cxs - cxr*czs;
	dfloat Gz = cxr*cys - cyr*cxs;
	dfloat cJ = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);
	dfloat volJ = cx*Gx + cy*Gy + cz*Gz; // xij*tx + yij*ty + zij*tz;
	dfloat nx, ny, nz;

	if(f==0){
	  nx = cyr*cz - czr*cy;
	  ny = czr*cx - cxr*cz;
	  nz = cxr*cy - cyr*cx;
	}

	if(f==1){
	  nx = cys*cz - czs*cy;
	  ny = czs*cx - cxs*cz;
	  nz = cxs*cy - cys*cx;
	}

	if(f==2){
	  nx = -cyr*cz + czr*cy;
	  ny = -czr*cx + cxr*cz;
	  nz = -cxr*cy + cyr*cx;
	}

	if(f==3){
	  nx = -cys*cz + czs*cy;
	  ny = -czs*cx + cxs*cz;
	  nz = -cxs*cy + cys*cx;
	}

	dfloat R = sqrt(cx*cx+cy*cy+cz*cz);

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

	int base = Nsgeo*(e*cubNq*Nfaces + n + f*cubNq);

	cubsgeo[base+NXID] = nx;
	cubsgeo[base+NYID] = ny;
	cubsgeo[base+NZID] = nz;
	cubsgeo[base+SJID] = sJ;
	cubsgeo[base+IHID] = sJ/volJ;
	//	cubsgeo[base+WSJID] = sJ*cubw[n];
      }
    }
  }


#if 0
  for(int e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      for(int n=0;n<Nq;++n){
	int idM = n+f*Nq+e*Nfaces*Nq;
	int idP = mapP[idM];
	int eP = idP/(Nq*Nfaces);
	int fP = (idP%(Nq*Nfaces))/Nq;
	int nP = (idP%Nq);
	int baseM = e*Nq*Nfaces*Nsgeo + f*Nq*Nsgeo + n;
	int baseP = eP*Nq*Nfaces*Nsgeo + fP*Nq*Nsgeo + nP;
	printf("e,f,n=(%d,%d,%d)-(%d,%d,%d): xP-xM=(%g,%g,%g) : norP+norM=%g,%g,%g\n",
	       e,f,n,eP,fP,nP,
	       x[vmapP[idM]]-x[vmapM[idM]],
	       y[vmapP[idM]]-y[vmapM[idM]],
	       z[vmapP[idM]]-z[vmapM[idM]],
	       sgeo[baseM+NXID*Nq]+sgeo[baseP+NXID*Nq],
	       sgeo[baseM+NYID*Nq]+sgeo[baseP+NYID*Nq],
	       sgeo[baseM+NZID*Nq]+sgeo[baseP+NZID*Nq]);

      }
    }
  }
#endif
  // TW: omit 1/min(h) calculation

  for(dlong e=0;e<Nelements;++e){ /* for each non-halo element */
    for(int n=0;n<Nfp*Nfaces;++n){
      dlong baseM = e*Nfp*Nfaces + n;
      dlong baseP = mapP[baseM];
      if(baseP<0) baseP = baseM;

      // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)
      dfloat hinvM = sgeo[baseM*Nsgeo + SJID]*sgeo[baseM*Nsgeo + IJID];
      dfloat hinvP = sgeo[baseP*Nsgeo + SJID]*sgeo[baseP*Nsgeo + IJID];

      //      printf("hinvM/P = %g,%g\n", hinvM, hinvP);

      sgeo[baseM*Nsgeo+IHID] = mymax(hinvM,hinvP);
      sgeo[baseP*Nsgeo+IHID] = mymax(hinvM,hinvP);
    }
  }

  for(dlong e=0;e<Nelements;++e){ /* for each non-halo element */
    for(int f=0;f<Nfaces;++f){
      dlong eP = EToE[e*Nfaces+f];
      dlong fP = EToF[e*Nfaces+f];

      dfloat maxhinv  = 0;
      for(int n=0;n<cubNq;++n){
	dlong idM = e*cubNq*Nfaces+f*cubNq+n;
	dfloat cxM = cubx[idM];
	dfloat cyM = cuby[idM];
	dfloat czM = cubz[idM];

	dfloat mindist2;
	int minidP = 0;
	// jump through hoops to find neighbor cubature node
	// [ not needed elsewhere since we interpolate consistently ]
	dlong idP;
	for(int m=0;m<cubNq;++m){
	  idP = eP*cubNq*Nfaces+fP*cubNq+m;

	  dfloat cxP = cubx[idP];
	  dfloat cyP = cuby[idP];
	  dfloat czP = cubz[idP];

	  dfloat dist2 = pow(cxP-cxM,2)+pow(cyP-cyM,2)+pow(czP-czM,2);

	  if(m==0 || dist2<mindist2){
	    mindist2 = dist2;
	    minidP = m;
	  }
	}

	if(mindist2>1e-12)
	printf("mindist2 = %g\n", mindist2);

	idM = Nsgeo*( e*cubNq*Nfaces+ f*cubNq+n)+IHID;
	idP = Nsgeo*(eP*cubNq*Nfaces+fP*cubNq+minidP)+IHID;

	dfloat hinv = mymax(cubsgeo[idM],cubsgeo[idP]);
	cubsgeo[idM] = hinv;
	cubsgeo[idP] = hinv;

      }
    }
  }


}
