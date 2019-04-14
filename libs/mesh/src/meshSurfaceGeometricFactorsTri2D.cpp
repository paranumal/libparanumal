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

void meshTri2D::SurfaceGeometricFactors(){

  /* unified storage array for geometric factors */
  Nsgeo = 6;
  sgeo = (dfloat*) calloc((Nelements+totalHaloPairs)*
				Nsgeo*Nfaces,
				sizeof(dfloat));

  for(dlong e=0;e<Nelements+totalHaloPairs;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*Nverts;

    dfloat xe1 = EX[id+0];
    dfloat xe2 = EX[id+1];
    dfloat xe3 = EX[id+2];

    dfloat ye1 = EY[id+0];
    dfloat ye2 = EY[id+1];
    dfloat ye3 = EY[id+2];

    /* compute geometric factors for affine coordinate transform*/
    dfloat J = 0.25*((xe2-xe1)*(ye3-ye1) - (xe3-xe1)*(ye2-ye1));
    if(J<0) {
      stringstream ss;
      ss << "Negative J found at element " << e << "\n";
      LIBP_ABORT(ss.str())
    }

    /* face 1 */
    dlong base = Nsgeo*Nfaces*e;
    dfloat nx1 = ye2-ye1;
    dfloat ny1 = -(xe2-xe1);
    dfloat  d1 = norm2(nx1,ny1);

    sgeo[base+NXID] = nx1/d1;
    sgeo[base+NYID] = ny1/d1;
    sgeo[base+SJID] = d1/2.;
    sgeo[base+IJID] = 1./J;

    /* face 2 */
    base += Nsgeo;
    dfloat nx2 = ye3-ye2;
    dfloat ny2 = -(xe3-xe2);
    dfloat  d2 = norm2(nx2,ny2);

    sgeo[base+NXID] = nx2/d2;
    sgeo[base+NYID] = ny2/d2;
    sgeo[base+SJID] = d2/2.; // TW fixed bug d1=>d2
    sgeo[base+IJID] = 1./J;

    /* face 3 */
    base += Nsgeo;
    dfloat nx3 = ye1-ye3;
    dfloat ny3 = -(xe1-xe3);
    dfloat  d3 = norm2(nx3,ny3);

    sgeo[base+NXID] = nx3/d3;
    sgeo[base+NYID] = ny3/d3;
    sgeo[base+SJID] = d3/2.;
    sgeo[base+IJID] = 1./J;
  }


  dfloat href = 0.;
  dfloat tol  = 1.;
  for(dlong e=0;e<Nelements;++e){ /* for each non-halo element */
    for(int f=0;f<Nfaces;++f){
      dlong baseM = e*Nfaces + f;

      // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)  A = L*h/2 => (J*2) = (sJ*2)*h/2 => h  = 2*J/sJ
      dfloat hinvM = sgeo[baseM*Nsgeo + SJID]*sgeo[baseM*Nsgeo + IJID];

      href = mymax(hinvM,href);
    }
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

      // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)  A = L*h/2 => (J*2) = (sJ*2)*h/2 => h  = 2*J/sJ
      dfloat hinvM = sgeo[baseM*Nsgeo + SJID]*sgeo[baseM*Nsgeo + IJID];
      dfloat hinvP = sgeo[baseP*Nsgeo + SJID]*sgeo[baseP*Nsgeo + IJID];

      sgeo[baseM*Nsgeo+IHID] = mymax(hinvM,hinvP);
      sgeo[baseP*Nsgeo+IHID] = mymax(hinvM,hinvP);

      if (EToB[f+e*Nfaces] > 0) { //enforce a stronger penalty on boundaries
        sgeo[baseM*Nsgeo+IHID] = mymax(sgeo[baseM*Nsgeo+IHID],tol*href);
        sgeo[baseP*Nsgeo+IHID] = mymax(sgeo[baseP*Nsgeo+IHID],tol*href);
      }
#if 0
      printf("e=%d f=%d (eP=%d,fP=%d) nx=%5.4f, ny=%5.4f, sJ=%5.4f, invJ=%5.4f, hinv=%f\n"
	     ,e,f,eP,fP,
	     sgeo[baseM*Nsgeo+NXID],
	     sgeo[baseM*Nsgeo+NYID],
	     sgeo[baseM*Nsgeo+SJID],
	     sgeo[baseM*Nsgeo+IJID],
	     sgeo[baseM*Nsgeo+IHID]);
#endif
    }
  }

}
