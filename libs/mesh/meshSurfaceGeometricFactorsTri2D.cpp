/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

namespace libp {

void mesh_t::SurfaceGeometricFactorsTri2D(){

  /* unified storage array for geometric factors */
  Nsgeo = 5;

  NXID  = 0;
  NYID  = 1;
  SJID  = 2;
  IJID  = 3;
  IHID  = 4;

  props["defines/" "p_Nsgeo"]= Nsgeo;
  props["defines/" "p_NXID"]= NXID;
  props["defines/" "p_NYID"]= NYID;
  props["defines/" "p_SJID"]= SJID;
  props["defines/" "p_IJID"]= IJID;
  props["defines/" "p_IHID"]= IHID;

  sgeo.malloc(Nelements*Nsgeo*Nfaces);

  memory<dfloat> hinv((Nelements+totalHaloPairs)*Nfaces);

  for(dlong e=0;e<Nelements;++e){ /* for each element */

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

    LIBP_ABORT("Negative J found at element " << e,
               J<0);

    /* face 1 */
    dlong base = Nsgeo*Nfaces*e;
    dfloat nx1 = ye2-ye1;
    dfloat ny1 = -(xe2-xe1);
    dfloat  d1 = sqrt((nx1)*(nx1)+(ny1)*(ny1));

    sgeo[base+NXID] = nx1/d1;
    sgeo[base+NYID] = ny1/d1;
    sgeo[base+SJID] = d1/2.;
    sgeo[base+IJID] = 1./J;

    hinv[Nfaces*e+0] = 0.25*d1/J;

    /* face 2 */
    base += Nsgeo;
    dfloat nx2 = ye3-ye2;
    dfloat ny2 = -(xe3-xe2);
    dfloat  d2 = sqrt((nx2)*(nx2)+(ny2)*(ny2));

    sgeo[base+NXID] = nx2/d2;
    sgeo[base+NYID] = ny2/d2;
    sgeo[base+SJID] = d2/2.; // TW fixed bug d1=>d2
    sgeo[base+IJID] = 1./J;

    hinv[Nfaces*e+1] = 0.25*d2/J;

    /* face 3 */
    base += Nsgeo;
    dfloat nx3 = ye1-ye3;
    dfloat ny3 = -(xe1-xe3);
    dfloat  d3 = sqrt((nx3)*(nx3)+(ny3)*(ny3));

    sgeo[base+NXID] = nx3/d3;
    sgeo[base+NYID] = ny3/d3;
    sgeo[base+SJID] = d3/2.;
    sgeo[base+IJID] = 1./J;

    hinv[Nfaces*e+2] = 0.25*d3/J;
  }

  halo.Exchange(hinv, Nfaces);

  // dfloat href = 0.;
  // dfloat tol  = 1.;
  // for(dlong e=0;e<Nelements;++e){ /* for each non-halo element */
  //   for(int f=0;f<Nfaces;++f){
  //     dlong baseM = e*Nfaces + f;

  //     // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)  A = L*h/2 => (J*2) = (sJ*2)*h/2 => h  = 2*J/sJ
  //     dfloat hinvM = sgeo[baseM*Nsgeo + SJID]*sgeo[baseM*Nsgeo + IJID];

  //     href = mymax(hinvM,href);
  //   }
  // }

  for(dlong eM=0;eM<Nelements;++eM){ /* for each non-halo element */
    for(int fM=0;fM<Nfaces;++fM){
      dlong eP = EToE[eM*Nfaces+fM];

      if (eP<0) eP = eM;

      int fP = EToF[eM*Nfaces+fM];
      if (fP<0) fP = fM;

      dlong baseM = eM*Nfaces + fM;
      dlong baseP = eP*Nfaces + fP;

      // rescaling - A = L*h/2 => (J*2) = (sJ*2)*h/2 => h  = 2*J/sJ
      dfloat hinvM = hinv[baseM];
      dfloat hinvP = hinv[baseP];
      sgeo[baseM*Nsgeo+IHID] = std::max(hinvM,hinvP);

      // if (EToB[fM+eM*Nfaces] > 0) { //enforce a stronger penalty on boundaries
      //   sgeo[baseM*Nsgeo+IHID] *= 2;
      // }
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

  o_sgeo = platform.malloc<dfloat>(sgeo);
}

} //namespace libp
