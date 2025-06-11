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
  
  dfloat volumeTet3D(const memory<dfloat> &EX,
		     const memory<dfloat> &EY,
		     const memory<dfloat> &EZ) {

    // Form edge vectors from vertex 0 to vertices 1, 2, and 3:
    dfloat v1x = EX[1] - EX[0], v1y = EY[1] - EY[0], v1z = EZ[1] - EZ[0];
    dfloat v2x = EX[2] - EX[0], v2y = EY[2] - EY[0], v2z = EZ[2] - EZ[0];
    dfloat v3x = EX[3] - EX[0], v3y = EY[3] - EY[0], v3z = EZ[3] - EZ[0];
    
    // Compute cross product v2 × v3:
    dfloat cx = v2y*v3z - v2z*v3y;
    dfloat cy = v2z*v3x - v2x*v3z;
    dfloat cz = v2x*v3y - v2y*v3x;
    
    // Scalar triple product v1 · (v2 × v3)
    dfloat triple = v1x*cx + v1y*cy + v1z*cz;
    
    // Volume = |triple| / 6
    return std::fabs(triple) / (dfloat)6.0;
  }

  dfloat faceAreaTet3D(const memory<int> &faceVertices, int face,
		       const memory<dfloat> &EX,
		       const memory<dfloat> &EY,
		       const memory<dfloat> &EZ){

    int i = faceVertices[face*3+0];
    int j = faceVertices[face*3+1];
    int k = faceVertices[face*3+2];

    /* edge vectors v1 = Pj – Pi, v2 = Pk – Pi */
    dfloat v1x = EX[j] - EX[i], v1y = EY[j] - EY[i], v1z = EZ[j] - EZ[i];
    dfloat v2x = EX[k] - EX[i], v2y = EY[k] - EY[i], v2z = EZ[k] - EZ[i];

    /* cross product v1 × v2 */
    dfloat cx = v1y*v2z - v1z*v2y;
    dfloat cy = v1z*v2x - v1x*v2z;
    dfloat cz = v1x*v2y - v1y*v2x;

    /* triangle area = 0.5 * ||cross|| */
    return (dfloat)0.5 * sqrt(cx*cx + cy*cy + cz*cz);
}

  void mesh_t::GeometricFactorsTet3D(){

  /*Set offsets*/
  Nvgeo = 11;

  RXID  = 0;
  RYID  = 1;
  RZID  = 2;
  SXID  = 3;
  SYID  = 4;
  SZID  = 5;
  TXID  = 6;
  TYID  = 7;
  TZID  = 8;
  JID   = 9;
  VOLHID = 10;
  
  props["defines/" "p_Nvgeo"]= Nvgeo;
  props["defines/" "p_RXID"]= RXID;
  props["defines/" "p_SXID"]= SXID;
  props["defines/" "p_TXID"]= TXID;

  props["defines/" "p_RYID"]= RYID;
  props["defines/" "p_SYID"]= SYID;
  props["defines/" "p_TYID"]= TYID;

  props["defines/" "p_RZID"]= RZID;
  props["defines/" "p_SZID"]= SZID;
  props["defines/" "p_TZID"]= TZID;

  props["defines/" "p_JID"]= JID;
  props["defines/" "p_VOLHID"]= VOLHID;

  /* unified storage array for geometric factors */
  vgeo.malloc((Nelements+totalHaloPairs)*Nvgeo);

  Nggeo = 6;

  G00ID=0;
  G01ID=1;
  G02ID=2;
  G11ID=3;
  G12ID=4;
  G22ID=5;

  props["defines/" "p_Nggeo"]= Nggeo;
  props["defines/" "p_G00ID"]= G00ID;
  props["defines/" "p_G01ID"]= G01ID;
  props["defines/" "p_G02ID"]= G02ID;
  props["defines/" "p_G11ID"]= G11ID;
  props["defines/" "p_G12ID"]= G12ID;
  props["defines/" "p_G22ID"]= G22ID;

  /* number of second order geometric factors */
  ggeo.malloc(Nelements*Nggeo);

  wJ.malloc(Nelements);
  
  dfloat maxvolh = 0;
  // dfloat minJ = 1e9, maxJ = -1e9;

  #pragma omp parallel for reduction(max:maxvolh)
  for(dlong e=0;e<Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*Nverts;

    /* vertex coordinates */
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

    LIBP_ABORT("Negative J found at element " << e, J<0);

    // minJ = mymin(minJ,J);
    // maxJ = mymax(maxJ,J);

    /* store geometric factors */
    vgeo[Nvgeo*e + RXID] = rx;
    vgeo[Nvgeo*e + RYID] = ry;
    vgeo[Nvgeo*e + RZID] = rz;
    vgeo[Nvgeo*e + SXID] = sx;
    vgeo[Nvgeo*e + SYID] = sy;
    vgeo[Nvgeo*e + SZID] = sz;
    vgeo[Nvgeo*e + TXID] = tx;
    vgeo[Nvgeo*e + TYID] = ty;
    vgeo[Nvgeo*e + TZID] = tz;
    vgeo[Nvgeo*e +  JID] = J;
    //    printf("geo: %g,%g,%g - %g,%g,%g - %g,%g,%g\n",
    //     rx,ry,rz, sx,sy,sz, tx,ty,tz);

    dfloat A[4], area = 0;
    dfloat minh = 1e9;
    dfloat vol = volumeTet3D(EX+id,EY+id,EZ+id);
    for(int f=0;f<Nverts;++f){
      A[f] = faceAreaTet3D(faceVertices, f, EX+id,EY+id,EZ+id);
      area += A[f];
      minh = std::min(minh, dim*vol/A[f]);
    }

    dfloat h = dim*vol/area;
    //    vgeo[Nvgeo*e + VOLHID] = h;
    vgeo[Nvgeo*e + VOLHID] = minh; 
    maxvolh = std::max(maxvolh, h);
    
    /* store second order geometric factors */
    ggeo[Nggeo*e + G00ID] = J*(rx*rx + ry*ry + rz*rz);
    ggeo[Nggeo*e + G01ID] = J*(rx*sx + ry*sy + rz*sz);
    ggeo[Nggeo*e + G02ID] = J*(rx*tx + ry*ty + rz*tz);
    ggeo[Nggeo*e + G11ID] = J*(sx*sx + sy*sy + sz*sz);
    ggeo[Nggeo*e + G12ID] = J*(sx*tx + sy*ty + sz*tz);
    ggeo[Nggeo*e + G22ID] = J*(tx*tx + ty*ty + tz*tz);

    wJ[e] = J;
  }
  //printf("minJ = %g, maxJ = %g\n", minJ, maxJ);

  std::cout << "maxvolh=" << maxvolh << std::endl;
  
  halo.Exchange(vgeo, Nvgeo);

  o_wJ   = platform.malloc<dfloat>(wJ);
  o_vgeo = platform.malloc<dfloat>(vgeo);
  o_ggeo = platform.malloc<dfloat>(ggeo);

  if constexpr (std::is_same_v<dfloat,pfloat>) {
    o_pfloat_ggeo = o_ggeo;
    o_pfloat_vgeo = o_vgeo;
    o_pfloat_wJ   = o_wJ;
  } else {
    memory<pfloat> pfloat_wJ(Nelements);
    memory<pfloat> pfloat_ggeo(Nggeo*Nelements);
    memory<pfloat> pfloat_vgeo(Nvgeo*Nelements);

    for(int n=0;n<Nggeo*Nelements;++n)
      pfloat_ggeo[n] = ggeo[n];
    for(int n=0;n<Nvgeo*Nelements;++n)
      pfloat_vgeo[n] = vgeo[n];
    for(int n=0;n<Nelements;++n)
      pfloat_wJ[n] = wJ[n];

    o_pfloat_ggeo = platform.malloc<pfloat>(pfloat_ggeo);
    o_pfloat_vgeo = platform.malloc<pfloat>(pfloat_vgeo);
    o_pfloat_wJ   = platform.malloc<pfloat>(pfloat_wJ);
  }


}

} //namespace libp
