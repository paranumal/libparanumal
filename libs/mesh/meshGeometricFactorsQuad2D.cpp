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

void mesh_t::GeometricFactorsQuad2D(){

  /*Set offsets*/
  Nvgeo = 8;

  RXID  = 0;
  RYID  = 1;
  SXID  = 2;
  SYID  = 3;
  JID   = 4;
  JWID  = 5;
  IJWID = 6;
  VOLHID = 7;
  
  props["defines/" "p_Nvgeo"]= Nvgeo;
  props["defines/" "p_RXID"]= RXID;
  props["defines/" "p_SXID"]= SXID;
  props["defines/" "p_RYID"]= RYID;
  props["defines/" "p_SYID"]= SYID;
  props["defines/" "p_JID"]= JID;
  props["defines/" "p_JWID"]= JWID;
  props["defines/" "p_IJWID"]= IJWID;
  props["defines/" "p_VOLHID"]= VOLHID;

  /* unified storage array for geometric factors */
  /* note that we have volume geometric factors for each node */
  vgeo.malloc((Nelements+totalHaloPairs)*Nvgeo*Np);

  Nggeo = 3;

  G00ID=0;
  G01ID=1;
  G11ID=2;

  props["defines/" "p_Nggeo"]= Nggeo;
  props["defines/" "p_G00ID"]= G00ID;
  props["defines/" "p_G01ID"]= G01ID;
  props["defines/" "p_G11ID"]= G11ID;

  /* number of second order geometric factors */
  ggeo.malloc(Nelements*Nggeo*Np);

  wJ.malloc(Nelements*Np);

  #pragma omp parallel for
  for(dlong e=0;e<Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*Nverts+0;

    dfloat xe1 = EX[id+0];
    dfloat xe2 = EX[id+1];
    dfloat xe3 = EX[id+2];
    dfloat xe4 = EX[id+3];

    dfloat ye1 = EY[id+0];
    dfloat ye2 = EY[id+1];
    dfloat ye3 = EY[id+2];
    dfloat ye4 = EY[id+3];
    
    for(int j=0;j<Nq;++j){
      for(int i=0;i<Nq;++i){

        int n = i + j*Nq;

        //differentiate physical coordinates
        dfloat xr = 0.0;
        dfloat xs = 0.0;
        dfloat yr = 0.0;
        dfloat ys = 0.0;

        for(int m=0;m<Nq;++m){
          int idr = e*Np + j*Nq + m;
          int ids = e*Np + m*Nq + i;
          xr += D[i*Nq+m]*x[idr];
          xs += D[j*Nq+m]*x[ids];
          yr += D[i*Nq+m]*y[idr];
          ys += D[j*Nq+m]*y[ids];
        }

        /* compute geometric factors for affine coordinate transform*/
        dfloat J = xr*ys - xs*yr;

        LIBP_ABORT("Negative J found at element " << e,
                   J<1e-8);

        dfloat rx =  ys/J;
        dfloat ry = -xs/J;
        dfloat sx = -yr/J;
        dfloat sy =  xr/J;
        dfloat JW = J*gllw[i]*gllw[j];

        /* store geometric factors */
        vgeo[Nvgeo*Np*e + n + Np*RXID] = rx;
        vgeo[Nvgeo*Np*e + n + Np*RYID] = ry;
        vgeo[Nvgeo*Np*e + n + Np*SXID] = sx;
        vgeo[Nvgeo*Np*e + n + Np*SYID] = sy;
        vgeo[Nvgeo*Np*e + n + Np*JID]  = J;
        vgeo[Nvgeo*Np*e + n + Np*JWID] = JW;
        vgeo[Nvgeo*Np*e + n + Np*IJWID] = 1./JW;

	/* compute a penalty scaling */
	dfloat h1 = sqrt( (xe2-xe1)*(xe2-xe1)  + (ye2-ye1)*(ye2-ye1) );
	dfloat h2 = sqrt( (xe3-xe2)*(xe3-xe2)  + (ye3-ye2)*(ye3-ye2) );
	dfloat h3 = sqrt( (xe4-xe3)*(xe4-xe3)  + (ye4-ye3)*(ye4-ye3) );
	dfloat h4 = sqrt( (xe1-xe4)*(xe1-xe4)  + (ye1-ye4)*(ye1-ye4) );
	
	dfloat h = dim*(4.*J)/(h1+h2+h3+h4); // ratio of area to perimeter (assumes bilinear geo map)
	vgeo[Nvgeo*Np*e +  n+ Np*VOLHID] = h;
	
        /* store second order geometric factors */
        ggeo[Nggeo*Np*e + n + Np*G00ID] = JW*(rx*rx + ry*ry);
        ggeo[Nggeo*Np*e + n + Np*G01ID] = JW*(rx*sx + ry*sy);
        ggeo[Nggeo*Np*e + n + Np*G11ID] = JW*(sx*sx + sy*sy);

        wJ[Np*e + n] = JW;
      }
    }
  }

  halo.Exchange(vgeo, Nvgeo*Np);

  o_wJ   = platform.malloc<dfloat>(wJ);
  o_vgeo = platform.malloc<dfloat>(vgeo);
  o_ggeo = platform.malloc<dfloat>(ggeo);

  if constexpr (std::is_same_v<dfloat,pfloat>) {
    o_pfloat_ggeo = o_ggeo;
    o_pfloat_vgeo = o_vgeo;
    o_pfloat_wJ   = o_wJ;
  } else {
    memory<pfloat> pfloat_wJ(Nelements*Np);
    memory<pfloat> pfloat_ggeo(Nggeo*Nelements*Np);
    memory<pfloat> pfloat_vgeo(Nvgeo*Nelements*Np);

    for(int n=0;n<Nggeo*Nelements*Np;++n)
      pfloat_ggeo[n] = ggeo[n];
    for(int n=0;n<Nvgeo*Nelements*Np;++n)
      pfloat_vgeo[n] = vgeo[n];
    for(int n=0;n<Nelements*Np;++n)
      pfloat_wJ[n] = wJ[n];

    o_pfloat_ggeo = platform.malloc<pfloat>(pfloat_ggeo);
    o_pfloat_vgeo = platform.malloc<pfloat>(pfloat_vgeo);
    o_pfloat_wJ   = platform.malloc<pfloat>(pfloat_wJ);
  }
}

} //namespace libp
