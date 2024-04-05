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

void mesh_t::GeometricFactorsLine1D(){

  /*Set offsets*/
  Nvgeo = 4;

  RXID  = 0;
  JID   = 1;
  JWID  = 2;
  IJWID = 3;

  props["defines/" "p_Nvgeo"]= Nvgeo;
  props["defines/" "p_RXID"]= RXID;
  props["defines/" "p_JID"]= JID;
  props["defines/" "p_JWID"]= JWID;
  props["defines/" "p_IJWID"]= IJWID;
  
  /* unified storage array for geometric factors */
  /* note that we have volume geometric factors for each node */
  vgeo.malloc((Nelements+totalHaloPairs)*Nvgeo*Np);
  
  Nggeo = 1;

  G00ID=0;

  props["defines/" "p_Nggeo"]= Nggeo;
  props["defines/" "p_G00ID"]= G00ID;

  /* number of second order geometric factors */
  ggeo.malloc(Nelements*Nggeo*Np);

  wJ.malloc(Nelements*Np);

  #pragma omp parallel for
  for(dlong e=0;e<Nelements;++e){ /* for each element */
    for(int i=0;i<Nq;++i){
      
      int n = i;
      
      //differentiate physical coordinates
      dfloat xr = 0.0;
      for(int m=0;m<Nq;++m){
	int idr = e*Np + m;
	xr += D[i*Nq+m]*x[idr];
      }
      
      /* compute geometric factors for affine coordinate transform*/
      dfloat J = xr;
      
      LIBP_ABORT("Negative J found at element " << e,
		 J<1e-8);
      
      dfloat rx =  1/xr;
      dfloat JW = J*gllw[i];
      
      /* store geometric factors */
      vgeo[Nvgeo*Np*e + n + Np*RXID] = rx;
      vgeo[Nvgeo*Np*e + n + Np*JID]  = J;
      vgeo[Nvgeo*Np*e + n + Np*JWID] = JW;
      vgeo[Nvgeo*Np*e + n + Np*IJWID] = 1./JW;
      
      /* store second order geometric factors */
      ggeo[Nggeo*Np*e + n + Np*G00ID] = JW*(rx*rx);
      
      wJ[Np*e + n] = JW;
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
