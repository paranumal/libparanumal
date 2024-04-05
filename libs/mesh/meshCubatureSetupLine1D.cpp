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

void mesh_t::CubatureSetupLine1D(){

  /* Quadrature data */
  cubN = N+1;
  cubNq = cubN+1;
  cubNp = cubNq;
  cubNfp = 1;
  intNfp = 1;

  // cubN+1 point Gauss-Legendre quadrature
  JacobiGQ(0, 0, cubN, cubr, cubw);

  // GLL to GL interpolation matrix
  InterpolationMatrix1D(N, gllz, cubr, cubInterp);

  //cubature project cubProject = cubInterp^T
  cubProject.malloc(cubNq*Nq);
  linAlg_t::matrixTranspose(cubNq, Nq, cubInterp, Nq, cubProject, cubNq);

  //cubature derivates matrix, cubD: differentiate on cubature nodes
  Dmatrix1D(cubN, cubr, cubr, cubD);

  // weak cubature derivative cubPDT = cubProject * cubD^T
  CubatureWeakDmatrix1D(Nq, cubNq, cubProject, cubD, cubPDT);

  // add compile time constants to kernels
  props["defines/" "p_cubNq"]= cubNq;
  props["defines/" "p_cubNp"]= cubNp;
  props["defines/" "p_intNfp"]= intNfp;
  props["defines/" "p_intNfpNfaces"]= intNfp*Nfaces;
  props["defines/" "p_cubNfp"]= cubNfp;

  // build transposes (we hold matrices as column major on device)
  memory<dfloat> cubProjectT(cubNq*Nq);
  memory<dfloat> cubInterpT(cubNq*Nq);
  linAlg_t::matrixTranspose(cubNq, Nq, cubInterp, Nq, cubInterpT, cubNq);
  linAlg_t::matrixTranspose(Nq, cubNq, cubProject, cubNq, cubProjectT, Nq);

  memory<dfloat> cubPDTT(cubNq*Nq);
  linAlg_t::matrixTranspose(Nq, cubNq, cubPDT, cubNq, cubPDTT, Nq);

  o_cubInterp  = platform.malloc<dfloat>(Nq*cubNq, cubInterpT);
  o_cubProject = platform.malloc<dfloat>(Nq*cubNq, cubProjectT);

  o_cubPDT = platform.malloc<dfloat>(Nq*cubNq, cubPDTT);
  o_cubD   = platform.malloc<dfloat>(cubNq*cubNq, cubD);

  o_intInterp = o_cubInterp;
  o_intLIFT = o_cubProject;

  cubwJ.malloc(Nelements*cubNp);
  cubvgeo.malloc(Nelements*Nvgeo*cubNp);
  cubggeo.malloc(Nelements*Nggeo*cubNp);
  cubsgeo.malloc(Nelements*Nsgeo*cubNq*Nfaces);

  //temp arrays
  memory<dfloat> xre(Np);
  memory<dfloat> xre1(cubNq);

  //geometric data for quadrature
  for(dlong e=0;e<Nelements;++e){ /* for each element */
    for(int i=0;i<Nq;++i){
      int n = i;
      
      //differentiate physical coordinates
      xre[n] = 0.0;
      
      for(int m=0;m<Nq;++m){
	int idr = e*Np + m;
	xre[n] += D[i*Nq+m]*x[idr];
      }
    }
    
    //interpolate derivaties to cubature
    for(int i=0;i<cubNq;++i){
      xre1[i] = 0.0; 
      for(int n=0;n<Nq;++n){
	xre1[i] += cubInterp[i*Nq + n]*xre[n];
      }
    }

    for(int i=0;i<cubNq;++i){
      dfloat xr = xre1[i];
      
      /* compute geometric factors for affine coordinate transform*/
      dfloat J = xr;
      
      LIBP_ABORT("Negative J found at element " << e, J<1e-8);
      
      dfloat rx =  1./xr;
      dfloat JW = J*cubw[i];
      
      /* store geometric factors */
      dlong base = Nvgeo*cubNp*e + i;
      cubvgeo[base + cubNp*RXID] = rx;
      cubvgeo[base + cubNp*JID]  = J;
      cubvgeo[base + cubNp*JWID] = JW;
      cubvgeo[base + cubNp*IJWID] = 1./JW;
      
      /* store second order geometric factors */
      base = Nggeo*cubNp*e + i;
      cubggeo[base + cubNp*G00ID] = JW*(rx*rx);
      
      cubwJ[cubNp*e + i] = JW;
    }
    
    for(int f=0;f<Nfaces;++f){ // for each face
      for(int m=0;m<1;++m){  // for each node on face
	
        //interpolate derivatives of physical coordinates
        dfloat xr = 0.0;
        for(int n=0;n<Nfp;++n){  // for each node on face
          /* volume index of face node */
          int idn = faceNodes[f*Nfp+n];
          xr += cubInterp[m*Nq + n]*xre[idn];
        }
	
        /* compute geometric factors for affine coordinate transform*/
        dfloat J = xr;
	
        dfloat rx =  1/xr;
	
        /* face f normal and length */
        dfloat nx=0.0;
        switch(f){
        case 0: nx = +rx;  break;
        case 1: nx = -rx;  break;
        }
        dfloat  sJ = sqrt((nx)*(nx));
        nx /= sJ; 
        sJ *= J;

        /* output index */
        dlong base = Nsgeo*(Nfaces*cubNq*e + cubNq*f + m);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        cubsgeo[base+NXID] = nx;
        cubsgeo[base+SJID] = sJ;
        cubsgeo[base+IJID] = 1./J;

        cubsgeo[base+WIJID] = 1./(J*cubw[0]);
        cubsgeo[base+WSJID] = sJ*cubw[m];
      }
    }
  }

  o_cubwJ   = platform.malloc<dfloat>(Nelements*cubNp, cubwJ);
  o_cubvgeo = platform.malloc<dfloat>(Nelements*Nvgeo*cubNp, cubvgeo);
  o_cubggeo = platform.malloc<dfloat>(Nelements*Nggeo*cubNp, cubggeo);
  o_cubsgeo = platform.malloc<dfloat>(Nelements*Nfaces*Nsgeo, cubsgeo);
}
  
} //namespace libp
