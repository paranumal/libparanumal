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

#include "esdg.hpp"

void esdg_t::SetupQuad2D(){

  // set up cubature grid
  cubature   = (settings.compareSetting("ADVECTION TYPE", "CUBATURE")) ? 1:0;
  entropyStable  = (settings.compareSetting("ADVECTION TYPE", "ENTROPYSTABLE")) ? 1:0;

  mesh.CubatureSetup();
  mesh.CubaturePhysicalNodes();

  // TURNED OFF FOR A MOMENT
  Nfields   = mesh.dim + 2;
  Ngrads = mesh.dim*mesh.dim; // fix this
  
  // build element centers
  memory<dfloat> cx, cy, cz;
  cx.malloc(mesh.Nelements);
  cy.malloc(mesh.Nelements);
  cz.malloc(mesh.Nelements);

  for(dlong e=0;e<mesh.Nelements;++e){
    dfloat cxe = 0, cye = 0, cze = 0;
    for(int n=0;n<mesh.Nverts;++n){
      cxe += mesh.EX[e*mesh.Nverts+n];
      cye += mesh.EY[e*mesh.Nverts+n];
    }
    cx[e] = cxe/mesh.Nverts;
    cy[e] = cye/mesh.Nverts;
    cz[e] = 0;
  }
  o_cx = platform.malloc<dfloat>(mesh.Nelements, cx);
  o_cy = platform.malloc<dfloat>(mesh.Nelements, cy);
  o_cz = platform.malloc<dfloat>(mesh.Nelements, cz);


  // build ESDG stuff specific to quads

  // face quadrature
  int Nq = mesh.Nq;
  int cubNq = mesh.cubNq;
  
  // to match Jesse Chan code we use: GQ(N) x GQ(N) (i.e. (N+1)x(N+1) Gauss-Legendre node)

  memory<dfloat> rout1D(2);
  rout1D[0] = -1;
  rout1D[1] = +1;
  
  memory<dfloat>  r1D(Nq);
  memory<dfloat>  w1D(Nq);
  memory<dfloat>  D1D(Nq*Nq);
  memory<dfloat> If1D(Nq*2);
  mesh.JacobiGLL(mesh.N, r1D, w1D); //Gauss-Legendre-Lobatto nodes
  mesh.Dmatrix1D(mesh.N, r1D, r1D, D1D);
  mesh.InterpolationMatrix1D(mesh.N, r1D, rout1D, If1D);
  
  // 1D quadrature operators
  memory<dfloat> rq1D(Nq);
  memory<dfloat> wq1D(Nq);
  memory<dfloat> Dq1D(Nq*Nq);

  memory<dfloat> Iqf1D(Nq*2);
  memory<dfloat> Iq1D(Nq*Nq);

  mesh.JacobiGQ (0., 0., mesh.N, rq1D, wq1D); //Gauss-Legendre-Lobatto nodes
  mesh.Dmatrix1D(mesh.N, rq1D, rq1D, Dq1D);
  mesh.InterpolationMatrix1D(mesh.N, rq1D, rout1D, Iqf1D);

  mesh.InterpolationMatrix1D(mesh.N, r1D, rq1D, Iq1D);

  //
  memory<dfloat> Dq1D_skewT(Nq*Nq);
  memory<dfloat> Lqf1DT(Nq*2);
  memory<dfloat> Iqf1DT(Nq*2);

  // differentiate at GQ nodes
  for(int i=0;i<Nq;++i){
    for(int j=0;j<Nq;++j){
      // store in column major
      Dq1D_skewT[i + j*Nq] = (wq1D[i]*Dq1D[i*Nq+j] - Dq1D[j*Nq+i]*wq1D[j])/wq1D[i];
    }
  }

  for(int i=0;i<Nq;++i){
    for(int j=0;j<2;++j){
      // lift from end points to GQ nodes
      Lqf1DT[i+j*Nq] = (Iqf1D[j*Nq + i])/wq1D[i];
      // interpolate from GQ to end points
      Iqf1DT[i+j*Nq] = Iqf1D[j*Nq+i];
    }
  }

  o_esDq1D_skew = platform.malloc<dfloat>(Nq*Nq, Dq1D_skewT);
  o_esLqf1D = platform.malloc<dfloat>(Nq*2, Lqf1DT);
  o_esIqf1D = platform.malloc<dfloat>(Nq*2, Iqf1DT);
  
  int Nq1 = mesh.intNfp;
  int N1 = Nq1-1;
  memory<dfloat> ir, is, it, iw;
  ir.malloc(Nq1*mesh.Nfaces);
  is.malloc(Nq1*mesh.Nfaces);
  it.malloc(Nq1*mesh.Nfaces);
  iw.malloc(Nq1*mesh.Nfaces);

  for(int n=0;n<Nq1;++n){
    ir[0*Nq1 + n] =  mesh.intr[n];
    ir[1*Nq1 + n] = +1.0;
    ir[2*Nq1 + n] = -mesh.intr[n]; // sign
    ir[3*Nq1 + n] = -1.0;
    
    is[0*Nq1 + n] = -1.0;
    is[1*Nq1 + n] =  mesh.intr[n];
    is[2*Nq1 + n] = +1.0;
    is[3*Nq1 + n] = -mesh.intr[n];
    
    iw[0*Nq1 + n] =  mesh.intw[n];
    iw[1*Nq1 + n] =  mesh.intw[n];
    iw[2*Nq1 + n] =  mesh.intw[n];
    iw[3*Nq1 + n] =  mesh.intw[n];
  }

  // geometric factors (maybe multiplied by J) at volume GQ nodes
  // esVgeo <= vgeo*J
  memory<dfloat> esVgeo(mesh.Nvgeo*mesh.Np*mesh.Nelements);
  memory<dfloat> tmp(Nq*Nq);

  for(dlong e=0;e<mesh.Nelements;++e){

    for(int g=0;g<mesh.Nvgeo;++g){

      for(int i=0;i<Nq;++i){
	for(int j=0;j<Nq;++j){
	  tmp[i*Nq+j] = 0;
	  for(int n=0;n<Nq;++n){
	    dlong  id = e*mesh.Nvgeo*Nq*Nq +         g*Nq*Nq + i*Nq + n;
	    dlong Jid = e*mesh.Nvgeo*Nq*Nq + mesh.JID*Nq*Nq + i*Nq + n;
	    tmp[i*Nq+j] += Iq1D[j*Nq+n]*mesh.vgeo[id]*mesh.vgeo[Jid];
	  }
	}
      }

      for(int i=0;i<Nq;++i){
	for(int j=0;j<Nq;++j){
	  dfloat Itmp = 0;
	  for(int n=0;n<Nq;++n){
	    Itmp += Iq1D[i*Nq+n]*tmp[n*Nq+j];
	  }
	  dlong  id = e*mesh.Nvgeo*Nq*Nq +         g*Nq*Nq + i*Nq + j;
	  esVgeo[id] = Itmp;
	}
      }
    }
  }

  o_esVgeo = platform.malloc<dfloat>(mesh.Nvgeo*mesh.Np*mesh.Nelements, esVgeo);

#if 0
  // 1D operators used for quadrature
  o_esMapPq = ; // positive trace index of nodes in surface GQ nodes
  o_esVfgeo = ; // geometric factors (maybe multiplied by J) at surface GQ nodes
  o_esSgeo = ;  // normal and surface Jacobian at surface GQ nodes
#endif
  
  // to do interpolate to GLL for output
  // to do evaluate initial 
  
}

