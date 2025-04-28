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

#include "stress.hpp"

// TW need to replace with stress and non-constant viscosity

void stress_t::BuildOperatorDiagonal(memory<dfloat>& diagA){

  if(comm_t::world().rank()==0) {printf("Building diagonal...");fflush(stdout);}

  // assume C0
  memory<dfloat> diagAL(mesh.Np*mesh.Nelements*Nfields);
  
  switch(mesh.elementType){
  case Mesh::TRIANGLES:
    BuildOperatorDiagonalContinuousTri2D(diagAL);
    break;
  case Mesh::QUADRILATERALS:
    BuildOperatorDiagonalContinuousQuad2D(diagAL);
    break;
  case Mesh::TETRAHEDRA:
    BuildOperatorDiagonalContinuousTet3D(diagAL);
    break;
  case Mesh::HEXAHEDRA:
    BuildOperatorDiagonalContinuousHex3D(diagAL);
    break;
  }
  
  //gather the diagonal to assemble it
  ogsMasked.Gather(diagA, diagAL, Nfields, ogs::Add, ogs::Trans);

  if(comm_t::world().rank()==0) printf("done.\n");
}


void stress_t::BuildOperatorDiagonalContinuousTri2D(memory<dfloat>& A) {

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    dlong gbase = eM*mesh.Nggeo;
    dfloat Grr = mesh.ggeo[gbase + mesh.G00ID];
    dfloat Grs = mesh.ggeo[gbase + mesh.G01ID];
    dfloat Gss = mesh.ggeo[gbase + mesh.G11ID];
    dfloat J   = mesh.wJ[eM];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      if (mapB[n+eM*mesh.Np]!=1) { //dont fill rows for masked nodes
        A[eM*mesh.Np+n] = J*lambda*mesh.MM[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Grr*mesh.Srr[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Grs*mesh.Srs[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Gss*mesh.Sss[n+n*mesh.Np];
      } else {
        A[eM*mesh.Np+n] = 1; //just put a 1 so A is invertable
      }
    }

    //add the rank boost for the allNeumann Poisson problem
    if (allNeumann) {
      for(int n=0;n<mesh.Np;++n){
        if (mapB[n+eM*mesh.Np]!=1) { //dont fill rows for masked nodes
          A[eM*mesh.Np+n] += allNeumannPenalty*allNeumannScale*allNeumannScale;
        }
      }
    }
  }
}

void stress_t::BuildOperatorDiagonalContinuousQuad2D(memory<dfloat>& A) {

  /*
    \sum_j d_j ( nu (d_j u_i + d_i u_j ))
    => -\sum_j (d_j phi, nu (d_j u_i + d_i u_j ))
    
    [2.*Dx'*nu*Dx*u + Dy'*nu*Dy*u] + Dy'*nu*Dx*v
    Dx'*nu*Dy*u + [Dx'*nu*Dx*v + 2.*Dy'*nu*Dy*v]
    
  */

  
  for(dlong e=0;e<mesh.Nelements;++e){
    for (int m=0;m<mesh.Nq;++m) {
      for (int n=0;n<mesh.Nq;++n) {
        dlong iid = n+m*mesh.Nq;
	dlong lid = iid + e*mesh.Np;
	dlong uid = 2*lid + 0;
	dlong vid = 2*lid + 1;
	
        if (mapB[n+m*mesh.Nq+e*mesh.Np]!=1) {

	  dlong vbase = e*mesh.Np*mesh.Nvgeo;
          A[uid] = 0;
	  A[vid] = 0;

          for (int k=0;k<mesh.Nq;k++) {
            int id = k+m*mesh.Nq;
	    dfloat Dkn = mesh.D[n+k*mesh.Nq];
	    dfloat rx = mesh.vgeo[vbase + id + mesh.RXID*mesh.Np];
	    dfloat ry = mesh.vgeo[vbase + id + mesh.RYID*mesh.Np];
	    dfloat wJ = mesh.vgeo[vbase + id + mesh.JWID*mesh.Np];	    
	    dfloat nut_km = nut[e*mesh.Np+id];

	    dfloat uGrr = (2.*rx*rx + ry*ry)*nut_km*wJ;
            A[uid] += uGrr*Dkn*Dkn; // strided for gather
	    
	    dfloat vGrr = (2.*ry*ry + rx*rx)*nut_km*wJ;
	    A[vid] += vGrr*Dkn*Dkn;
          }


          for (int k=0;k<mesh.Nq;k++) {
            int id = n+k*mesh.Nq;
	    dfloat Dkm = mesh.D[m+k*mesh.Nq];
	    dfloat sx = mesh.vgeo[vbase + id + mesh.SXID*mesh.Np];
	    dfloat sy = mesh.vgeo[vbase + id + mesh.SYID*mesh.Np];
	    dfloat wJ = mesh.vgeo[vbase + id + mesh.JWID*mesh.Np];	    
	    dfloat nut_km = nut[e*mesh.Np+id];

	    dfloat uGss = (2.*sx*sx + sy*sy)*nut_km*wJ;
            A[uid] += uGss*Dkm*Dkm; // strided for gather

	    dfloat vGss = (2.*sy*sy + sx*sx)*nut_km*wJ;
	    A[vid] += vGss*Dkm*Dkm;
          }

	  {
            int id = n+m*mesh.Nq;
	    dfloat Dnn = mesh.D[n+n*mesh.Nq];
	    dfloat Dmm = mesh.D[m+m*mesh.Nq];
	    dfloat rx = mesh.vgeo[vbase + id + mesh.RXID*mesh.Np];
	    dfloat ry = mesh.vgeo[vbase + id + mesh.RYID*mesh.Np];
	    dfloat sx = mesh.vgeo[vbase + id + mesh.SXID*mesh.Np];
	    dfloat sy = mesh.vgeo[vbase + id + mesh.SYID*mesh.Np];
	    dfloat wJ = mesh.vgeo[vbase + id + mesh.JWID*mesh.Np];	    
	    dfloat nut_nm = nut[e*mesh.Np+id];
	    
	    dfloat uGrs = 2.*(2.*rx*sx + ry*sy)*nut_nm*wJ;
            A[uid] += uGrs*Dnn*Dmm; // strided for gather
	    
	    dfloat vGrs = 2.*(2.*ry*sy + rx*sx)*nut_nm*wJ;
	    A[vid] += vGrs*Dnn*Dmm;

	    // do not need off diagonal blocks
	  }
	  
	  dfloat JW = mesh.vgeo[vbase + n + m*mesh.Nq + mesh.JWID*mesh.Np];
          A[uid] += JW*lambda;
	  A[vid] += JW*lambda;

        } else {
          A[uid] = 1; //just put a 1 so A is invertable
	  A[vid] = 1;
        }
      }
    }

    //add the rank boost for the allNeumann Poisson problem
    if (allNeumann) {
      for(int n=0;n<mesh.Np;++n){
        if (mapB[n+e*mesh.Np]!=1) { //dont fill rows for masked nodes
	  dlong id = e*mesh.Np+n;
	  dfloat fac = allNeumannPenalty*allNeumannScale*allNeumannScale;
          A[Nfields*id+0] += fac;
	  A[Nfields*id+1] += fac;
        }
      }
    }
  }
}


void stress_t::BuildOperatorDiagonalContinuousTet3D(memory<dfloat>& A) {

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    dlong gbase = eM*mesh.Nggeo;
    dfloat Grr = mesh.ggeo[gbase + mesh.G00ID];
    dfloat Grs = mesh.ggeo[gbase + mesh.G01ID];
    dfloat Grt = mesh.ggeo[gbase + mesh.G02ID];
    dfloat Gss = mesh.ggeo[gbase + mesh.G11ID];
    dfloat Gst = mesh.ggeo[gbase + mesh.G12ID];
    dfloat Gtt = mesh.ggeo[gbase + mesh.G22ID];
    dfloat J   = mesh.wJ[eM];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      if (mapB[n+eM*mesh.Np]!=1) { //dont fill rows for masked nodes
        A[eM*mesh.Np+n] = J*lambda*mesh.MM[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Grr*mesh.Srr[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Grs*mesh.Srs[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Grt*mesh.Srt[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Gss*mesh.Sss[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Gst*mesh.Sst[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Gtt*mesh.Stt[n+n*mesh.Np];
      } else {
        A[eM*mesh.Np+n] = 1; //just put a 1 so A is invertable
      }
    }

    //add the rank boost for the allNeumann Poisson problem
    if (allNeumann) {
      for(int n=0;n<mesh.Np;++n){
        if (mapB[n+eM*mesh.Np]!=1) { //dont fill rows for masked nodes
          A[eM*mesh.Np+n] += allNeumannPenalty*allNeumannScale*allNeumannScale;
        }
      }
    }
  }
}

void stress_t::BuildOperatorDiagonalContinuousHex3D(memory<dfloat>& A) {

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    for (int nz=0;nz<mesh.Nq;nz++) {
    for (int ny=0;ny<mesh.Nq;ny++) {
    for (int nx=0;nx<mesh.Nq;nx++) {
      int idn = nx+ny*mesh.Nq+nz*mesh.Nq*mesh.Nq;
      if (mapB[idn+eM*mesh.Np]!=1) {
        A[eM*mesh.Np+idn] = 0;

        int id = nx+ny*mesh.Nq+nz*mesh.Nq*mesh.Nq;
        dlong base = eM*mesh.Np*mesh.Nggeo;


        dfloat Grs = mesh.ggeo[base + id + mesh.G01ID*mesh.Np];
        A[eM*mesh.Np+idn] += 2*Grs*mesh.D[nx+nx*mesh.Nq]*mesh.D[ny+ny*mesh.Nq];

        dfloat Grt = mesh.ggeo[base + id + mesh.G02ID*mesh.Np];
        A[eM*mesh.Np+idn] += 2*Grt*mesh.D[nx+nx*mesh.Nq]*mesh.D[nz+nz*mesh.Nq];

        dfloat Gst = mesh.ggeo[base + id + mesh.G12ID*mesh.Np];
        A[eM*mesh.Np+idn] += 2*Gst*mesh.D[ny+ny*mesh.Nq]*mesh.D[nz+nz*mesh.Nq];

        for (int k=0;k<mesh.Nq;k++) {
          int iid = k+ny*mesh.Nq+nz*mesh.Nq*mesh.Nq;
          dfloat Grr = mesh.ggeo[base + iid + mesh.G00ID*mesh.Np];
          A[eM*mesh.Np+idn] += Grr*mesh.D[nx+k*mesh.Nq]*mesh.D[nx+k*mesh.Nq];
        }

        for (int k=0;k<mesh.Nq;k++) {
          int iid = nx+k*mesh.Nq+nz*mesh.Nq*mesh.Nq;
          dfloat Gss = mesh.ggeo[base + iid + mesh.G11ID*mesh.Np];
          A[eM*mesh.Np+idn] += Gss*mesh.D[ny+k*mesh.Nq]*mesh.D[ny+k*mesh.Nq];
        }

        for (int k=0;k<mesh.Nq;k++) {
          int iid = nx+ny*mesh.Nq+k*mesh.Nq*mesh.Nq;
          dfloat Gtt = mesh.ggeo[base + iid + mesh.G22ID*mesh.Np];
          A[eM*mesh.Np+idn] += Gtt*mesh.D[nz+k*mesh.Nq]*mesh.D[nz+k*mesh.Nq];
        }

        dfloat JW = mesh.wJ[eM*mesh.Np + idn];
        A[eM*mesh.Np+idn] += JW*lambda;
      } else {
        A[eM*mesh.Np+idn] = 1; //just put a 1 so A is invertable
      }
    }
    }
    }

    //add the rank boost for the allNeumann Poisson problem
    if (allNeumann) {
      for(int n=0;n<mesh.Np;++n){
        if (mapB[n+eM*mesh.Np]!=1) { //dont fill rows for masked nodes
          A[eM*mesh.Np+n] += allNeumannPenalty*allNeumannScale*allNeumannScale;
        }
      }
    }
  }
}
