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
  
  std::cout << "BuildOperatorDiagonalContinuousTri2D not implemented" << std::endl;
  exit(-1);
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

	dfloat fac = 2; 
	
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

	    dfloat uGrr = (fac*rx*rx + ry*ry)*nut_km*wJ;
            A[uid] += uGrr*Dkn*Dkn; // strided for gather
	    
	    dfloat vGrr = (fac*ry*ry + rx*rx)*nut_km*wJ;
	    A[vid] += vGrr*Dkn*Dkn;
          }


          for (int k=0;k<mesh.Nq;k++) {
            int id = n+k*mesh.Nq;
	    dfloat Dkm = mesh.D[m+k*mesh.Nq];
	    dfloat sx = mesh.vgeo[vbase + id + mesh.SXID*mesh.Np];
	    dfloat sy = mesh.vgeo[vbase + id + mesh.SYID*mesh.Np];
	    dfloat wJ = mesh.vgeo[vbase + id + mesh.JWID*mesh.Np];	    
	    dfloat nut_nk = nut[e*mesh.Np+id];
	    
	    dfloat uGss = (fac*sx*sx + sy*sy)*nut_nk*wJ;
            A[uid] += uGss*Dkm*Dkm; // strided for gather
	    
	    dfloat vGss = (fac*sy*sy + sx*sx)*nut_nk*wJ;
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
	    
	    dfloat uGrs = 2.*(fac*rx*sx + ry*sy)*nut_nm*wJ;
            A[uid] += uGrs*Dnn*Dmm; // strided for gather
	    
	    dfloat vGrs = 2.*(fac*ry*sy + rx*sx)*nut_nm*wJ;
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

  std::cout << "BuildOperatorDiagonalContinuousTet3D not implemented" << std::endl;
  exit(-1);
}

void stress_t::BuildOperatorDiagonalContinuousHex3D(memory<dfloat>& A) {

  std::cout << "BuildOperatorDiagonalContinuousHex3D not implemented" << std::endl;
  exit(-1);
}
