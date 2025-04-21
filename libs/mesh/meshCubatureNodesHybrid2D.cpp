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

  void mesh_t::ElementCubaturePhysicalNodesTri2D(dlong e, memory<dfloat> &cubx, memory<dfloat> &cuby){
    
    dlong id = e*maxNverts+0;

    dfloat xe1 = EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = EX[id+1];
    dfloat xe3 = EX[id+2];

    dfloat ye1 = EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = EY[id+1];
    dfloat ye3 = EY[id+2];

    Mesh::ElementType etype = elementTypes[e];
    dlong cubNpe = cubatureNp(etype);
    
    dlong cnt = 0;
    for(int n=0;n<cubNpe;++n){ /* for each node */

      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = cubrHybrid[etype][n];
      dfloat sn = cubsHybrid[etype][n];
      
      /* physical coordinate of interpolation node */
      cubx[cnt] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
      cuby[cnt] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
      ++cnt;
    }
  }
}

void mesh_t::ElementIntNodesTri2D(dlong e, memory<dfloat> &intx, memory<dfloat> &inty){

  dlong cnt = 0, icnt = 0, Icnt = 0;
  dlong ElementNfaces = ElementNfaces(etype);
  for(int f=0;f<ElementNfaces;++f){
    dlong intNfpFace = ElementIntNfp(etype, f);
    dlong NfpFace = ElementNfp(etype, f);
    for(int n=0;n<intNfpFace;++n){
      dfloat ix = 0, iy = 0;
      for(int m=0;m<NfpFace;++m){
	dlong vid = vmapM[cnt+m];
	dfloat xm = x[vid];
	dfloat ym = y[vid];
	dfloat Inm = intInterpHybrid[etype][m+n*NfpFace+Icnt]; // could make this matrix of matrices
	ix += Inm*xm;
	iy += Inm*ym;
      }
      intx[icnt] = ix;
      inty[icnt] = iy;
      ++icnt;
    }

    Icnt += intNfpFace*NfpFace;
    cnt += NfpFace;
  }
}

  
void mesh_t::CubaturePhysicalNodesHybrid2D(){

  dlong cubNpTotal = 0;
  for(dlong e=0;e<Nelements;++e){ /* for each element */
    Mesh::ElementType etype = elementTypes[e];
    cubNpTotal += cubatureNp(etype);
  }
    
  cubx.malloc(cubNpTotal);
  cuby.malloc(cubNpTotal);

  dlong cnt = 0;
  for(dlong e=0;e<Nelements;++e){ /* for each element */
    
    Mesh::ElementType etype = elementTypes[e];
    
    switch(etype){
    case TRIANGLES:
      ElementCubaturePhysicalNodesTri2D(e, cubx+cnt, cuby+cnt); break;
    case QUADRILATERALS:
      ElementCubaturePhysicalNodesQuad2D(e, cubx+cnt, cuby+cnt); break;
    default:
      std::cout << "mesh::CubatureNodesHybrid2D UNKNOWN ELEMENT TYPE" << std::endl;
      exit(-1); // correct exit call ?
    }

    cnt += ElementCubatureNp(etype);
  }

  o_cubx = platform.malloc<dfloat>(cubNpTotal, cubx);
  o_cuby = platform.malloc<dfloat>(cubNpTotal, cuby);

  dlong intTotal = 0;
  for(dlong e=0;e<Nelements;++e){ /* for each element */
    Mesh::ElementType etype = elementTypes[e];
    intTotal += ElementIntNfp(etype); // one arg => all nodes for all faces
  }
  
  //Face cubature
  intx.malloc(intTotal);
  inty.malloc(intTotal);
  cnt = 0;
  for(dlong e=0;e<Nelements;++e){

    Mesh::ElementType etype = elementTypes[e];
    
    switch(etype){
    case TRIANGLES:
      ElementIntNodesTri2D(e, intx+cnt, inty+cnt); break;
    case QUADRILATERALS:
      ElementIntNodesQuad2D(e, intx+cnt, inty+cnt); break;
    default:
      std::cout << "mesh::CubatureNodesHybrid2D UNKNOWN ELEMENT TYPE" << std::endl;
      exit(-1); // correct exit call ?
    }

    cnt += ElementIntNfp(etype);
  }
  
  o_intx = platform.malloc<dfloat>(intTotal, intx);
  o_inty = platform.malloc<dfloat>(intTotal, inty);
}

} //namespace libp
