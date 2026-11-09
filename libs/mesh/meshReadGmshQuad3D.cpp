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
#include "mesh/meshGmsh.hpp"

namespace libp {

/*
   purpose: read gmsh quadrilateral mesh (.msh format 2.2 or 4.1)
*/
void mesh_t::ReadGmshQuad3D(const std::string fileName){

  gmsh::reader_t gmshFile(fileName);

  /* node coordinates */
  Nnodes = gmshFile.Nnodes;
  memory<dfloat> VX = gmshFile.VX;
  memory<dfloat> VY = gmshFile.VY;
  memory<dfloat> VZ = gmshFile.VZ;

  /* find # of quadrilaterals and boundary faces */
  std::unordered_map<int, hlong> counts = gmshFile.CountElements();
  hlong Nquadrilaterals = counts[gmsh::QUADRILATERAL];
  hlong gNboundaryFaces = counts[gmsh::LINE];

  hlong chunk = (hlong) Nquadrilaterals/size;
  int remainder = (int) (Nquadrilaterals - chunk*size);

  hlong NquadrilateralsLocal = chunk + (rank<remainder);

  /* where do these elements start ? */
  hlong start = rank*chunk + std::min(rank, remainder);
  hlong end   = start + NquadrilateralsLocal-1;

  /* allocate space for Element node index data */
  EToV.malloc(NquadrilateralsLocal*Nverts);
  elementInfo.malloc(NquadrilateralsLocal);

  /* allocate space for boundary faces (tag + face vertices) */
  boundaryInfo.malloc(gNboundaryFaces*(NfaceVertices+1));

  /* scan through file looking for quadrilaterals and boundary faces */
  hlong cnt=0, bcnt=0;
  Nquadrilaterals = 0;

  gmshFile.ForEachElement([&](const int gmshType,
                              const hlong physicalTag,
                              const hlong* v) {
    if(gmshType==gmsh::LINE){ // boundary face
      boundaryInfo[bcnt*(NfaceVertices+1)] = physicalTag;
      for(int n=0;n<NfaceVertices;++n)
        boundaryInfo[bcnt*(NfaceVertices+1)+n+1] = v[n];
      ++bcnt;
    }
    if(gmshType==gmsh::QUADRILATERAL){ // quadrilateral
      if(start<=Nquadrilaterals && Nquadrilaterals<=end){
        elementInfo[cnt] = physicalTag;
        for(int n=0;n<Nverts;++n)
          EToV[cnt*Nverts+n] = v[n];

        ++cnt;
      }
      ++Nquadrilaterals;
    }
  });

  /* record number of boundary faces found */
  NboundaryFaces = bcnt;

  /* record number of found quadrilaterals */
  Nelements = (dlong) NquadrilateralsLocal;

  /* collect vertices for each element */
  EX.malloc(Nverts*Nelements);
  EY.malloc(Nverts*Nelements);
  EZ.malloc(Nverts*Nelements);
  for(dlong e=0;e<Nelements;++e){
    for(int n=0;n<Nverts;++n){
      hlong vid = EToV[e*Nverts+n];
      EX[e*Nverts+n] = VX[vid];
      EY[e*Nverts+n] = VY[vid];
      EZ[e*Nverts+n] = VZ[vid];
    }
  }
}

} //namespace libp
