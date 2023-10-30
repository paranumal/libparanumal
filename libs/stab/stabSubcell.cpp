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

#include "stab.hpp"

namespace libp {

void stab_t::stabSetupSubcell(){
  /*

  settings.getSetting("SUBCELL NUMBER", N);

  if(N<mesh.N){
     LIBP_FORCE_ABORT("Subcell number can not be less than N");
   }

   switch (mesh.elementType) {
      case Mesh::TRIANGLES:
        stabSetupSubcellTri2D();
        break;
      case Mesh::QUADRILATERALS:
        // stabSetupSubcellQuad2D();
         // LIBP_FORCE_ABORT("Limiter is not implemented yet");
        break;
      case Mesh::TETRAHEDRA:
        // stabSetupSubcellTet3D();
        break;
      case Mesh::HEXAHEDRA:
        // stabSetupSubcellHex3D();
        // LIBP_FORCE_ABORT("FV-Subcell is not implemented yet");
        break;
    } 


   sq.malloc((mesh.Nelements + mesh.totalHaloPairs)*Nsubcells*sNfields, 0.0); 
   o_sq = platform.malloc<dfloat>(sq);

   sqf.malloc((mesh.Nelements + mesh.totalHaloPairs)*Nsubcells*sNfields*Nfaces, 0.0); 
   o_sqf= platform.malloc<dfloat>(sqf);

   srhsq.malloc((mesh.Nelements + mesh.totalHaloPairs)*Nsubcells*sNfields, 0.0); 
   o_srhsq = platform.malloc<dfloat>(srhsq);

    */
}


} //namespace libp
