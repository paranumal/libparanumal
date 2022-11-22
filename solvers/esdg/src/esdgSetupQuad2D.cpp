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
  
  
}

