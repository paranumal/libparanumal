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

#include "mesh.hpp"
#include "mesh3D.hpp"

void mesh3D::OccaSetup(){

  this->mesh_t::OccaSetup();

  o_x = device.malloc(Nelements*Np*sizeof(dfloat), x);
  o_y = device.malloc(Nelements*Np*sizeof(dfloat), y);
  o_z = device.malloc(Nelements*Np*sizeof(dfloat), z);

  props["defines/" "p_NXID"]= NXID;
  props["defines/" "p_NYID"]= NYID;
  props["defines/" "p_NZID"]= NZID;
  props["defines/" "p_SJID"]= SJID;
  props["defines/" "p_IJID"]= IJID;
  props["defines/" "p_IHID"]= IHID;
  props["defines/" "p_WSJID"]= WSJID;
  props["defines/" "p_WIJID"]= WIJID;
  props["defines/" "p_STXID"]= STXID;
  props["defines/" "p_STYID"]= STYID;
  props["defines/" "p_STZID"]= STZID;
  props["defines/" "p_SBXID"]= SBXID;
  props["defines/" "p_SBYID"]= SBYID;
  props["defines/" "p_SBZID"]= SBZID;

  props["defines/" "p_G00ID"]= G00ID;
  props["defines/" "p_G01ID"]= G01ID;
  props["defines/" "p_G02ID"]= G02ID;
  props["defines/" "p_G11ID"]= G11ID;
  props["defines/" "p_G12ID"]= G12ID;
  props["defines/" "p_G22ID"]= G22ID;
  props["defines/" "p_GWJID"]= GWJID;


  props["defines/" "p_RXID"]= RXID;
  props["defines/" "p_SXID"]= SXID;
  props["defines/" "p_TXID"]= TXID;

  props["defines/" "p_RYID"]= RYID;
  props["defines/" "p_SYID"]= SYID;
  props["defines/" "p_TYID"]= TYID;

  props["defines/" "p_RZID"]= RZID;
  props["defines/" "p_SZID"]= SZID;
  props["defines/" "p_TZID"]= TZID;

  props["defines/" "p_JID"]= JID;
  props["defines/" "p_JWID"]= JWID;
  props["defines/" "p_IJWID"]= IJWID;
}
