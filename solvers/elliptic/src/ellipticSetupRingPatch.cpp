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

#include "elliptic.hpp"
#include "ellipticPrecon.hpp"

elliptic_t elliptic_t::SetupRingPatch(mesh_t& meshPatch){

  //just reuse the current solver if there are no neighbors
  if (mesh.size == 1) return *this;

  //shallow copy
  elliptic_t elliptic = *this;

  elliptic.mesh = meshPatch;
  elliptic.comm = meshPatch.comm;

  //buffer for gradient
  if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    dlong Ntotal = meshPatch.Np*meshPatch.Nelements;
    elliptic.grad.malloc(Ntotal*4, 0.0);
    elliptic.o_grad = platform.malloc<dfloat>(elliptic.grad);
  } else {
    //buffer for local Ax
    dlong Ntotal = meshPatch.Np*meshPatch.Nelements;
    elliptic.o_AqL = platform.malloc<dfloat>(Ntotal);
  }

  /*setup trace halo exchange */
  elliptic.traceHalo = meshPatch.HaloTraceSetup(Nfields);

  //setup boundary flags and make mask and masked ogs
  elliptic.BoundarySetup();

  if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    elliptic.Ndofs = elliptic.ogsMasked.Ngather*Nfields;
  } else {
    elliptic.Ndofs = meshPatch.Nelements*meshPatch.Np*Nfields;
  }

  elliptic.precon = precon_t();

  return elliptic;
}
