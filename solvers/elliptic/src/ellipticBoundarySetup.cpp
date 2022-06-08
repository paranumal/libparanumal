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
#include <limits>

void elliptic_t::BoundarySetup(){

  //check all the bounaries for a Dirichlet
  allNeumann = (lambda==0) ? 1 : 0; //if lambda>0 we don't care about all Neumann problem
  allNeumannPenalty = 1.;

  //translate the mesh's element-to-boundaryflag mapping
  EToB.malloc(mesh.Nelements*mesh.Nfaces, 0);
  for (dlong e=0;e<mesh.Nelements;e++) {
    for (int f=0;f<mesh.Nfaces;f++) {
      int bc = mesh.EToB[e*mesh.Nfaces+f];
      if (bc>0) {
        int BC = BCType[bc];         //translate mesh's boundary flag
        EToB[e*mesh.Nfaces+f] = BC;  //record it
        if (BC!=2) allNeumann = 0;   //check if its a Dirchlet
      }
    }
  }
  o_EToB = platform.malloc<int>(EToB);

  //collect the allNeumann flags from other ranks
  mesh.comm.Allreduce(allNeumann, Comm::Min);

  //translate the mesh's node-wise bc flag
  Nmasked = 0;
  mapB.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np, 0);
  for (int n=0;n<mesh.Nelements*mesh.Np;n++) {
    int bc = mesh.mapB[n];
    if (bc>0) {
      int BC = BCType[bc];     //translate mesh's boundary flag
      mapB[n] = BC;  //record it

      if (mapB[n] == 1) Nmasked++;   //Dirichlet boundary
    }
  }
  o_mapB = platform.malloc<int>(mapB);

  maskIds.malloc(Nmasked);
  Nmasked =0; //reset
  for (dlong n=0;n<mesh.Nelements*mesh.Np;n++) {
    if (mapB[n] == 1) maskIds[Nmasked++] = n;
  }
  o_maskIds = platform.malloc<int>(maskIds);

  //make a masked version of the global id numbering
  maskedGlobalIds.malloc(mesh.Nelements*mesh.Np);
  maskedGlobalIds.copyFrom(mesh.globalIds);
  for (dlong n=0;n<Nmasked;n++) {
    maskedGlobalIds[maskIds[n]] = 0;
  }

  //use the masked ids to make another gs handle (signed so the gather is defined)
  bool verbose = settings.compareSetting("VERBOSE", "TRUE") ? true : false;
  bool unique = true; //flag a unique node in every gather node
  ogsMasked.Setup(mesh.Nelements*mesh.Np, maskedGlobalIds,
                  mesh.comm, ogs::Signed, ogs::Auto,
                  unique, verbose, platform);

  //setup normalization constant
  if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    allNeumannScale = 1./sqrt((dfloat)mesh.Np*mesh.NelementsGlobal);
  } else {
    //note that we can use the mesh ogs, since there are no masked nodes
    allNeumannScale = 1./sqrt((dfloat)ogsMasked.NgatherGlobal);
  }

  /* use the masked gs handle to define a global ordering */
  dlong Ntotal  = mesh.Np*mesh.Nelements; // number of degrees of freedom on this rank (before gathering)
  hlong Ngather = ogsMasked.Ngather;     // number of degrees of freedom on this rank (after gathering)

  // build inverse degree vectors
  // used for the weight in linear solvers (used in C0)
  weight.malloc(Ntotal, 1.0);

  weightG.malloc(Ngather);
  ogsMasked.Gather(weightG, weight, 1, ogs::Add, ogs::Trans);

  for(dlong n=0;n<ogsMasked.Ngather;++n) {
    if (weightG[n]>0.0) weightG[n] = 1./weightG[n];
  }

  ogsMasked.Scatter(weight, weightG, 1, ogs::NoTrans);

  o_weight  = platform.malloc<dfloat>(weight);
  o_weightG = platform.malloc<dfloat>(weightG);

  // create a global numbering system
  memory<hlong> globalIds(Ngather);

  // every gathered degree of freedom has its own global id
  hlong globalOffset=static_cast<hlong>(Ngather);
  comm.Scan(Ngather, globalOffset);
  globalOffset = globalOffset-Ngather;

  //use the offsets to set a consecutive global numbering
  for (dlong n =0;n<ogsMasked.Ngather;n++) {
    globalIds[n] = n + globalOffset;
  }

  //scatter this numbering to the original nodes
  maskedGlobalNumbering.malloc(Ntotal, -1);
  ogsMasked.Scatter(maskedGlobalNumbering, globalIds, 1, ogs::NoTrans);

  /* Build halo exchange for gathered ordering */
  gHalo.SetupFromGather(ogsMasked);

  GlobalToLocal.malloc(mesh.Nelements*mesh.Np);
  ogsMasked.SetupGlobalToLocalMapping(GlobalToLocal);

  o_GlobalToLocal = platform.malloc<dlong>(GlobalToLocal);
}
