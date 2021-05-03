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

#include "elliptic.hpp"

void elliptic_t::BoundarySetup(){

  //check all the bounaries for a Dirichlet
  int localAllNeumann = (lambda==0) ? 1 : 0; //if lambda>0 we don't care about all Neumann problem

  //setup a custom element-to-boundaryflag mapping
  EToB = (int *) calloc(mesh.Nelements*mesh.Nfaces,sizeof(int));
  for (dlong e=0;e<mesh.Nelements;e++) {
    for (int f=0;f<mesh.Nfaces;f++) {
      int bc = mesh.EToB[e*mesh.Nfaces+f];
      if (bc>0) {
        int BC = BCType[bc];  //translate mesh's boundary flag
        EToB[e*mesh.Nfaces+f] = BC;    //record it
        if (BC!=2) localAllNeumann = 0;     //check if its a Dirchlet
      }
    }
  }
  o_EToB = platform.malloc(mesh.Nelements*mesh.Nfaces*sizeof(int), EToB);

  //collect the allNeumann flags from other ranks
  MPI_Allreduce(&localAllNeumann, &allNeumann, 1, MPI_INT, MPI_MIN, mesh.comm);


  //make a node-wise bc flag (prioritize Dirichlet boundaries over Neumann)
  mapB = (int *) calloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np,sizeof(int));
  const int largeNumber = 1<<20;
  for (dlong e=0;e<mesh.Nelements;e++) {
    for (int n=0;n<mesh.Np;n++) mapB[n+e*mesh.Np] = largeNumber;
    for (int f=0;f<mesh.Nfaces;f++) {
      int bc = EToB[f+e*mesh.Nfaces];
      if (bc>0) {
        for (int n=0;n<mesh.Nfp;n++) {
          int fid = mesh.faceNodes[n+f*mesh.Nfp];
          mapB[fid+e*mesh.Np] = mymin(bc,mapB[fid+e*mesh.Np]);
        }
      }
    }
  }

  hlong localChange = 0, gatherChange = 1;

  // keep comparing numbers on positive and negative traces until convergence
  while(gatherChange>0){

    // reset change counter
    localChange = 0;

    // send halo data and recv into extension of buffer
    mesh.halo->Exchange(mapB, mesh.Np, ogs_int);

    // compare trace vertices
    for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.Nfaces*mesh.Nfp;++n){
        dlong id  = e*mesh.Nfaces*mesh.Nfp + n;
        dlong idM = mesh.vmapM[id];
        dlong idP = mesh.vmapP[id];

        int mapBM = mapB[idM];
        int mapBP = mapB[idP];

        if(mapBP<mapBM){
          localChange=1;
          mapB[idM] = mapBP;
        }
      }
    }

    // sum up changes
    MPI_Allreduce(&localChange, &gatherChange, 1, MPI_HLONG, MPI_MAX, mesh.comm);
  }

  //use the bc flags to find masked ids
  Nmasked = 0;
  for (dlong n=0;n<mesh.Nelements*mesh.Np;n++) {
    if (mapB[n] == largeNumber) {//no boundary
      mapB[n] = 0;
    } else if (mapB[n] == 1) {   //Dirichlet boundary
      Nmasked++;
    }
  }
  o_mapB = platform.malloc(mesh.Nelements*mesh.Np*sizeof(int), mapB);


  maskIds = (dlong *) calloc(Nmasked, sizeof(dlong));
  Nmasked =0; //reset
  for (dlong n=0;n<mesh.Nelements*mesh.Np;n++)
    if (mapB[n] == 1) maskIds[Nmasked++] = n;

  if (Nmasked) o_maskIds = platform.malloc(Nmasked*sizeof(dlong), maskIds);

  //make a masked version of the global id numbering
  maskedGlobalIds = (hlong *) calloc(mesh.Nelements*mesh.Np,sizeof(hlong));
  memcpy(maskedGlobalIds, mesh.globalIds, mesh.Nelements*mesh.Np*sizeof(hlong));
  for (dlong n=0;n<Nmasked;n++)
    maskedGlobalIds[maskIds[n]] = 0;

  //use the masked ids to make another gs handle (signed so the gather is defined)
  int verbose = 0;
  ogs_t::Unique(maskedGlobalIds, mesh.Nelements*mesh.Np, mesh.comm);     //flag a unique node in every gather node
  ogsMasked = ogs_t::Setup(mesh.Nelements*mesh.Np, maskedGlobalIds,
                           mesh.comm, verbose, platform);

  /* use the masked gs handle to define a global ordering */
  dlong Ntotal  = mesh.Np*mesh.Nelements; // number of degrees of freedom on this rank (before gathering)
  hlong Ngather = ogsMasked->Ngather;     // number of degrees of freedom on this rank (after gathering)

  //setup normalization constant
  allNeumannPenalty = 1.;
  if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    allNeumannScale = 1./sqrt((dfloat)mesh.Np*mesh.NelementsGlobal);
  } else {
    //note that we can use the mesh ogs, since there are no masked nodes
    allNeumannScale = 1./sqrt((dfloat)ogsMasked->NgatherGlobal);
  }

  // build inverse degree vectors
  // used for the weight in linear solvers (used in C0)
  weight  = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  weightG = (dfloat*) calloc(ogsMasked->Ngather, sizeof(dfloat));
  for(dlong n=0;n<Ntotal;++n) weight[n] = 1.0;

  ogsMasked->Gather(weightG, weight, ogs_dfloat, ogs_add, ogs_trans);
  for(dlong n=0;n<ogsMasked->Ngather;++n)
    if (weightG[n]) weightG[n] = 1./weightG[n];

  ogsMasked->Scatter(weight, weightG, ogs_dfloat, ogs_add, ogs_notrans);

  o_weight  = platform.malloc(Ntotal*sizeof(dfloat), weight);
  o_weightG = platform.malloc(ogsMasked->Ngather*sizeof(dfloat), weightG);

  // create a global numbering system
  hlong *globalIds = (hlong *) calloc(Ngather,sizeof(hlong));

  // every gathered degree of freedom has its own global id
  hlong *globalStarts = (hlong*) calloc(mesh.size+1,sizeof(hlong));
  MPI_Allgather(&Ngather, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh.comm);
  for(int r=0;r<mesh.size;++r)
    globalStarts[r+1] = globalStarts[r] + globalStarts[r+1];

  //use the offsets to set a consecutive global numbering
  for (dlong n =0;n<ogsMasked->Ngather;n++) {
    globalIds[n] = n + globalStarts[mesh.rank];
  }

  //scatter this numbering to the original nodes
  maskedGlobalNumbering = (hlong *) calloc(Ntotal,sizeof(hlong));
  for (dlong n=0;n<Ntotal;n++) maskedGlobalNumbering[n] = -1;
  ogsMasked->Scatter(maskedGlobalNumbering, globalIds, ogs_hlong, ogs_add, ogs_notrans);
  free(globalIds);

  /* Build halo exchange for gathered ordering */
  ogsMasked->GatheredHaloExchangeSetup();
}
