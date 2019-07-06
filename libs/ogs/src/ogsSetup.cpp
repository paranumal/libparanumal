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

#include "ogs.hpp"
#include "ogsKernels.hpp"

typedef struct{

  dlong localId;    // local node id
  hlong baseId;     // original global index

  dlong newId;         // new global id
  int owned;

}parallelNode_t;

// compare on baseId then by localId
static int compareBaseId(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  if(abs(fa->baseId) < abs(fb->baseId)) return -1; //group by abs(baseId)
  if(abs(fa->baseId) > abs(fb->baseId)) return +1;

  if(fa->localId < fb->localId) return -1; //sort by local id
  if(fa->localId > fb->localId) return +1;

  return 0;
}

// compare on haloOwned then localId
static int compareLocalId(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  if(fa->localId < fb->localId) return -1;
  if(fa->localId > fb->localId) return +1;

  return 0;
}

ogs_t *ogs_t::Setup(dlong N, hlong *ids, MPI_Comm &comm,
                    int verbose, occa::device &device){

  ogs_t *ogs = new ogs_t(comm, device);

  ogs->NhaloGather = 0;
  ogs->Ngather = 0;
  ogs->Nlocal = 0;
  ogs->NlocalGather = 0;
  ogs->Nhalo = 0;
  ogs->NhaloGather = 0;
  ogs->NownedHalo = 0;

  ogs->localGatherOffsets = NULL;
  ogs->localGatherIds = NULL;

  ogs->haloGatherOffsets = NULL;
  ogs->haloGatherIds = NULL;

  ogs->gshSym = NULL;
  ogs->gshNonSym = NULL;

  //Keep track of how many gs handles we've created, and
  // build kernels if this is the first
  if (!ogs::Nrefs) ogs::initKernels(comm, device);
  ogs::Nrefs++;

  ogs->N = N;

  int rank, size;
  MPI_Comm_rank(ogs->comm, &rank);
  MPI_Comm_size(ogs->comm, &size);

  //use the host gs to find what nodes are local to this rank
  int *minRank = (int *) calloc(N,sizeof(int));
  int *maxRank = (int *) calloc(N,sizeof(int));
  hlong *flagIds   = (hlong *) calloc(N,sizeof(hlong));
  for (dlong i=0;i<N;i++) {
    minRank[i] = rank;
    maxRank[i] = rank;
    flagIds[i] = abs(ids[i]); //note that we ignore negative ids in setup and flag our own
  }

  //make a host gs handle (calls gslib)
  void *gsHandle = ogs::gsSetup(comm, N, flagIds, 0, 0);

  ogs::gsGatherScatter(minRank, 1, 1, 0, ogs_int, ogs_min, gsHandle); //minRank[n] contains the smallest rank taking part in the gather of node n
  ogs::gsGatherScatter(maxRank, 1, 1, 0, ogs_int, ogs_max, gsHandle); //maxRank[n] contains the largest rank taking part in the gather of node n
  ogs::gsUnique(flagIds, N, comm); //one unique node in each group is 'flagged' kept positive while others are turned negative.

  //discard the large gs handle
  ogs::gsFree(gsHandle);

  //count local and halo nodes
  ogs->Nlocal=0; ogs->Nhalo=0; ogs->NownedHalo=0;
  for (dlong i=0;i<N;i++) {
    if (flagIds[i]==0) continue;

    if ((minRank[i]!=rank)||(maxRank[i]!=rank)) {
      ogs->Nhalo++;
      if (flagIds[i]>0) ogs->NownedHalo++;
    } else {
      ogs->Nlocal++;
    }
  }

  //set up the local gatherScatter
  parallelNode_t *localNodes = (parallelNode_t*) calloc(ogs->Nlocal,sizeof(parallelNode_t));

  dlong cnt=0;
  for (dlong i=0;i<N;i++) {
    if (flagIds[i]==0) continue;

    if ((minRank[i]==rank)&&(maxRank[i]==rank)) {
      localNodes[cnt].localId = i;
      localNodes[cnt].baseId  = flagIds[i];
      localNodes[cnt].owned   = 0;
      cnt++;
    }
  }

  // sort based on base ids then local id
  qsort(localNodes, ogs->Nlocal, sizeof(parallelNode_t), compareBaseId);

  if (ogs->Nlocal) {
    ogs->NlocalGather = 0;
    localNodes[0].newId = 0;
    localNodes[0].owned = 1;
    for (dlong i=1;i<ogs->Nlocal;i++) {
      int s = 0;
      if (abs(localNodes[i].baseId)!=abs(localNodes[i-1].baseId)) {
        ogs->NlocalGather++;
        s = 1;
      }
      localNodes[i].newId = ogs->NlocalGather;
      localNodes[i].owned = s;
    }
    ogs->NlocalGather++;
  }

  // sort based on local ids
  qsort(localNodes, ogs->Nlocal, sizeof(parallelNode_t), compareLocalId);

  //tally up how many nodes are being gathered to each gatherNode and
  //  map to a local ordering
  dlong *localGatherCounts = (dlong*) calloc(ogs->NlocalGather,sizeof(dlong));
  dlong *localGatherMap    = (dlong*) calloc(ogs->NlocalGather,sizeof(dlong));
  cnt = 0;
  for (dlong i=0;i<ogs->Nlocal;i++) {
    dlong newId = localNodes[i].newId; //get the ordered id

    if (localNodes[i].owned)
      localGatherMap[newId] = cnt++; //record a new index if this is a new gatherNode

    localNodes[i].newId = localGatherMap[newId]; //reorder
    localGatherCounts[localGatherMap[newId]]++;  //tally
  }
  free(localGatherMap);

  ogs->localGatherOffsets = (dlong*) calloc(ogs->NlocalGather+1,sizeof(dlong));
  for (dlong i=0;i<ogs->NlocalGather;i++) {
    ogs->localGatherOffsets[i+1] = ogs->localGatherOffsets[i] + localGatherCounts[i];
    localGatherCounts[i] = 0;
  }

  ogs->localGatherIds = (dlong*) calloc(ogs->Nlocal+1,sizeof(dlong)); //extra entry so the occa buffer will actually exist
  for (dlong i=0;i<ogs->Nlocal;i++) {
    dlong gatherId = localNodes[i].newId;
    dlong offset = ogs->localGatherOffsets[gatherId];
    int index  = localGatherCounts[gatherId];

    ogs->localGatherIds[offset+index] = localNodes[i].localId;
    localGatherCounts[gatherId]++;
  }
  free(localGatherCounts);

  ogs->o_localGatherOffsets = device.malloc((ogs->NlocalGather+1)*sizeof(dlong), ogs->localGatherOffsets);
  ogs->o_localGatherIds     = device.malloc((ogs->Nlocal+1)*sizeof(dlong), ogs->localGatherIds);

  free(localNodes);

  //set up the halo gatherScatter
  parallelNode_t *haloNodes = (parallelNode_t*) calloc(ogs->Nhalo+1,sizeof(parallelNode_t));

  cnt=0;
  for (dlong i=0;i<N;i++) {
    if (flagIds[i]==0) continue;

    if ((minRank[i]!=rank)||(maxRank[i]!=rank)) {
      haloNodes[cnt].localId = i;
      haloNodes[cnt].baseId  = flagIds[i];
      haloNodes[cnt].owned   = 0;
      cnt++;
    }
  }

  // sort based on base ids then local id
  qsort(haloNodes, ogs->Nhalo, sizeof(parallelNode_t), compareBaseId);

  cnt = 0;
  ogs->NhaloGather=0;

  if(ogs->Nhalo){
    //move the flagged node to the lowest local index if present
    haloNodes[0].newId = 0;
    haloNodes[0].owned = 1;

    for (dlong i=1;i<ogs->Nhalo;i++) {
      int s = 0;
      if (abs(haloNodes[i].baseId)!=abs(haloNodes[i-1].baseId)) { //new gather node
        s = 1;
        cnt = i;
        ogs->NhaloGather++;
      }

      haloNodes[i].owned = s;
      haloNodes[i].newId = ogs->NhaloGather;
      if (haloNodes[i].baseId>0) {
        haloNodes[i].baseId   = -abs(haloNodes[i].baseId);
        haloNodes[cnt].baseId =  abs(haloNodes[cnt].baseId);
      }
    }
    ogs->NhaloGather++;
  }

  // sort based on local ids
  qsort(haloNodes, ogs->Nhalo, sizeof(parallelNode_t), compareLocalId);

  //tally up how many nodes are being gathered to each gatherNode and
  //  map to a local ordering
  dlong *haloGatherCounts = (dlong*) calloc(ogs->NhaloGather+1,sizeof(dlong));
  dlong *haloGatherMap    = (dlong*) calloc(ogs->NhaloGather+1,sizeof(dlong));
  hlong *symIds    = (hlong *) calloc(ogs->NhaloGather+1,sizeof(hlong));
  hlong *nonSymIds = (hlong *) calloc(ogs->NhaloGather+1,sizeof(hlong));

  cnt = 0;
  dlong cnt2 = ogs->NownedHalo;
  for (dlong i=0;i<ogs->Nhalo;i++) {
    dlong newId = haloNodes[i].newId; //get the ordered id

    if (haloNodes[i].owned) {
      dlong c;
      if (haloNodes[i].baseId>0)
        c = cnt++;  //all the owned nodes are packed at the beginning
      else
        c = cnt2++; //all the un-owned nodes are at the end

      symIds[c]    = abs(haloNodes[i].baseId); //record the base id
      nonSymIds[c] = haloNodes[i].baseId;      //record the base id
      haloGatherMap[newId] = c; //record a new index if this is a new gatherNode
    }

    haloNodes[i].newId = haloGatherMap[newId];  //reorder
    haloGatherCounts[haloGatherMap[newId]]++;  //tally
  }
  free(haloGatherMap);

  ogs->haloGatherOffsets = (dlong*) calloc(ogs->NhaloGather+1,sizeof(dlong));
  for (dlong i=0;i<ogs->NhaloGather;i++) {
    ogs->haloGatherOffsets[i+1] = ogs->haloGatherOffsets[i] + haloGatherCounts[i];
    haloGatherCounts[i] = 0;
  }

  ogs->haloGatherIds = (dlong*) calloc(ogs->Nhalo+1,sizeof(dlong));
  for (dlong i=0;i<ogs->Nhalo;i++) {
    dlong gatherId = haloNodes[i].newId;
    dlong offset = ogs->haloGatherOffsets[gatherId];
    int index  = haloGatherCounts[gatherId];

    ogs->haloGatherIds[offset+index] = haloNodes[i].localId;
    haloGatherCounts[gatherId]++;
  }
  free(haloGatherCounts);

  ogs->o_haloGatherOffsets = device.malloc((ogs->NhaloGather+1)*sizeof(dlong), ogs->haloGatherOffsets);
  ogs->o_haloGatherIds     = device.malloc((ogs->Nhalo+1)*sizeof(dlong), ogs->haloGatherIds);

  //make a host gs handle
  ogs->gshSym    = ogs::gsSetup(comm, ogs->NhaloGather, symIds,    0,0);
  ogs->gshNonSym = ogs::gsSetup(comm, ogs->NhaloGather, nonSymIds, 0,0);

  free(symIds); free(nonSymIds);
  free(haloNodes);

  free(minRank); free(maxRank); free(flagIds);

  //total number of owned gathered nodes
  ogs->Ngather = ogs->NlocalGather+ogs->NownedHalo;

  ogs->hostBuf = NULL;
  ogs->haloBuf = NULL;
  ogs->hostBufSize = 0;

  return ogs;
}


void ogs_t::Free() {

  if (Nlocal) {
    free(localGatherOffsets);
    free(localGatherIds);
    o_localGatherOffsets.free();
    o_localGatherIds.free();
    Nlocal = 0;
  }

  if (Nhalo) {
    free(haloGatherOffsets);
    free(haloGatherIds);
    o_haloGatherOffsets.free();
    o_haloGatherIds.free();
    Nhalo = 0;
  }

  ogs::gsFree(gshSym);
  ogs::gsFree(gshNonSym);

  ogs::Nrefs--;
  if (!ogs::Nrefs) ogs::freeKernels();
}

void ogs_t::reallocHostBuffer(size_t Nbytes) {
  if (NhaloGather) {
    if (hostBufSize < NhaloGather*Nbytes) {
      if (hostBufSize) free(hostBuf);
      hostBuf = (void *) malloc(NhaloGather*Nbytes);
      hostBufSize = NhaloGather*Nbytes;
    }
  }
}

void ogs_t::reallocOccaBuffer(size_t Nbytes) {
  if (NhaloGather) {
    if (o_haloBuf.size() < NhaloGather*Nbytes) {
      if (o_haloBuf.size()) o_haloBuf.free();
      haloBuf = occaHostMallocPinned(device, NhaloGather*Nbytes, NULL,
                                     o_haloBuf, h_haloBuf);
    }
  }
}