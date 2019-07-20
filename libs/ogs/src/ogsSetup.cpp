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
  int sign;

}parallelNode_t;

// compare on baseId then by localId
static int compareBaseId(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  if(abs(fa->baseId) < abs(fb->baseId)) return -1; //group by abs(baseId)
  if(abs(fa->baseId) > abs(fb->baseId)) return +1;

  if(fa->baseId > fb->baseId) return -1; //positive ids first
  if(fa->baseId < fb->baseId) return +1;

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

  ogs->Ngather = 0;

  ogs->Nlocal = 0;
  ogs->NlocalGather = 0;

  ogs->Nhalo = 0;
  ogs->NhaloGather = 0;
  ogs->NhaloScatter = 0;

  ogs->localGatherOffsets = NULL;
  ogs->localGatherIds = NULL;
  ogs->localScatterOffsets = NULL;
  ogs->localScatterIds = NULL;

  ogs->haloGatherOffsets = NULL;
  ogs->haloGatherIds = NULL;
  ogs->haloScatterOffsets = NULL;
  ogs->haloScatterIds = NULL;

  ogs->gsh = NULL;

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
    flagIds[i] = abs(ids[i]); //ignore negative ids for this
  }

  //make a host gs handle (calls gslib)
  void *gsHandle = ogs::gsSetup(comm, N, flagIds, 0, 0);
  ogs::gsGatherScatter(minRank, 1, 1, 0, ogs_int, ogs_min, ogs_notrans, gsHandle);
  ogs::gsGatherScatter(maxRank, 1, 1, 0, ogs_int, ogs_max, ogs_notrans, gsHandle);
  ogs::gsFree(gsHandle); //discard the large gs handle
  free(flagIds);

  //minRank[n] contains the smallest rank taking part in the gatherScatter of node n
  //maxRank[n] contains the largest rank taking part in the gatherScatter of node n

  //count local and halo nodes
  ogs->Nlocal=0; ogs->Nhalo=0;
  for (dlong i=0;i<N;i++) {
    if (ids[i]==0) continue;

    if ((minRank[i]!=rank)||(maxRank[i]!=rank)) {
      ogs->Nhalo++;
    } else {
      ogs->Nlocal++;
    }
  }

  //set up the local gatherScatter
  parallelNode_t *localNodes = (parallelNode_t*) calloc(ogs->Nlocal,sizeof(parallelNode_t));

  dlong cnt=0;
  for (dlong i=0;i<N;i++) {
    if (ids[i]==0) continue;

    if ((minRank[i]==rank)&&(maxRank[i]==rank)) {
      localNodes[cnt].localId = i;
      localNodes[cnt].baseId  = ids[i];
      cnt++;
    }
  }

  // sort based on base ids (putting positive ids first) then local id
  qsort(localNodes, ogs->Nlocal, sizeof(parallelNode_t), compareBaseId);

  //flag each set of ids by whether there is at least one positive id
  // and count how many local gather/scatter nodes we have
  ogs->NlocalGather = 0;
  ogs->NlocalScatter = 0;
  if (ogs->Nlocal) {
    localNodes[0].newId = 0;
    int sign = (localNodes[0].baseId > 0) ? 1 : -1;
    localNodes[0].sign = sign;
    if (sign > 0) ogs->NlocalGather++;

    for (dlong i=1;i<ogs->Nlocal;i++) {
      if (abs(localNodes[i].baseId)!=abs(localNodes[i-1].baseId)) {
        sign = (localNodes[i].baseId > 0) ? 1 : -1;
        ogs->NlocalScatter++;
        if (sign > 0) ogs->NlocalGather++;
      }

      localNodes[i].newId = ogs->NlocalScatter;
      localNodes[i].sign = sign;
    }
    ogs->NlocalScatter++;
  }

  // sort back to local ids
  qsort(localNodes, ogs->Nlocal, sizeof(parallelNode_t), compareLocalId);

  //tally up how many nodes are being gathered to each gatherNode and
  //  map to a local ordering
  dlong *localGatherCounts  = (dlong*) calloc(ogs->NlocalScatter,sizeof(dlong));
  dlong *localScatterCounts = (dlong*) calloc(ogs->NlocalScatter,sizeof(dlong));
  dlong *localMap = (dlong*) calloc(ogs->NlocalScatter,sizeof(dlong));

  for (dlong i=0;i<ogs->NlocalScatter;i++) localMap[i] = -1; //initialize map

  cnt = 0;
  dlong cnt2 = ogs->NlocalGather;
  for (dlong i=0;i<ogs->Nlocal;i++) {
    dlong newId = localNodes[i].newId; //get the ordered id

    //record a new index if this is a new gatherNode (pure negative nodes appended at the end)
    if (localMap[newId]==-1) {
      if (localNodes[i].sign > 0)
        localMap[newId] = cnt++;
      else
        localMap[newId] = cnt2++;
    }

    dlong gid = localMap[newId];
    localNodes[i].newId = gid; //reorder
    localScatterCounts[gid]++;  //tally
    if (localNodes[i].baseId > 0)
      localGatherCounts[gid]++;  //tally
  }
  free(localMap);

  ogs->localGatherOffsets  = (dlong*) calloc(ogs->NlocalScatter+1,sizeof(dlong));
  ogs->localScatterOffsets = (dlong*) calloc(ogs->NlocalScatter+1,sizeof(dlong));
  for (dlong i=0;i<ogs->NlocalScatter;i++) {
    ogs->localGatherOffsets[i+1]  = ogs->localGatherOffsets[i]  + localGatherCounts[i];
    ogs->localScatterOffsets[i+1] = ogs->localScatterOffsets[i] + localScatterCounts[i];

    //reset counters
    localGatherCounts[i] = 0;
    localScatterCounts[i] = 0;
  }

  dlong totalLocalGather  = ogs->localGatherOffsets[ogs->NlocalScatter];
  dlong totalLocalScatter = ogs->localScatterOffsets[ogs->NlocalScatter];

  ogs->localGatherIds  = (dlong*) calloc(totalLocalGather+1,sizeof(dlong)); //extra entry so the occa buffer will actually exist
  ogs->localScatterIds = (dlong*) calloc(totalLocalScatter+1,sizeof(dlong)); //extra entry so the occa buffer will actually exist
  for (dlong i=0;i<ogs->Nlocal;i++) {
    dlong gid = localNodes[i].newId;

    dlong soffset = ogs->localScatterOffsets[gid];
    int sindex  = localScatterCounts[gid];
    ogs->localScatterIds[soffset+sindex] = localNodes[i].localId;
    localScatterCounts[gid]++;

    if (localNodes[i].baseId > 0) {
      dlong goffset = ogs->localGatherOffsets[gid];
      int gindex  = localGatherCounts[gid];
      ogs->localGatherIds[goffset+gindex] = localNodes[i].localId;
      localGatherCounts[gid]++;
    }
  }
  free(localGatherCounts);
  free(localScatterCounts);

  ogs->o_localGatherOffsets  = device.malloc((ogs->NlocalScatter+1)*sizeof(dlong), ogs->localGatherOffsets);
  ogs->o_localScatterOffsets = device.malloc((ogs->NlocalScatter+1)*sizeof(dlong), ogs->localScatterOffsets);

  ogs->o_localGatherIds  = device.malloc((totalLocalGather+1)*sizeof(dlong), ogs->localGatherIds);
  ogs->o_localScatterIds = device.malloc((totalLocalScatter+1)*sizeof(dlong), ogs->localScatterIds);

  free(localNodes);

  //make some compressed versions of the gather/scatter ids for the fused gs kernel
  ogs->NlocalFused=0;
  ogs->NlocalFusedSym=0;

  dlong totalGatherIdsFused = 0;
  dlong totalScatterIdsFused = 0;
  dlong totalIdsFusedSym = 0;

  for (dlong n=0;n<ogs->NlocalScatter;n++) {
    int gatherCnt  = ogs->localGatherOffsets[n+1] -ogs->localGatherOffsets[n];
    int scatterCnt = ogs->localScatterOffsets[n+1]-ogs->localScatterOffsets[n];

    //only include this node if both the gather and scatter interact with mulitple nodes
    // otherwise the op is identity and ignored
    if ((gatherCnt>1)||(scatterCnt>1)) {
      ogs->NlocalFused++;
      totalGatherIdsFused  += gatherCnt;
      totalScatterIdsFused += scatterCnt;
    }

    //for the sym op only the scatter ids are used
    if (scatterCnt>1) {
      ogs->NlocalFusedSym++;
      totalIdsFusedSym += scatterCnt;
    }
  }

  ogs->localFusedGatherOffsets  = (dlong*) calloc(ogs->NlocalFused+1,sizeof(dlong));
  ogs->localFusedScatterOffsets = (dlong*) calloc(ogs->NlocalFused+1,sizeof(dlong));
  ogs->localFusedOffsets = (dlong*) calloc(ogs->NlocalFusedSym+1,sizeof(dlong));
  ogs->localFusedGatherIds  = (dlong*) calloc(totalGatherIdsFused+1,sizeof(dlong));
  ogs->localFusedScatterIds = (dlong*) calloc(totalScatterIdsFused+1,sizeof(dlong));
  ogs->localFusedIds = (dlong*) calloc(totalIdsFusedSym+1,sizeof(dlong));

  //reset counters
  ogs->NlocalFused=0;
  ogs->NlocalFusedSym=0;
  totalGatherIdsFused = 0;
  totalScatterIdsFused = 0;
  totalIdsFusedSym = 0;
  for (dlong n=0;n<ogs->NlocalScatter;n++) {
    int gatherCnt  = ogs->localGatherOffsets[n+1] -ogs->localGatherOffsets[n];
    int scatterCnt = ogs->localScatterOffsets[n+1]-ogs->localScatterOffsets[n];

    //only include this node if both the gather and scatter interact with mulitple nodes
    // otherwise the op is identity and ignored
    if ((gatherCnt>1)||(scatterCnt>1)) {
      ogs->NlocalFused++;
      ogs->localFusedGatherOffsets[ogs->NlocalFused] = gatherCnt + ogs->localFusedGatherOffsets[ogs->NlocalFused-1];
      ogs->localFusedScatterOffsets[ogs->NlocalFused] = scatterCnt + ogs->localFusedScatterOffsets[ogs->NlocalFused-1];

      for (int i=ogs->localGatherOffsets[n];i<ogs->localGatherOffsets[n+1];i++)
        ogs->localFusedGatherIds[totalGatherIdsFused++] = ogs->localGatherIds[i];

      for (int i=ogs->localScatterOffsets[n];i<ogs->localScatterOffsets[n+1];i++)
        ogs->localFusedScatterIds[totalScatterIdsFused++] = ogs->localScatterIds[i];
    }

    //for the sym op only the scatter ids are used
    if (scatterCnt>1) {
      ogs->NlocalFusedSym++;
      ogs->localFusedOffsets[ogs->NlocalFusedSym] = scatterCnt + ogs->localFusedOffsets[ogs->NlocalFusedSym-1];

      for (int i=ogs->localScatterOffsets[n];i<ogs->localScatterOffsets[n+1];i++)
        ogs->localFusedIds[totalIdsFusedSym++] = ogs->localScatterIds[i];
    }
  }

  ogs->o_localFusedGatherOffsets  = device.malloc((ogs->NlocalFused+1)*sizeof(dlong), ogs->localFusedGatherOffsets);
  ogs->o_localFusedScatterOffsets = device.malloc((ogs->NlocalFused+1)*sizeof(dlong), ogs->localFusedScatterOffsets);
  ogs->o_localFusedOffsets = device.malloc((ogs->NlocalFusedSym+1)*sizeof(dlong), ogs->localFusedScatterOffsets);
  ogs->o_localFusedGatherIds  = device.malloc((totalGatherIdsFused+1)*sizeof(dlong),  ogs->localFusedGatherIds);
  ogs->o_localFusedScatterIds = device.malloc((totalScatterIdsFused+1)*sizeof(dlong), ogs->localFusedScatterIds);
  ogs->o_localFusedIds = device.malloc((totalIdsFusedSym+1)*sizeof(dlong), ogs->localFusedScatterIds);

  //set up the halo gatherScatter
  parallelNode_t *haloNodes = (parallelNode_t*) calloc(ogs->Nhalo+1,sizeof(parallelNode_t));

  cnt=0;
  for (dlong i=0;i<N;i++) {
    if (ids[i]==0) continue;

    if ((minRank[i]!=rank)||(maxRank[i]!=rank)) {
      haloNodes[cnt].localId = i;
      haloNodes[cnt].baseId  = ids[i];
      cnt++;
    }
  }

  // sort based on base ids (putting positive ids first) then local id
  qsort(haloNodes, ogs->Nhalo, sizeof(parallelNode_t), compareBaseId);

  ogs->NhaloGather=0;
  ogs->NhaloScatter=0;
  if (ogs->Nhalo) {
    haloNodes[0].newId = 0;
    int sign = (haloNodes[0].baseId > 0) ? 1 : -1;
    haloNodes[0].sign = sign;
    if (sign > 0) ogs->NhaloGather++;

    for (dlong i=1;i<ogs->Nhalo;i++) {
      if (abs(haloNodes[i].baseId)!=abs(haloNodes[i-1].baseId)) {
        sign = (haloNodes[i].baseId > 0) ? 1 : -1;
        ogs->NhaloScatter++;
        if (sign > 0) ogs->NhaloGather++;
      }

      haloNodes[i].newId = ogs->NhaloScatter;
      haloNodes[i].sign = sign;
    }
    ogs->NhaloScatter++;
  }

  // sort based on local ids
  qsort(haloNodes, ogs->Nhalo, sizeof(parallelNode_t), compareLocalId);

  //tally up how many nodes are being gathered to each gatherNode and
  //  map to a local ordering
  dlong *haloGatherCounts  = (dlong*) calloc(ogs->NhaloGather+1,sizeof(dlong));
  dlong *haloScatterCounts = (dlong*) calloc(ogs->NhaloScatter+1,sizeof(dlong));
  dlong *haloMap = (dlong*)  calloc(ogs->NhaloScatter+1,sizeof(dlong));
  hlong *haloIds = (hlong *) calloc(ogs->NhaloScatter+1,sizeof(hlong));
  hlong *haloIdsSym = (hlong *) calloc(ogs->NhaloScatter+1,sizeof(hlong));

  for (dlong i=0;i<ogs->NhaloScatter;i++) haloMap[i] = -1; //initialize map

  cnt = 0;
  cnt2 = ogs->NhaloGather;
  for (dlong i=0;i<ogs->Nhalo;i++) {
    dlong newId = haloNodes[i].newId; //get the ordered id

    if (haloMap[newId] == -1) {
      if (haloNodes[i].sign > 0)
        haloMap[newId] = cnt++;
      else
        haloMap[newId] = cnt2++;

      //record the base id of the gathered node
      haloIds[haloMap[newId]] = haloNodes[i].sign*abs(haloNodes[i].baseId);
      haloIdsSym[haloMap[newId]] = abs(haloNodes[i].baseId);
    }

    dlong gid = haloMap[newId];
    haloNodes[i].newId = gid;  //reorder
    haloScatterCounts[gid]++;  //tally
    if (haloNodes[i].baseId>0)
      haloGatherCounts[gid]++;  //tally
  }
  free(haloMap);

  ogs->haloGatherOffsets  = (dlong*) calloc(ogs->NhaloGather+1,sizeof(dlong));
  ogs->haloScatterOffsets = (dlong*) calloc(ogs->NhaloScatter+1,sizeof(dlong));
  for (dlong i=0;i<ogs->NhaloGather;i++) {
    ogs->haloGatherOffsets[i+1] = ogs->haloGatherOffsets[i] + haloGatherCounts[i];
    haloGatherCounts[i] = 0;
  }
  for (dlong i=0;i<ogs->NhaloScatter;i++) {
    ogs->haloScatterOffsets[i+1] = ogs->haloScatterOffsets[i] + haloScatterCounts[i];
    haloScatterCounts[i] = 0;
  }

  dlong totalHaloGather  = ogs->haloGatherOffsets[ogs->NhaloGather];
  dlong totalHaloScatter = ogs->haloScatterOffsets[ogs->NhaloScatter];

  ogs->haloGatherIds  = (dlong*) calloc(totalHaloGather+1,sizeof(dlong));
  ogs->haloScatterIds = (dlong*) calloc(totalHaloScatter+1,sizeof(dlong));
  for (dlong i=0;i<ogs->Nhalo;i++) {
    dlong gid = haloNodes[i].newId;

    dlong soffset = ogs->haloScatterOffsets[gid];
    int sindex  = haloScatterCounts[gid];
    ogs->haloScatterIds[soffset+sindex] = haloNodes[i].localId;
    haloScatterCounts[gid]++;

    if (haloNodes[i].baseId > 0) {
      dlong goffset = ogs->haloGatherOffsets[gid];
      int gindex  = haloGatherCounts[gid];
      ogs->haloGatherIds[goffset+gindex] = haloNodes[i].localId;
      haloGatherCounts[gid]++;
    }
  }
  free(haloGatherCounts);
  free(haloScatterCounts);

  ogs->o_haloGatherOffsets  = device.malloc((ogs->NhaloGather+1)*sizeof(dlong), ogs->haloGatherOffsets);
  ogs->o_haloScatterOffsets = device.malloc((ogs->NhaloScatter+1)*sizeof(dlong), ogs->haloScatterOffsets);

  ogs->o_haloGatherIds  = device.malloc((totalHaloGather+1)*sizeof(dlong), ogs->haloGatherIds);
  ogs->o_haloScatterIds = device.malloc((totalHaloScatter+1)*sizeof(dlong), ogs->haloScatterIds);

  free(haloNodes);

  //make a host gs handle
  ogs->Nlocal = ogs->NlocalScatter;
  ogs->Nhalo = ogs->NhaloScatter;
  ogs->gsh    = ogs::gsSetup(comm, ogs->Nhalo, haloIds, 0,0);
  ogs->gshSym = ogs::gsSetup(comm, ogs->Nhalo, haloIdsSym, 0,0);

  free(haloIds);
  free(haloIdsSym);

  free(minRank); free(maxRank);

  //total number of owned gathered nodes
  ogs->Ngather = ogs->NlocalGather+ogs->NhaloGather;

  ogs->hostBuf = NULL;
  ogs->haloBuf = NULL;
  ogs->hostBufSize = 0;

  return ogs;
}


void ogs_t::Free() {

  if (Nlocal) {
    free(localGatherOffsets);
    free(localScatterOffsets);
    free(localGatherIds);
    free(localScatterIds);
    o_localGatherOffsets.free();
    o_localScatterOffsets.free();
    o_localGatherIds.free();
    o_localScatterIds.free();
    Nlocal = 0;
  }

  if (Nhalo) {
    free(haloGatherOffsets);
    free(haloScatterOffsets);
    free(haloGatherIds);
    free(haloScatterIds);
    o_haloGatherOffsets.free();
    o_haloScatterOffsets.free();
    o_haloGatherIds.free();
    o_haloScatterIds.free();
    Nhalo = 0;
  }

  ogs::gsFree(gsh);

  ogs::Nrefs--;
  if (!ogs::Nrefs) ogs::freeKernels();
}

void ogs_t::reallocHostBuffer(size_t Nbytes) {
  if (Nhalo) {
    if (hostBufSize < Nhalo*Nbytes) {
      if (hostBufSize) free(hostBuf);
      hostBuf = (void *) malloc(Nhalo*Nbytes);
      hostBufSize = Nhalo*Nbytes;
    }
  }
}

void ogs_t::reallocOccaBuffer(size_t Nbytes) {
  if (Nhalo) {
    if (o_haloBuf.size() < Nhalo*Nbytes) {
      if (o_haloBuf.size()) o_haloBuf.free();
      haloBuf = occaHostMallocPinned(device, Nhalo*Nbytes, NULL,
                                     o_haloBuf, h_haloBuf);
    }
  }
}