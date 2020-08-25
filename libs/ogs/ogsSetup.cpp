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
#include "ogs/ogsKernels.hpp"

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

void setupRowBlocks(ogsData_t &A, platform_t &platform);

ogs_t *ogs_t::Setup(dlong N, hlong *ids, MPI_Comm &comm,
                    int verbose, platform_t& platform){

  ogs_t *ogs = new ogs_t(platform, comm);

  //Keep track of how many gs handles we've created, and
  // build kernels if this is the first
  if (!ogs::Nrefs) ogs::initKernels(platform);
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
  ogs->localGather.Nrows = 0;
  ogs->localScatter.Nrows = 0;
  if (ogs->Nlocal) {
    localNodes[0].newId = 0;
    int sign = (localNodes[0].baseId > 0) ? 1 : -1;
    localNodes[0].sign = sign;
    if (sign > 0) ogs->localGather.Nrows++;

    for (dlong i=1;i<ogs->Nlocal;i++) {
      if (abs(localNodes[i].baseId)!=abs(localNodes[i-1].baseId)) {
        sign = (localNodes[i].baseId > 0) ? 1 : -1;
        ogs->localScatter.Nrows++;
        if (sign > 0) ogs->localGather.Nrows++;
      }

      localNodes[i].newId = ogs->localScatter.Nrows;
      localNodes[i].sign = sign;
    }
    ogs->localScatter.Nrows++;
  }

  // sort back to local ids
  qsort(localNodes, ogs->Nlocal, sizeof(parallelNode_t), compareLocalId);

  //tally up how many nodes are being gathered to each gatherNode and
  //  map to a local ordering
  dlong *localGatherCounts  = (dlong*) calloc(ogs->localScatter.Nrows,sizeof(dlong));
  dlong *localScatterCounts = (dlong*) calloc(ogs->localScatter.Nrows,sizeof(dlong));

  dlong *localMap = (dlong*) calloc(ogs->localScatter.Nrows,sizeof(dlong));

  for (dlong i=0;i<ogs->localScatter.Nrows;i++) localMap[i] = -1; //initialize map

  cnt = 0;
  dlong cnt2 = ogs->localGather.Nrows;
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

  ogs->localGather.rowStarts  = (dlong*) calloc(ogs->localScatter.Nrows+1,sizeof(dlong));
  ogs->localScatter.rowStarts = (dlong*) calloc(ogs->localScatter.Nrows+1,sizeof(dlong));
  for (dlong i=0;i<ogs->localScatter.Nrows;i++) {
    ogs->localGather.rowStarts[i+1]  = ogs->localGather.rowStarts[i]  + localGatherCounts[i];
    ogs->localScatter.rowStarts[i+1] = ogs->localScatter.rowStarts[i] + localScatterCounts[i];

    //reset counters
    localScatterCounts[i] = 0;
    localGatherCounts[i] = 0;
  }

  ogs->localGather.nnz  = ogs->localGather.rowStarts[ogs->localGather.Nrows];
  ogs->localScatter.nnz = ogs->localScatter.rowStarts[ogs->localScatter.Nrows];

  ogs->localGather.colIds  = (dlong*) calloc(ogs->localGather.nnz+1,sizeof(dlong)); //extra entry so the occa buffer will actually exist
  ogs->localScatter.colIds = (dlong*) calloc(ogs->localScatter.nnz+1,sizeof(dlong)); //extra entry so the occa buffer will actually exist
  for (dlong i=0;i<ogs->Nlocal;i++) {
    dlong gid = localNodes[i].newId;

    dlong soffset = ogs->localScatter.rowStarts[gid];
    int sindex  = localScatterCounts[gid];
    ogs->localScatter.colIds[soffset+sindex] = localNodes[i].localId;
    localScatterCounts[gid]++;

    if (localNodes[i].baseId > 0) {
      dlong goffset = ogs->localGather.rowStarts[gid];
      int gindex  = localGatherCounts[gid];
      ogs->localGather.colIds[goffset+gindex] = localNodes[i].localId;
      localGatherCounts[gid]++;
    }
  }
  free(localGatherCounts);
  free(localScatterCounts);

  ogs->localGather.o_rowStarts  = platform.malloc((ogs->localScatter.Nrows+1)*sizeof(dlong), ogs->localGather.rowStarts);
  ogs->localScatter.o_rowStarts = platform.malloc((ogs->localScatter.Nrows+1)*sizeof(dlong), ogs->localScatter.rowStarts);

  ogs->localGather.o_colIds  = platform.malloc((ogs->localGather.nnz+1)*sizeof(dlong), ogs->localGather.colIds);
  ogs->localScatter.o_colIds = platform.malloc((ogs->localScatter.nnz+1)*sizeof(dlong), ogs->localScatter.colIds);

  //divide the list of colIds into roughly equal sized blocks so that each
  // threadblock loads approxiamtely an equal amount of data
  setupRowBlocks(ogs->localGather, platform);
  setupRowBlocks(ogs->localScatter, platform);

  free(localNodes);

  //make some compressed versions of the gather/scatter ids for the fused gs kernel
  ogs->fusedGather.Nrows=0;
  ogs->fusedScatter.Nrows=0;
  ogs->symGatherScatter.Nrows=0;

  ogs->fusedGather.nnz=0;
  ogs->fusedScatter.nnz=0;
  ogs->symGatherScatter.nnz=0;

  for (dlong n=0;n<ogs->localScatter.Nrows;n++) {
    int gatherCnt  = ogs->localGather.rowStarts[n+1] -ogs->localGather.rowStarts[n];
    int scatterCnt = ogs->localScatter.rowStarts[n+1]-ogs->localScatter.rowStarts[n];

    //only include this node if either the gather or scatter interact with mulitple nodes
    // otherwise the op is identity and ignored
    if ((gatherCnt>1)||(scatterCnt>1)) {
      ogs->fusedGather.Nrows++;
      ogs->fusedScatter.Nrows++;
      ogs->fusedGather.nnz  += gatherCnt;
      ogs->fusedScatter.nnz += scatterCnt;
    }

    //for the sym op only the scatter ids are used
    if (scatterCnt>1) {
      ogs->symGatherScatter.Nrows++;
      ogs->symGatherScatter.nnz += scatterCnt;
    }
  }

  ogs->fusedGather.rowStarts  = (dlong*) calloc(ogs->fusedScatter.Nrows+1,sizeof(dlong));
  ogs->fusedScatter.rowStarts = (dlong*) calloc(ogs->fusedScatter.Nrows+1,sizeof(dlong));
  ogs->symGatherScatter.rowStarts  = (dlong*) calloc(ogs->symGatherScatter.Nrows+1,sizeof(dlong));

  ogs->fusedGather.colIds  = (dlong*) calloc(ogs->fusedGather.nnz+1,sizeof(dlong));
  ogs->fusedScatter.colIds  = (dlong*) calloc(ogs->fusedScatter.nnz+1,sizeof(dlong));
  ogs->symGatherScatter.colIds  = (dlong*) calloc(ogs->symGatherScatter.nnz+1,sizeof(dlong));

  //reset counters
  ogs->fusedGather.Nrows=0;
  ogs->fusedScatter.Nrows=0;
  ogs->symGatherScatter.Nrows=0;

  ogs->fusedGather.nnz=0;
  ogs->fusedScatter.nnz=0;
  ogs->symGatherScatter.nnz=0;
  for (dlong n=0;n<ogs->localScatter.Nrows;n++) {
    int gatherCnt  = ogs->localGather.rowStarts[n+1] -ogs->localGather.rowStarts[n];
    int scatterCnt = ogs->localScatter.rowStarts[n+1]-ogs->localScatter.rowStarts[n];

    //only include this node if either the gather and scatter interact with mulitple nodes
    // otherwise the op is identity and ignored
    if ((gatherCnt>1)||(scatterCnt>1)) {
      ogs->fusedGather.Nrows++;
      ogs->fusedScatter.Nrows++;
      ogs->fusedGather.rowStarts[ogs->fusedGather.Nrows]   = gatherCnt  + ogs->fusedGather.rowStarts[ogs->fusedGather.Nrows-1];
      ogs->fusedScatter.rowStarts[ogs->fusedScatter.Nrows] = scatterCnt + ogs->fusedScatter.rowStarts[ogs->fusedScatter.Nrows-1];

      for (int i=ogs->localGather.rowStarts[n];i<ogs->localGather.rowStarts[n+1];i++)
        ogs->fusedGather.colIds[ogs->fusedGather.nnz++] = ogs->localGather.colIds[i];

      for (int i=ogs->localScatter.rowStarts[n];i<ogs->localScatter.rowStarts[n+1];i++)
        ogs->fusedScatter.colIds[ogs->fusedScatter.nnz++] = ogs->localScatter.colIds[i];
    }

    //for the sym op only the scatter ids are used
    if (scatterCnt>1) {
      ogs->symGatherScatter.Nrows++;
      ogs->symGatherScatter.rowStarts[ogs->symGatherScatter.Nrows] = scatterCnt + ogs->symGatherScatter.rowStarts[ogs->symGatherScatter.Nrows-1];

      for (int i=ogs->localScatter.rowStarts[n];i<ogs->localScatter.rowStarts[n+1];i++)
        ogs->symGatherScatter.colIds[ogs->symGatherScatter.nnz++] = ogs->localScatter.colIds[i];
    }
  }

  ogs->fusedGather.o_rowStarts  = platform.malloc((ogs->fusedScatter.Nrows+1)*sizeof(dlong), ogs->fusedGather.rowStarts);
  ogs->fusedScatter.o_rowStarts = platform.malloc((ogs->fusedScatter.Nrows+1)*sizeof(dlong), ogs->fusedScatter.rowStarts);
  ogs->symGatherScatter.o_rowStarts = platform.malloc((ogs->symGatherScatter.Nrows+1)*sizeof(dlong), ogs->symGatherScatter.rowStarts);

  ogs->fusedGather.o_colIds  = platform.malloc((ogs->fusedGather.nnz+1)*sizeof(dlong), ogs->fusedGather.colIds);
  ogs->fusedScatter.o_colIds = platform.malloc((ogs->fusedScatter.nnz+1)*sizeof(dlong), ogs->fusedScatter.colIds);
  ogs->symGatherScatter.o_colIds = platform.malloc((ogs->symGatherScatter.nnz+1)*sizeof(dlong), ogs->symGatherScatter.colIds);

  setupRowBlocks(ogs->fusedGather, platform);
  setupRowBlocks(ogs->fusedScatter, platform);
  setupRowBlocks(ogs->symGatherScatter, platform);

  //use the blocking from the fused scatter for the fusded gather as well
  if (ogs->fusedGather.blockRowStarts) free(ogs->fusedGather.blockRowStarts);
  ogs->fusedGather.o_blockRowStarts.free();
  ogs->fusedGather.NrowBlocks = ogs->fusedScatter.NrowBlocks;
  ogs->fusedGather.blockRowStarts = ogs->fusedScatter.blockRowStarts;
  ogs->fusedGather.o_blockRowStarts = ogs->fusedScatter.o_blockRowStarts;

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

  ogs->haloGather.Nrows = 0;
  ogs->haloScatter.Nrows = 0;

  if (ogs->Nhalo) {
    haloNodes[0].newId = 0;
    int sign = (haloNodes[0].baseId > 0) ? 1 : -1;
    haloNodes[0].sign = sign;
    if (sign > 0) ogs->haloGather.Nrows++;

    for (dlong i=1;i<ogs->Nhalo;i++) {
      if (abs(haloNodes[i].baseId)!=abs(haloNodes[i-1].baseId)) {
        sign = (haloNodes[i].baseId > 0) ? 1 : -1;
        ogs->haloScatter.Nrows++;
        if (sign > 0) ogs->haloGather.Nrows++;
      }

      haloNodes[i].newId = ogs->haloScatter.Nrows;
      haloNodes[i].sign = sign;
    }
    ogs->haloScatter.Nrows++;
  }

  // sort based on local ids
  qsort(haloNodes, ogs->Nhalo, sizeof(parallelNode_t), compareLocalId);

  //tally up how many nodes are being gathered to each gatherNode and
  //  map to a local ordering
  dlong *haloGatherCounts  = (dlong*) calloc(ogs->haloGather.Nrows+1,sizeof(dlong));
  dlong *haloScatterCounts = (dlong*) calloc(ogs->haloScatter.Nrows+1,sizeof(dlong));
  dlong *haloMap = (dlong*)  calloc(ogs->haloScatter.Nrows+1,sizeof(dlong));
  hlong *haloIds = (hlong *) calloc(ogs->haloScatter.Nrows+1,sizeof(hlong));
  hlong *haloIdsSym = (hlong *) calloc(ogs->haloScatter.Nrows+1,sizeof(hlong));

  for (dlong i=0;i<ogs->haloScatter.Nrows;i++) haloMap[i] = -1; //initialize map

  cnt = 0;
  cnt2 = ogs->haloGather.Nrows;
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

  ogs->haloGather.rowStarts  = (dlong*) calloc(ogs->haloGather.Nrows+1,sizeof(dlong));
  ogs->haloScatter.rowStarts = (dlong*) calloc(ogs->haloScatter.Nrows+1,sizeof(dlong));
  for (dlong i=0;i<ogs->haloGather.Nrows;i++) {
    ogs->haloGather.rowStarts[i+1] = ogs->haloGather.rowStarts[i] + haloGatherCounts[i];
    haloGatherCounts[i] = 0;
  }
  for (dlong i=0;i<ogs->haloScatter.Nrows;i++) {
    ogs->haloScatter.rowStarts[i+1] = ogs->haloScatter.rowStarts[i] + haloScatterCounts[i];
    haloScatterCounts[i] = 0;
  }

  ogs->haloGather.nnz  = ogs->haloGather.rowStarts[ogs->haloGather.Nrows];
  ogs->haloScatter.nnz = ogs->haloScatter.rowStarts[ogs->haloScatter.Nrows];

  ogs->haloGather.colIds  = (dlong*) calloc(ogs->haloGather.nnz+1,sizeof(dlong));
  ogs->haloScatter.colIds = (dlong*) calloc(ogs->haloScatter.nnz+1,sizeof(dlong));
  for (dlong i=0;i<ogs->Nhalo;i++) {
    dlong gid = haloNodes[i].newId;

    dlong soffset = ogs->haloScatter.rowStarts[gid];
    int sindex  = haloScatterCounts[gid];
    ogs->haloScatter.colIds[soffset+sindex] = haloNodes[i].localId;
    haloScatterCounts[gid]++;

    if (haloNodes[i].baseId > 0) {
      dlong goffset = ogs->haloGather.rowStarts[gid];
      int gindex  = haloGatherCounts[gid];
      ogs->haloGather.colIds[goffset+gindex] = haloNodes[i].localId;
      haloGatherCounts[gid]++;
    }
  }
  free(haloGatherCounts);
  free(haloScatterCounts);

  ogs->haloGather.o_rowStarts  = platform.malloc((ogs->haloGather.Nrows+1)*sizeof(dlong), ogs->haloGather.rowStarts);
  ogs->haloScatter.o_rowStarts = platform.malloc((ogs->haloScatter.Nrows+1)*sizeof(dlong), ogs->haloScatter.rowStarts);

  ogs->haloGather.o_colIds  = platform.malloc((ogs->haloGather.nnz+1)*sizeof(dlong), ogs->haloGather.colIds);
  ogs->haloScatter.o_colIds = platform.malloc((ogs->haloScatter.nnz+1)*sizeof(dlong), ogs->haloScatter.colIds);

  setupRowBlocks(ogs->haloGather, platform);
  setupRowBlocks(ogs->haloScatter, platform);

  free(haloNodes);

  //make a host gs handle
  ogs->Nlocal = ogs->localScatter.Nrows;
  ogs->Nhalo = ogs->haloScatter.Nrows;
  ogs->gsh    = ogs::gsSetup(comm, ogs->Nhalo, haloIds, 0,0);
  ogs->gshSym = ogs::gsSetup(comm, ogs->Nhalo, haloIdsSym, 0,0);

  free(haloIds);
  free(haloIdsSym);

  free(minRank); free(maxRank);

  //total number of owned gathered nodes
  ogs->Ngather = ogs->localGather.Nrows+ogs->haloGather.Nrows;

  hlong NgatherLocal = (hlong) ogs->Ngather;
  MPI_Allreduce(&NgatherLocal, &(ogs->NgatherGlobal), 1, MPI_HLONG, MPI_SUM, comm);

  ogs->hostBuf = nullptr;
  ogs->haloBuf = nullptr;
  ogs->hostBufSize = 0;

  return ogs;
}

void ogs_t::Free() {

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
      haloBuf = platform.hostMalloc(Nhalo*Nbytes, nullptr, h_haloBuf);
      o_haloBuf = platform.malloc(Nhalo*Nbytes);
    }
  }
}

void setupRowBlocks(ogsData_t &A, platform_t &platform) {

  dlong blockSum=0;
  A.NrowBlocks=0;
  if (A.Nrows) A.NrowBlocks++;
  for (dlong i=0;i<A.Nrows;i++) {
    dlong rowSize = A.rowStarts[i+1]-A.rowStarts[i];

    if (rowSize > ogs::gatherNodesPerBlock) {
      //this row is pathalogically big. We can't currently run this
      stringstream ss;
      ss << "Multiplicity of global node id: " << i << "in ogsSetup is too large.";
      LIBP_ABORT(ss.str())
    }

    if (blockSum+rowSize > ogs::gatherNodesPerBlock) { //adding this row will exceed the nnz per block
      A.NrowBlocks++; //count the previous block
      blockSum=rowSize; //start a new row block
    } else {
      blockSum+=rowSize; //add this row to the block
    }
  }

  A.blockRowStarts  = (dlong*) calloc(A.NrowBlocks+1,sizeof(dlong));

  blockSum=0;
  A.NrowBlocks=0;
  if (A.Nrows) A.NrowBlocks++;
  for (dlong i=0;i<A.Nrows;i++) {
    dlong rowSize = A.rowStarts[i+1]-A.rowStarts[i];

    if (blockSum+rowSize > ogs::gatherNodesPerBlock) { //adding this row will exceed the nnz per block
      A.blockRowStarts[A.NrowBlocks++] = i; //mark the previous block
      blockSum=rowSize; //start a new row block
    } else {
      blockSum+=rowSize; //add this row to the block
    }
  }
  A.blockRowStarts[A.NrowBlocks] = A.Nrows;

  A.o_blockRowStarts = platform.malloc((A.NrowBlocks+1)*sizeof(dlong), A.blockRowStarts);
}