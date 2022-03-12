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

#include "ogs.hpp"
#include "ogs/ogsUtils.hpp"
#include "ogs/ogsOperator.hpp"
#include "ogs/ogsExchange.hpp"
#include "timer.hpp"

#ifdef GLIBCXX_PARALLEL
#include <parallel/algorithm>
using __gnu_parallel::sort;
#else
using std::sort;
#endif

namespace libp {

namespace ogs {

void ogs_t::Setup(const dlong _N,
                  memory<hlong> ids,
                  comm_t _comm,
                  const Kind _kind,
                  const Method method,
                  const bool _unique,
                  const bool verbose,
                  platform_t& _platform){
  ogsBase_t::Setup(_N, ids, _comm, _kind, method, _unique, verbose, _platform);
}

void halo_t::Setup(const dlong _N,
                  memory<hlong> ids,
                  comm_t _comm,
                  const Method method,
                  const bool verbose,
                  platform_t& _platform){
  ogsBase_t::Setup(_N, ids, _comm, Halo, method, false, verbose, _platform);

  Nhalo = NhaloT - NhaloP; //number of extra recieved nodes
}

/********************************
 * Setup
 ********************************/
void ogsBase_t::Setup(const dlong _N,
                      memory<hlong> ids,
                      comm_t _comm,
                      const Kind _kind,
                      const Method method,
                      const bool _unique,
                      const bool verbose,
                      platform_t& _platform){

  //release resources if this ogs was setup before
  Free();

  timePoint_t start = Time();

  platform = _platform;

  if (!dataStream.isInitialized())
      dataStream = platform.device.createStream();

  N = _N;
  comm = _comm;
  kind = _kind;
  unique = _unique;

  int rank, size;
  rank = comm.rank();
  size = comm.size();

  //sanity check options
  LIBP_ABORT("Invalid ogs setup requested",
             (kind==Unsigned && unique==true)
              || (kind==Halo && unique==true));

  //count how many ids are non-zero
  dlong Nids=0;
  for (dlong n=0;n<N;n++)
    if (ids[n]!=0) Nids++;

  // make list of nodes
  memory<parallelNode_t> nodes(Nids);

  //fill the data (squeezing out zero ids)
  Nids=0;
  for (dlong n=0;n<N;n++) {
    if (ids[n]!=0) {
      nodes[Nids].localId = Nids; //record a compressed id first (useful for ordering)
      nodes[Nids].baseId = (kind==Unsigned) ?
                            abs(ids[n]) : ids[n]; //record global id
      nodes[Nids].rank = rank;
      nodes[Nids].destRank = abs(ids[n]) % size;
      Nids++;
    }
  }

  //flag which nodes are shared via MPI
  FindSharedNodes(Nids, nodes, verbose);

  //Index the local and halo baseIds on this rank and
  // construct sharedNodes which contains all the info
  // we need to setup the MPI exchange.
  dlong Nshared=0;
  memory<parallelNode_t> sharedNodes;
  ConstructSharedNodes(Nids, nodes, Nshared, sharedNodes);

  Nids=0;
  for (dlong n=0;n<N;n++) {
    if (ids[n]!=0) {
      nodes[Nids].localId = n; //record the real id now

      //if we altered the signs of ids, write them back
      if (unique)
        ids[n] = nodes[Nids].baseId;

      Nids++;
    }
  }

  //setup local gather operators
  if (kind==Signed)
    LocalSignedSetup(Nids, nodes);
  else if (kind==Unsigned)
    LocalUnsignedSetup(Nids, nodes);
  else
    LocalHaloSetup(Nids, nodes);

  //with that, we're done with the local nodes list
  nodes.free();

  // At this point, we've setup gs operators to gather/scatter the purely local nodes,
  // and gather/scatter the shared halo nodes to/from a coalesced ordering. We now
  // need gs operators to scatter/gather the coalesced halo nodes to/from the expected
  // orderings for MPI communications.

  if (method == AllToAll) {
    exchange = std::shared_ptr<ogsExchange_t>(
                  new ogsAllToAll_t(Nshared, sharedNodes,
                                    *gatherHalo, dataStream,
                                    comm, platform));
  } else if (method == Pairwise) {
    exchange = std::shared_ptr<ogsExchange_t>(
                  new ogsPairwise_t(Nshared, sharedNodes,
                                    *gatherHalo, dataStream,
                                    comm, platform));
  } else if (method == CrystalRouter) {
    exchange = std::shared_ptr<ogsExchange_t>(
                  new ogsCrystalRouter_t(Nshared, sharedNodes,
                                         *gatherHalo, dataStream,
                                         comm, platform));
  } else { //Auto
    exchange = std::shared_ptr<ogsExchange_t>(
                  AutoSetup(Nshared, sharedNodes,
                            *gatherHalo, comm,
                            platform, verbose));
  }

  timePoint_t end = GlobalPlatformTime(platform);
  double elapsedTime = ElapsedTime(start, end);

  if (!rank && verbose) {
    std::cout << "ogs Setup Time: " << elapsedTime << " seconds." << std::endl;
  }
}

void ogsBase_t::FindSharedNodes(const dlong Nids,
                                memory<parallelNode_t> &nodes,
                                const int verbose){

  int rank, size;
  rank = comm.rank();
  size = comm.size();

  memory<int> sendCounts(size,0);
  memory<int> recvCounts(size);
  memory<int> sendOffsets(size+1);
  memory<int> recvOffsets(size+1);

  //count number of ids we're sending
  for (dlong n=0;n<Nids;n++) {
    sendCounts[nodes[n].destRank]++;
  }

  comm.Alltoall(sendCounts, recvCounts);

  sendOffsets[0] = 0;
  recvOffsets[0] = 0;
  for (int r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];

    //reset counter
    sendCounts[r] = 0;
  }

  //write a send ordering into newIds
  for (dlong n=0;n<Nids;n++) {
    const int r = nodes[n].destRank;
    nodes[n].newId = sendOffsets[r]+sendCounts[r]++;
  }

  // permute the list to send ordering
  permute(Nids, nodes, [](const parallelNode_t& a) { return a.newId; } );

  dlong recvN = recvOffsets[size]; //total ids to recv

  memory<parallelNode_t> recvNodes(recvN);

  //Send all the nodes to their destination rank.
  comm.Alltoallv(    nodes, sendCounts, sendOffsets,
                 recvNodes, recvCounts, recvOffsets);

  //remember this ordering
  for (dlong n=0;n<recvN;n++) {
    recvNodes[n].newId = n;
  }

  // sort based on base ids
  sort(recvNodes.ptr(), recvNodes.ptr()+recvN,
       [](const parallelNode_t& a, const parallelNode_t& b) {
         return abs(a.baseId) < abs(b.baseId);
       });

  // We now have a collection of nodes associated with some subset of all global Ids
  // Our list is sorted by baseId to group nodes with the same globalId together
  // We now want to flag which nodes are shared via MPI

  int is_unique=1;

  dlong Nshared=0;

  dlong start=0;
  for (dlong n=0;n<recvN;n++) {
    if (n==recvN-1 || abs(recvNodes[n].baseId)!=abs(recvNodes[n+1].baseId)) {
      dlong end = n+1;

      int positiveCount=0;
      if (unique) {
        //Make a single node from each baseId group the sole positive node
        const hlong baseId = abs(recvNodes[start].baseId);

        //pick a random node in this group
        const int m = (rand() % (end-start));

        for (int i=start;i<end;i++)
          recvNodes[i].baseId = -baseId;

        recvNodes[start+m].baseId = baseId;
        positiveCount=1;
      } else {
        //count how many postive baseIds there are in this group
        for (int i=start;i<end;i++)
          if (recvNodes[i].baseId>0) positiveCount++;

        //if we didnt find a sole positive baseId, the gather is not well-defined
        if (positiveCount!=1) is_unique=0;
      }

      // When making a halo excahnge, check that we have a leading positive id
      LIBP_ABORT("Found " << positiveCount << " positive Ids for baseId: "
                 << abs(recvNodes[start].baseId)<< ".",
                 kind==Halo && positiveCount!=1);

      //determine if this node is shared via MPI,
      int shared=1;
      const int r = recvNodes[start].rank;
      for (int i=start+1;i<end;i++) {
        if (recvNodes[i].rank != r) {
          shared=2;
          Nshared++;
          break;
        }
      }

      //set shared flag.
      for (int i=start;i<end;i++) {
        recvNodes[i].sign = shared;
      }

      //set new baseId group start point
      start=n+1;
    }
  }

  //shared the unique node check so we know if the gather operation is well-defined
  comm.Allreduce(is_unique, comm_t::Min);
  gather_defined = (is_unique==1);

  hlong Nshared_global = Nshared;
  comm.Reduce(Nshared_global, 0);
  if (!rank && verbose) {
    std::cout << "ogs Setup: " << Nshared_global << " unique labels shared." << std::endl;
  }

  //at this point each collection of baseIds either has all nodes have
  // sign = 1, meaning all the nodes with this baseId are on the
  // same rank, or have sign=2, meaning that baseId must be communicated

  // permute recv nodes back to recv'd ordering
  permute(recvN, recvNodes, [](const parallelNode_t& a) { return a.newId; } );

  //Return all the nodes to their origin rank.
  comm.Alltoallv(recvNodes, recvCounts, recvOffsets,
                     nodes, sendCounts, sendOffsets);
}

void ogsBase_t::ConstructSharedNodes(const dlong Nids,
                                     memory<parallelNode_t> &nodes,
                                     dlong &Nshared,
                                     memory<parallelNode_t> &sharedNodes) {

  int size = comm.size();

  // sort based on abs(baseId)
  sort(nodes.ptr(), nodes.ptr()+Nids,
       [](const parallelNode_t& a, const parallelNode_t& b) {
         if(abs(a.baseId) < abs(b.baseId)) return true; //group by abs(baseId)
         if(abs(a.baseId) > abs(b.baseId)) return false;

         return a.baseId > b.baseId; //positive ids on a rank first
       });

  //count how many unique global Ids we have on this rank
  // and flag baseId groups that have a positive baseId somewhere on this rank
  dlong NbaseIds=0;
  NlocalT=0; NlocalP=0;
  NhaloT=0; NhaloP=0;
  dlong start=0;
  for (dlong n=0;n<Nids;n++) {
    if (n==Nids-1 || abs(nodes[n].baseId)!=abs(nodes[n+1].baseId)) {
      dlong end = n+1;

      //if there's no leading postive id, flag this baseId group as negative
      int sign = abs(nodes[start].sign);
      if (nodes[start].baseId<0) {
        sign = -sign;
        for (int i=start;i<end;i++) {
          nodes[i].sign = sign;
        }
      }

      //count the positive/negative local and halo gather nodes
      if (abs(sign)==1) {
        NlocalT++;
        if (sign==1) NlocalP++;
      } else {
        NhaloT++;
        if (sign==2) NhaloP++;
      }

      //record the new ordering
      for (int i=start;i<end;i++) {
        nodes[i].newId=NbaseIds;
      }

      NbaseIds++;
      start = end;
    }
  }

  //total number of positive owned gathered nodes
  Ngather = NlocalP+NhaloP;

  //global total
  NgatherGlobal = Ngather;
  comm.Allreduce(NgatherGlobal);

  //extract the leading node from each shared baseId
  memory<parallelNode_t> sendSharedNodes(NhaloT);

  NhaloT=0;
  for (dlong n=0;n<Nids;n++) {
    if (n==0 || abs(nodes[n].baseId)!=abs(nodes[n-1].baseId)) {
      if (abs(nodes[n].sign)==2) {
        sendSharedNodes[NhaloT++] = nodes[n];
      }
    }
  }

  // permute the list back to local id ordering
  permute(Nids, nodes, [](const parallelNode_t& a) { return a.localId; } );

  // Use the newId index to reorder the baseId groups based on
  // the order we encouter them in their original ordering.
  memory<dlong> indexMap(NbaseIds, -1);

  dlong localCntN = 0, localCntT = NlocalP;  //start point for local gather nodes
  dlong haloCntN  = 0, haloCntT  = NhaloP;   //start point for halo gather nodes
  for (dlong n=0;n<Nids;n++) {
    const dlong newId = nodes[n].newId; //get the new baseId group id

    //record a new index if we've not encoutered this baseId group before
    if (indexMap[newId]==-1) {
      if        (nodes[n].sign== 1) {
        indexMap[newId] = localCntN++;
      } else if (nodes[n].sign==-1) {
        indexMap[newId] = localCntT++;
      } else if (nodes[n].sign== 2) {
        indexMap[newId] = haloCntN++;
      } else { //nodes[n].sign==-2
        indexMap[newId] = haloCntT++;
      }
    }

    const dlong gid = indexMap[newId];
    nodes[n].newId = gid; //reorder
  }

  //re-order the shared node list
  for (dlong n=0;n<NhaloT;n++) {
    const dlong newId = sendSharedNodes[n].newId; //get the new baseId group id
    const dlong gid = indexMap[newId];
    sendSharedNodes[n].localId = gid; //reorder the localId to the compressed order
  }

  indexMap.free();

  memory<int> sendCounts(size,0);
  memory<int> recvCounts(size);
  memory<int> sendOffsets(size+1);
  memory<int> recvOffsets(size+1);

  // sort based on destination rank
  sort(sendSharedNodes.ptr(), sendSharedNodes.ptr()+NhaloT,
       [](const parallelNode_t& a, const parallelNode_t& b) {
         return a.destRank < b.destRank;
       });

  //count number of ids we're sending
  for (dlong n=0;n<NhaloT;n++) {
    sendCounts[sendSharedNodes[n].destRank]++;
  }

  comm.Alltoall(sendCounts, recvCounts);

  sendOffsets[0] = 0;
  recvOffsets[0] = 0;
  for (int r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }
  dlong recvN = recvOffsets[size]; //total ids to recv

  memory<parallelNode_t> recvSharedNodes(recvN);

  //Send all the nodes to their destination rank.
  comm.Alltoallv(sendSharedNodes, sendCounts, sendOffsets,
                 recvSharedNodes, recvCounts, recvOffsets);

  //free up some space
  sendSharedNodes.free();
  sendCounts.free();
  recvCounts.free();
  sendOffsets.free();
  recvOffsets.free();

  // sort based on base ids
  sort(recvSharedNodes.ptr(), recvSharedNodes.ptr()+recvN,
       [](const parallelNode_t& a, const parallelNode_t& b) {
         return abs(a.baseId) < abs(b.baseId);
       });

  //count number of shared nodes we will be sending
  memory<int> sharedSendCounts(size,0);
  memory<int> sharedRecvCounts(size);
  memory<int> sharedSendOffsets(size+1);
  memory<int> sharedRecvOffsets(size+1);

  start=0;
  for (dlong n=0;n<recvN;n++) {
    if (n==recvN-1 || abs(recvSharedNodes[n].baseId)!=abs(recvSharedNodes[n+1].baseId)) {
      dlong end = n+1;

      for (int i=start;i<end;i++) {
        //We'll be sending all the shared nodes to each rank involved
        sharedSendCounts[recvSharedNodes[i].rank] += end-start-1;
      }

      //set new baseId group start point
      start=n+1;
    }
  }

  // Each rank has a set of shared global Ids and for each global id that
  // rank knows what MPI ranks participate in gathering. We now send this
  // information to the involved ranks.

  //share counts
  comm.Alltoall(sharedSendCounts, sharedRecvCounts);

  //cumulative sum
  sharedSendOffsets[0] = 0;
  sharedRecvOffsets[0] = 0;
  for (int r=0;r<size;r++) {
    sharedSendOffsets[r+1] = sharedSendOffsets[r]+sharedSendCounts[r];
    sharedRecvOffsets[r+1] = sharedRecvOffsets[r]+sharedRecvCounts[r];
  }

  //make a send buffer
  memory<parallelNode_t> sharedSendNodes(sharedSendOffsets[size]);

  //reset sendCounts
  for (int r=0;r<size;r++) sharedSendCounts[r]=0;

  start=0;
  for (dlong n=0;n<recvN;n++) {
    if (n==recvN-1 || abs(recvSharedNodes[n].baseId)!=abs(recvSharedNodes[n+1].baseId)) {
      dlong end = n+1;

      //build the node list to send
      for (int i=start;i<end;i++) {
        const int r = recvSharedNodes[i].rank;
        const dlong id = recvSharedNodes[i].localId;
        const int sign = recvSharedNodes[i].sign;

        int sid = sharedSendCounts[r]+sharedSendOffsets[r];
        for (int j=start;j<end;j++) {
          if (j==i) continue; //dont bother sending this rank's own node
          sharedSendNodes[sid] = recvSharedNodes[j];
          sharedSendNodes[sid].newId = id;
          sharedSendNodes[sid].sign = sign;
          sid++;
        }
        sharedSendCounts[r] += end-start-1;
      }

      //set new baseId group start point
      start=n+1;
    }
  }
  recvSharedNodes.free();

  //make sharedNodes to hold the exchange data we recv
  Nshared = sharedRecvOffsets[size];
  sharedNodes = memory<parallelNode_t>(Nshared);

  //Share all the gathering info
  comm.Alltoallv(sharedSendNodes, sharedSendCounts, sharedSendOffsets,
                     sharedNodes, sharedRecvCounts, sharedRecvOffsets);
}

//Make local and halo gather operators using nodes list
void ogsBase_t::LocalSignedSetup(const dlong Nids, memory<parallelNode_t> &nodes){

  gatherLocal = std::make_shared<ogsOperator_t>(platform);
  gatherHalo  = std::make_shared<ogsOperator_t>(platform);

  gatherLocal->kind = Signed;
  gatherHalo->kind = Signed;

  gatherLocal->Ncols = N;
  gatherHalo->Ncols = N;

  gatherLocal->NrowsN = NlocalP;
  gatherLocal->NrowsT = NlocalT;
  gatherHalo->NrowsN = NhaloP;
  gatherHalo->NrowsT = NhaloT;

  //tally up how many nodes are being gathered to each gatherNode and
  //  map to a local ordering
  memory<dlong> localGatherNCounts(gatherLocal->NrowsT,0);
  memory<dlong> localGatherTCounts(gatherLocal->NrowsT,0);
  memory<dlong> haloGatherNCounts(gatherHalo->NrowsT,0);
  memory<dlong> haloGatherTCounts(gatherHalo->NrowsT,0);

  for (dlong i=0;i<Nids;i++) {
    const dlong gid = nodes[i].newId; //re-mapped baseId on this rank

    if (abs(nodes[i].sign)==1) { //local
      if (nodes[i].baseId>0) localGatherNCounts[gid]++;  //tally
      localGatherTCounts[gid]++;  //tally
    } else { //halo
      if (nodes[i].baseId>0) haloGatherNCounts[gid]++;  //tally
      haloGatherTCounts[gid]++;  //tally
    }
  }

  //make local row offsets
  gatherLocal->rowStartsN.malloc(gatherLocal->NrowsT+1);
  gatherLocal->rowStartsT.malloc(gatherLocal->NrowsT+1);
  gatherLocal->rowStartsN[0] = 0;
  gatherLocal->rowStartsT[0] = 0;
  for (dlong i=0;i<gatherLocal->NrowsT;i++) {
    gatherLocal->rowStartsN[i+1] = gatherLocal->rowStartsN[i] + localGatherNCounts[i];
    gatherLocal->rowStartsT[i+1] = gatherLocal->rowStartsT[i] + localGatherTCounts[i];
    localGatherNCounts[i] = 0; //reset counters
    localGatherTCounts[i] = 0; //reset counters
  }
  gatherLocal->nnzN = gatherLocal->rowStartsN[gatherLocal->NrowsT];
  gatherLocal->nnzT = gatherLocal->rowStartsT[gatherLocal->NrowsT];
  gatherLocal->colIdsN.malloc(gatherLocal->nnzN);
  gatherLocal->colIdsT.malloc(gatherLocal->nnzT);

  //make halo row offsets
  gatherHalo->rowStartsN.malloc(gatherHalo->NrowsT+1);
  gatherHalo->rowStartsT.malloc(gatherHalo->NrowsT+1);
  gatherHalo->rowStartsN[0] = 0;
  gatherHalo->rowStartsT[0] = 0;
  for (dlong i=0;i<gatherHalo->NrowsT;i++) {
    gatherHalo->rowStartsN[i+1] = gatherHalo->rowStartsN[i] + haloGatherNCounts[i];
    gatherHalo->rowStartsT[i+1] = gatherHalo->rowStartsT[i] + haloGatherTCounts[i];
    haloGatherNCounts[i] = 0;
    haloGatherTCounts[i] = 0;
  }
  gatherHalo->nnzN = gatherHalo->rowStartsN[gatherHalo->NrowsT];
  gatherHalo->nnzT = gatherHalo->rowStartsT[gatherHalo->NrowsT];
  gatherHalo->colIdsN.malloc(gatherHalo->nnzN);
  gatherHalo->colIdsT.malloc(gatherHalo->nnzT);


  for (dlong i=0;i<Nids;i++) {
    const dlong gid = nodes[i].newId;

    if (abs(nodes[i].sign)==1) { //local gather group
      if (nodes[i].baseId>0) {
        const dlong soffset = gatherLocal->rowStartsN[gid];
        const int sindex  = localGatherNCounts[gid];
        gatherLocal->colIdsN[soffset+sindex] = nodes[i].localId;
        localGatherNCounts[gid]++;
      }
      const dlong soffset = gatherLocal->rowStartsT[gid];
      const int sindex  = localGatherTCounts[gid];
      gatherLocal->colIdsT[soffset+sindex] = nodes[i].localId;
      localGatherTCounts[gid]++;
    } else {
      if (nodes[i].baseId>0) {
        const dlong soffset = gatherHalo->rowStartsN[gid];
        const int sindex  = haloGatherNCounts[gid];
        gatherHalo->colIdsN[soffset+sindex] = nodes[i].localId;
        haloGatherNCounts[gid]++;
      }
      const dlong soffset = gatherHalo->rowStartsT[gid];
      const int sindex  = haloGatherTCounts[gid];
      gatherHalo->colIdsT[soffset+sindex] = nodes[i].localId;
      haloGatherTCounts[gid]++;
    }
  }
  localGatherNCounts.free();
  localGatherTCounts.free();
  haloGatherNCounts.free();
  haloGatherTCounts.free();

  gatherLocal->o_rowStartsN = platform.malloc(gatherLocal->rowStartsN);
  gatherLocal->o_rowStartsT = platform.malloc(gatherLocal->rowStartsT);
  gatherLocal->o_colIdsN = platform.malloc(gatherLocal->colIdsN);
  gatherLocal->o_colIdsT = platform.malloc(gatherLocal->colIdsT);

  gatherHalo->o_rowStartsN = platform.malloc(gatherHalo->rowStartsN);
  gatherHalo->o_rowStartsT = platform.malloc(gatherHalo->rowStartsT);
  gatherHalo->o_colIdsN = platform.malloc(gatherHalo->colIdsN);
  gatherHalo->o_colIdsT = platform.malloc(gatherHalo->colIdsT);

  //divide the list of colIds into roughly equal sized blocks so that each
  // threadblock loads approximately an equal amount of data
  gatherLocal->setupRowBlocks();
  gatherHalo->setupRowBlocks();
}

//Make local and halo gather operators using nodes list
void ogsBase_t::LocalUnsignedSetup(const dlong Nids, memory<parallelNode_t> &nodes){

  gatherLocal = std::make_shared<ogsOperator_t>(platform);
  gatherHalo  = std::make_shared<ogsOperator_t>(platform);

  gatherLocal->kind = Unsigned;
  gatherHalo->kind = Unsigned;

  gatherLocal->Ncols = N;
  gatherHalo->Ncols = N;

  gatherLocal->NrowsN = NlocalP;
  gatherLocal->NrowsT = NlocalT;
  gatherHalo->NrowsN = NhaloP;
  gatherHalo->NrowsT = NhaloT;

  //tally up how many nodes are being gathered to each gatherNode and
  //  map to a local ordering
  memory<dlong> localGatherTCounts(gatherLocal->NrowsT,0);
  memory<dlong> haloGatherTCounts(gatherHalo->NrowsT,0);

  for (dlong i=0;i<Nids;i++) {
    const dlong gid = nodes[i].newId; //re-mapped baseId on this rank

    if (abs(nodes[i].sign)==1) { //local
      localGatherTCounts[gid]++;  //tally
    } else { //halo
      haloGatherTCounts[gid]++;  //tally
    }
  }

  //make local row offsets
  gatherLocal->rowStartsT.malloc(gatherLocal->NrowsT+1);
  gatherLocal->rowStartsN = gatherLocal->rowStartsT;
  gatherLocal->rowStartsT[0] = 0;
  for (dlong i=0;i<gatherLocal->NrowsT;i++) {
    gatherLocal->rowStartsT[i+1] = gatherLocal->rowStartsT[i] + localGatherTCounts[i];
    localGatherTCounts[i] = 0; //reset counters
  }
  gatherLocal->nnzT = gatherLocal->rowStartsT[gatherLocal->NrowsT];
  gatherLocal->nnzN = gatherLocal->nnzT;
  gatherLocal->colIdsT.malloc(gatherLocal->nnzT);
  gatherLocal->colIdsN = gatherLocal->colIdsT;

  //make halo row offsets
  gatherHalo->rowStartsT.malloc(gatherHalo->NrowsT+1);
  gatherHalo->rowStartsN = gatherHalo->rowStartsT;
  gatherHalo->rowStartsT[0] = 0;
  for (dlong i=0;i<gatherHalo->NrowsT;i++) {
    gatherHalo->rowStartsT[i+1] = gatherHalo->rowStartsT[i] + haloGatherTCounts[i];
    haloGatherTCounts[i] = 0;
  }
  gatherHalo->nnzT = gatherHalo->rowStartsT[gatherHalo->NrowsT];
  gatherHalo->nnzN = gatherHalo->nnzT;
  gatherHalo->colIdsT.malloc(gatherHalo->nnzT);
  gatherHalo->colIdsN = gatherHalo->colIdsT;


  for (dlong i=0;i<Nids;i++) {
    const dlong gid = nodes[i].newId;

    if (abs(nodes[i].sign)==1) { //local gather group
      const dlong soffset = gatherLocal->rowStartsT[gid];
      const int sindex  = localGatherTCounts[gid];
      gatherLocal->colIdsT[soffset+sindex] = nodes[i].localId;
      localGatherTCounts[gid]++;
    } else {
      const dlong soffset = gatherHalo->rowStartsT[gid];
      const int sindex  = haloGatherTCounts[gid];
      gatherHalo->colIdsT[soffset+sindex] = nodes[i].localId;
      haloGatherTCounts[gid]++;
    }
  }
  localGatherTCounts.free();
  haloGatherTCounts.free();

  gatherLocal->o_rowStartsT = platform.malloc(gatherLocal->rowStartsT);
  gatherLocal->o_rowStartsN = gatherLocal->o_rowStartsT;
  gatherLocal->o_colIdsT = platform.malloc(gatherLocal->colIdsT);
  gatherLocal->o_colIdsN = gatherLocal->o_colIdsT;

  gatherHalo->o_rowStartsT = platform.malloc(gatherHalo->rowStartsT);
  gatherHalo->o_rowStartsN = gatherHalo->o_rowStartsT;
  gatherHalo->o_colIdsT = platform.malloc(gatherHalo->colIdsT);
  gatherHalo->o_colIdsN = gatherHalo->o_colIdsT;

  //divide the list of colIds into roughly equal sized blocks so that each
  // threadblock loads approximately an equal amount of data
  gatherLocal->setupRowBlocks();
  gatherHalo->setupRowBlocks();
}

//Make local and halo gather operators using nodes list
void ogsBase_t::LocalHaloSetup(const dlong Nids, memory<parallelNode_t> &nodes){

  gatherHalo  = std::make_shared<ogsOperator_t>(platform);
  gatherHalo->kind = Signed;

  gatherHalo->Ncols = N;

  gatherHalo->NrowsN = NhaloP;
  gatherHalo->NrowsT = NhaloT;

  //tally up how many nodes are being gathered to each gatherNode and
  //  map to a local ordering
  memory<dlong> haloGatherNCounts(gatherHalo->NrowsT,0);
  memory<dlong> haloGatherTCounts(gatherHalo->NrowsT,0);

  for (dlong i=0;i<Nids;i++) {
    const dlong gid = nodes[i].newId; //re-mapped baseId on this rank

    if (abs(nodes[i].sign)==2) {//halo
      if (nodes[i].sign==2) haloGatherNCounts[gid]++;  //tally
      haloGatherTCounts[gid]++;  //tally
    }
  }

  //make halo row offsets
  gatherHalo->rowStartsN.malloc(gatherHalo->NrowsT+1);
  gatherHalo->rowStartsT.malloc(gatherHalo->NrowsT+1);
  gatherHalo->rowStartsN[0]=0;
  gatherHalo->rowStartsT[0]=0;
  for (dlong i=0;i<gatherHalo->NrowsT;i++) {
    gatherHalo->rowStartsN[i+1] = gatherHalo->rowStartsN[i] + haloGatherNCounts[i];
    gatherHalo->rowStartsT[i+1] = gatherHalo->rowStartsT[i] + haloGatherTCounts[i];
    haloGatherNCounts[i] = 0;
    haloGatherTCounts[i] = 0;
  }
  gatherHalo->nnzN = gatherHalo->rowStartsN[gatherHalo->NrowsT];
  gatherHalo->nnzT = gatherHalo->rowStartsT[gatherHalo->NrowsT];
  gatherHalo->colIdsN.malloc(gatherHalo->nnzN);
  gatherHalo->colIdsT.malloc(gatherHalo->nnzT);


  for (dlong i=0;i<Nids;i++) {
    const dlong gid = nodes[i].newId;

    if (abs(nodes[i].sign)==2) {
      if (nodes[i].sign==2) {
        const dlong soffset = gatherHalo->rowStartsN[gid];
        const int sindex  = haloGatherNCounts[gid];
        gatherHalo->colIdsN[soffset+sindex] = nodes[i].localId;
        haloGatherNCounts[gid]++;
      }
      const dlong soffset = gatherHalo->rowStartsT[gid];
      const int sindex  = haloGatherTCounts[gid];
      gatherHalo->colIdsT[soffset+sindex] = nodes[i].localId;
      haloGatherTCounts[gid]++;
    }
  }
  haloGatherNCounts.free();
  haloGatherTCounts.free();

  gatherHalo->o_rowStartsN = platform.malloc(gatherHalo->rowStartsN);
  gatherHalo->o_rowStartsT = platform.malloc(gatherHalo->rowStartsT);
  gatherHalo->o_colIdsN = platform.malloc(gatherHalo->colIdsN);
  gatherHalo->o_colIdsT = platform.malloc(gatherHalo->colIdsT);

  //divide the list of colIds into roughly equal sized blocks so that each
  // threadblock loads approximately an equal amount of data
  gatherHalo->setupRowBlocks();
}

void ogsBase_t::Free() {
  comm.Free();
  gatherLocal = nullptr;
  gatherHalo = nullptr;
  exchange = nullptr;
  N=0;
  NlocalT=0;
  NhaloT=0;
  Ngather=0;
  NgatherGlobal=0;
}

void ogsBase_t::AssertGatherDefined() {
  LIBP_ABORT("Gather operation not well-defined.",
             !gather_defined);
}

//Populate the local mapping of the original ids and the gathered ordering
void ogs_t::SetupGlobalToLocalMapping(memory<dlong> GlobalToLocal) {

  LIBP_ABORT("ogs handle is not set up.",
             NgatherGlobal==0);

  //Note: Must have GlobalToLocal have N entries.

  memory<dlong> ids(NlocalT+NhaloT);

  for (dlong n=0;n<NlocalT+NhaloT;n++)
    ids[n] = n;

  for (dlong n=0;n<N;n++)
    GlobalToLocal[n] = -1;

  gatherLocal->Scatter(GlobalToLocal, ids,
                       1, NoTrans);
  gatherHalo->Scatter(GlobalToLocal, ids+NlocalT,
                       1, NoTrans);
}

void halo_t::SetupFromGather(ogs_t& ogs) {

  ogs.AssertGatherDefined();

  platform = ogs.platform;
  comm = ogs.comm;

  N = ogs.NlocalT + ogs.NhaloT;

  Ngather = Ngather;
  Nhalo = ogs.NhaloT - ogs.NhaloP;

  NgatherGlobal = ogs.NgatherGlobal;

  kind = Halo;
  unique = ogs.unique;

  NlocalP = ogs.NlocalP;
  NlocalT  = ogs.NlocalT;

  NhaloP = ogs.NhaloP;
  NhaloT  = ogs.NhaloT;

  gather_defined=false;

  gathered_halo=true;

  exchange = ogs.exchange;
}

} //namespace ogs

} //namespace libp
