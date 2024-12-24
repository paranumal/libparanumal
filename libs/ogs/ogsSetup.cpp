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
#include "primitives.hpp"

namespace libp {

namespace ogs {


constexpr int FLAG_LOCAL=1;
constexpr int FLAG_SHARED=2;

static void EnumerateGatherGroups(const dlong N,
                                  const memory<hlong> baseIds,
                                  dlong& NgroupsP,
                                  dlong& NgroupsT,
                                  memory<dlong> gids);
static void SplitGroups(const dlong Nids,
                        const memory<dlong> colIds,
                        const memory<hlong> baseIds,
                        const memory<int> sharedFlag,
                        const dlong Nlocal,
                        const dlong Nhalo,
                        memory<dlong> localColIds,
                        memory<dlong> haloColIds,
                        memory<hlong> localBaseIds,
                        memory<hlong> haloBaseIds);
static void SortByRow(const dlong N,
                      memory<dlong>& rowIds,
                      memory<dlong>& colIds,
                      memory<hlong>& baseIds);

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

  // Seed RNG
  int rank = comm.rank();
  prim::seedRNG(rank);

  //sanity check options
  LIBP_ABORT("Invalid ogs setup requested",
             (kind==Unsigned && unique==true)
              || (kind==Halo && unique==true));


  memory<int> flags(N);

  //flag non-zero ids
  #pragma omp parallel for
  for (dlong n=0;n<N;++n) {
    flags[n] = (ids[n]!=0) ? 1 : 0;
  }

  //count how many ids are non-zero
  dlong Nids = prim::count(N, flags, 1);

  /*Make an array of locations for non-zero base ids*/
  memory<dlong> colIds(Nids);
  prim::select(N, flags, 1, colIds);
  flags.free();

  /*Compress out zero ids*/
  memory<hlong> baseIds(Nids);
  prim::transformGather(Nids, colIds, ids, baseIds);

  //flag which baseId groups are shared via MPI
  memory<int> sharedFlag(Nids);
  FindSharedGroups(Nids,
                  baseIds,
                  sharedFlag,
                  verbose);

  //if we altered the signs of ids, write them back
  if (unique) {
    prim::transformScatter(Nids, colIds, baseIds, ids);
  }

  /*
  Use the shared flag to split the lists of colIds and baseIds
  into local and halo lists
  */
  dlong Nlocal = prim::count(Nids, sharedFlag, FLAG_LOCAL);
  dlong Nhalo = Nids - Nlocal;
  memory<dlong> localColIds(Nlocal);
  memory<dlong> haloColIds(Nhalo);
  memory<hlong> localBaseIds(Nlocal);
  memory<hlong> haloBaseIds(Nhalo);

  SplitGroups(Nids,
              colIds,
              baseIds,
              sharedFlag,
              Nlocal,
              Nhalo,
              localColIds,
              haloColIds,
              localBaseIds,
              haloBaseIds);

  /*Number the gather groups to make the row indexing*/
  memory<dlong> localRowIds(Nlocal);
  memory<dlong> haloRowIds(Nhalo);

  EnumerateGatherGroups(Nlocal, localBaseIds, NlocalP, NlocalT, localRowIds);
  EnumerateGatherGroups(Nhalo, haloBaseIds, NhaloP, NhaloT, haloRowIds);

  /*Sort the sparse COO data by the row index*/
  SortByRow(Nlocal, localRowIds, localColIds, localBaseIds);
  SortByRow(Nhalo, haloRowIds, haloColIds, haloBaseIds);

  /*Make the local and halo gather operators*/
  ConstructGatherOperators(Nlocal, localRowIds, localColIds, localBaseIds,
                           Nhalo, haloRowIds, haloColIds, haloBaseIds);

  //total number of positive owned gathered nodes
  Ngather = NlocalP+NhaloP;

  //global total
  NgatherGlobal = Ngather;
  comm.Allreduce(NgatherGlobal);

  // Compress the list of baseIds in the halo into the list they appear after gathering
  memory<hlong> gatheredHaloBaseIds(NhaloT);
  #pragma omp parallel for
  for (dlong n=0;n<NhaloT;++n) {
    const hlong baseId = std::abs(haloBaseIds[gatherHalo->rowStartsT[n]]);
    gatheredHaloBaseIds[n] = (n<NhaloP) ? baseId : -baseId;
  }
  haloBaseIds = gatheredHaloBaseIds;

  //Construct sharedNodes which contains all the info
  // we need to setup the MPI exchange.
  dlong Nshared = 0;
  memory<int>   sharedRemoteRanks;
  memory<dlong> sharedLocalRows;
  memory<dlong> sharedRemoteRows;
  memory<hlong> sharedLocalBaseIds;
  memory<hlong> sharedRemoteBaseIds;
  ConstructSharedNodes(haloBaseIds,
                       Nshared,
                       sharedRemoteRanks,
                       sharedLocalRows,
                       sharedRemoteRows,
                       sharedLocalBaseIds,
                       sharedRemoteBaseIds);

  // At this point, we've setup gs operators to gather/scatter the purely local nodes,
  // and gather/scatter the shared halo nodes to/from a coalesced ordering. We now
  // need gs operators to scatter/gather the coalesced halo nodes to/from the expected
  // orderings for MPI communications.

  Kind knd = (kind == Unsigned) ? Unsigned : Signed;

  if (method == AllToAll) {
    exchange = std::shared_ptr<ogsExchange_t>(
                  new ogsAllToAll_t(knd,
                                    Nshared,
                                    sharedRemoteRanks,
                                    sharedLocalRows,
                                    sharedRemoteRows,
                                    sharedLocalBaseIds,
                                    sharedRemoteBaseIds,
                                    *gatherHalo,
                                    dataStream,
                                    comm,
                                    platform));
  } else if (method == Pairwise) {
    exchange = std::shared_ptr<ogsExchange_t>(
                  new ogsPairwise_t(knd,
                                    Nshared,
                                    sharedRemoteRanks,
                                    sharedLocalRows,
                                    sharedRemoteRows,
                                    sharedLocalBaseIds,
                                    sharedRemoteBaseIds,
                                    *gatherHalo,
                                    dataStream,
                                    comm,
                                    platform));
  } else if (method == CrystalRouter) {
    exchange = std::shared_ptr<ogsExchange_t>(
                  new ogsCrystalRouter_t(knd,
                                         Nshared,
                                         sharedRemoteRanks,
                                         sharedLocalRows,
                                         sharedRemoteRows,
                                         sharedLocalBaseIds,
                                         sharedRemoteBaseIds,
                                         haloBaseIds,
                                         *gatherHalo,
                                         dataStream,
                                         comm,
                                         platform));
  } else { //Auto
    exchange = std::shared_ptr<ogsExchange_t>(
                  AutoSetup(Nshared,
                            sharedRemoteRanks,
                            sharedLocalRows,
                            sharedRemoteRows,
                            sharedLocalBaseIds,
                            sharedRemoteBaseIds,
                            haloBaseIds,
                            *gatherHalo,
                            comm,
                            platform,
                            verbose));
  }

  timePoint_t end = GlobalPlatformTime(platform);
  double elapsedTime = ElapsedTime(start, end);

  if (!rank && verbose) {
    std::cout << "ogs Setup Time: " << elapsedTime << " seconds." << std::endl;
  }
}

void ogsBase_t::FindSharedGroups(const dlong Nids,
                                 memory<hlong> baseIds,
                                 memory<int> sharedFlag,
                                 const int verbose){

  int rank = comm.rank();
  int size = comm.size();

  /*Create list of destination ranks*/
  memory<int> destRanks(Nids);
  memory<dlong> destRankSortIds(Nids);

  #pragma omp parallel for
  for (dlong n=0;n<Nids;n++) { destRanks[n] = std::abs(baseIds[n]) % size; }

  /*Sort by destination to get the send ordering*/
  prim::sort(Nids, destRanks, destRankSortIds);

  memory<int> sendCounts(size);
  memory<int> recvCounts(size);
  memory<int> sendOffsets(size+1);
  memory<int> recvOffsets(size+1);

  /*Get length and offsets of groups of destinations*/
  prim::runLengthEncodeConsecutive(Nids, destRanks, size, sendOffsets);
  prim::adjacentDifference(size, sendOffsets+1, sendCounts);

  comm.Alltoall(sendCounts, recvCounts);
  recvOffsets[0] = 0;
  prim::inclusiveScan(size, recvCounts, recvOffsets+1);
  dlong recvN = recvOffsets[size]; //total ids to recv

  destRanks.free();

  //Arrange the list of baseIds to send order
  memory<hlong> sendBaseIds(Nids);
  prim::transformGather(Nids, destRankSortIds, baseIds, sendBaseIds);

  //Send the baseIds to their destination rank.
  memory<hlong> recvBaseIds(recvN);
  comm.Alltoallv(sendBaseIds, sendCounts, sendOffsets,
                 recvBaseIds, recvCounts, recvOffsets);

  // Send the source rank of each node to destination
  memory<hlong> sendRanks(Nids, rank); //initialize entries to 'rank'
  memory<hlong> recvRanks(recvN);
  comm.Alltoallv(sendRanks, sendCounts, sendOffsets,
                 recvRanks, recvCounts, recvOffsets);

  /*Group recieved nodes together by the abs(baseId)*/
  memory<hlong> recvAbsBaseIds(recvN);
  memory<dlong> baseSortIds(recvN);

  /*Sort by abs(baseId) and record ordering*/
  prim::abs(recvN, recvBaseIds, recvAbsBaseIds);
  prim::sort(recvN, recvAbsBaseIds, baseSortIds);

  /*Get offsets to groups of baseIds*/
  dlong NbaseIdGroups = 0;
  memory<dlong> baseIdGroupOffsets;
  prim::runLengthEncode(recvN, recvAbsBaseIds, NbaseIdGroups, baseIdGroupOffsets);

  recvAbsBaseIds.free();

  // We now have a collection of nodes associated with some subset of all global Ids
  // Our list is sorted by baseId to group nodes with the same globalId together

  int is_unique=1;

  if (unique) {
    //Make a single node from each baseId group the sole positive node

    memory<int> rands(NbaseIdGroups);
    prim::random(NbaseIdGroups, rands);

    #pragma omp parallel for
    for (dlong n=0;n<NbaseIdGroups;n++) {
      const dlong start = baseIdGroupOffsets[n];
      const dlong end = baseIdGroupOffsets[n+1];

      const hlong baseId = std::abs(recvBaseIds[baseSortIds[start]]);

      //pick a random node in this group
      const int m = (rands[n] % (end-start));

      for (dlong i=start;i<end;i++) {
        recvBaseIds[baseSortIds[i]] = -baseId;
      }

      recvBaseIds[baseSortIds[start+m]] = baseId;
    }

  } else {

    memory<dlong> uniqueFlag(NbaseIdGroups, 0);

    #pragma omp parallel for
    for (dlong n=0;n<NbaseIdGroups;n++) {
      const dlong start = baseIdGroupOffsets[n];
      const dlong end = baseIdGroupOffsets[n+1];

      int positiveCount=0;
      //count how many postive baseIds there are in this group
      for (dlong i=start;i<end;i++) {
        positiveCount += (recvBaseIds[baseSortIds[i]] > 0) ? 1 : 0;
      }

      //if we didnt find a sole positive baseId, flag this group. The gather is not well-defined
      uniqueFlag[n] = (positiveCount==1) ? 1 : 0;
    }

    is_unique = (prim::count(NbaseIdGroups, uniqueFlag, 0)==0) ? 1 : 0;

    if (kind==Halo && !is_unique) {
      // When making a halo exchange, we have to have a single positive id in each baseId group
      LIBP_FORCE_ABORT("Halo exchange not well-defined. Some baseId groups have no unique positive id.");
    }
  }

  //shared the unique node check so we know if the gather operation is well-defined
  comm.Allreduce(is_unique, comm_t::Min);
  gather_defined = (is_unique==1);


  memory<int> recvSharedFlag(recvN);
  memory<int> sendSharedFlag(Nids);
  memory<int> flags(NbaseIdGroups); //flag for whether baseId group is shared over MPI

  // We now want to flag which nodes are shared via MPI
  #pragma omp parallel for
  for (dlong n=0;n<NbaseIdGroups;n++) {
    const dlong start = baseIdGroupOffsets[n];
    const dlong end = baseIdGroupOffsets[n+1];

    //determine if this node is shared via MPI,
    int flag=FLAG_LOCAL;
    const int r = recvRanks[baseSortIds[start]];
    for (dlong i=start+1;i<end;i++) {
      if (recvRanks[baseSortIds[i]] != r) {
        flag=FLAG_SHARED;
        break;
      }
    }

    //set shared flag.
    for (dlong i=start;i<end;i++) {
      recvSharedFlag[baseSortIds[i]] = flag;
    }

    flags[n] = flag;
  }

  dlong Nshared = prim::count(NbaseIdGroups, flags, FLAG_SHARED);
  hlong Nshared_global = Nshared;
  comm.Reduce(Nshared_global, 0);

  //at this point each collection of baseIds either has all nodes have
  // sign = 1, meaning all the nodes with this baseId are on the
  // same rank, or have sign=2, meaning that baseId must be communicated

  //Share the group signs back to their source
  comm.Alltoallv(recvSharedFlag, recvCounts, recvOffsets,
                 sendSharedFlag, sendCounts, sendOffsets);

  prim::transformScatter(Nids, destRankSortIds, sendSharedFlag, sharedFlag);

  if (unique) {
    //Share the group signs back to their source
    comm.Alltoallv(recvBaseIds, recvCounts, recvOffsets,
                   sendBaseIds, sendCounts, sendOffsets);

    prim::transformScatter(Nids, destRankSortIds, sendBaseIds, baseIds);
  }

  if (!rank && verbose) {
    std::cout << "ogs Setup: " << Nshared_global << " unique labels shared." << std::endl;
  }
}

static void SplitGroups(const dlong Nids,
                        const memory<dlong> colIds,
                        const memory<hlong> baseIds,
                        const memory<int> sharedFlag,
                        const dlong Nlocal,
                        const dlong Nhalo,
                        memory<dlong> localColIds,
                        memory<dlong> haloColIds,
                        memory<hlong> localBaseIds,
                        memory<hlong> haloBaseIds) {
  memory<dlong> localIds(Nlocal);
  memory<dlong> haloIds(Nhalo);
  prim::select(Nids, sharedFlag, FLAG_LOCAL, localIds);
  prim::select(Nids, sharedFlag, FLAG_SHARED, haloIds);

  prim::transformGather(Nlocal, localIds, baseIds, localBaseIds);
  prim::transformGather(Nhalo,   haloIds, baseIds, haloBaseIds);
  prim::transformGather(Nlocal, localIds, colIds, localColIds);
  prim::transformGather(Nhalo,   haloIds, colIds, haloColIds);
}

static void EnumerateGatherGroups(const dlong N,
                                  const memory<hlong> baseIds,
                                  dlong& NgroupsP,
                                  dlong& NgroupsT,
                                  memory<dlong> gids) {

  /*
  Enumerate each baseId groups according to a gathered ordering,
  groups with at least 1 positive baseId appear first in the
  gathered ordering
  */
  memory<hlong> absBaseIds(N);
  memory<dlong> sortIds(N);

  /*Sort by abs(baseId) to group baseIds*/
  prim::abs(N, baseIds, absBaseIds);
  prim::stableSort(N, absBaseIds, sortIds);

  /*Compute how many baseId groups we have, and get offsets to groups of baseIds*/
  NgroupsT = 0;
  memory<dlong> groupOffsets;
  prim::runLengthEncode(N, absBaseIds, NgroupsT, groupOffsets);

  prim::set(N, 0, gids);

  /* Sign groups that have no positive baseId on this rank and mark the first appearance */
  #pragma omp parallel for
  for (dlong n=0;n<NgroupsT;++n) {
    const dlong start = groupOffsets[n];
    const dlong end = groupOffsets[n+1];

    // Check for a positive baseId in this group
    int sign = -1;
    for (dlong i=start;i<end;++i) {
      if (baseIds[sortIds[i]] > 0) {
        sign = 1; //Found a positive baseId in this group
        break;
      }
    }
    /*Since the sort was stable, 'start' is the first appearance of this baseId group*/
    gids[sortIds[start]] = sign;
  }

  NgroupsP = prim::count(N, gids, 1);

  /* Get the ids of the first appearances of each baseId group, in their original ordering */
  memory<dlong> gatherIds(NgroupsT);
  prim::select(N, gids,  1, gatherIds);
  prim::select(N, gids, -1, gatherIds+NgroupsP);

  /*Enumerate the first entry of the groups*/
  #pragma omp parallel for
  for (dlong n=0;n<NgroupsT;n++) {
    gids[gatherIds[n]] = n;
  }

  /* Propagate numbering to whole group */
  #pragma omp parallel for
  for (dlong n=0;n<NgroupsT;++n) {
    const dlong start = groupOffsets[n];
    const dlong end = groupOffsets[n+1];

    const dlong gid = gids[sortIds[start]];
    for (dlong i=start+1;i<end;i++) {
      gids[sortIds[i]] = gid;
    }
  }
}

static void SortByRow(const dlong N,
                      memory<dlong>& rowIds,
                      memory<dlong>& colIds,
                      memory<hlong>& baseIds) {

  memory<dlong> sortIds(N);

  /*Sort groups by their row*/
  prim::stableSort(N, rowIds, sortIds);

  memory<hlong> newBaseIds(N);
  prim::transformGather(N, sortIds, baseIds, newBaseIds);

  memory<dlong> newColIds(N);
  prim::transformGather(N, sortIds, colIds, newColIds);

  baseIds = newBaseIds;
  colIds = newColIds;
}

void ogsBase_t::ConstructGatherOperators(const dlong Nlocal,
                                         const memory<dlong> localRowIds,
                                         const memory<dlong> localColIds,
                                         const memory<hlong> localBaseIds,
                                         const dlong Nhalo,
                                         const memory<dlong> haloRowIds,
                                         const memory<dlong> haloColIds,
                                         const memory<hlong> haloBaseIds) {

  Kind knd = (kind == Unsigned) ? Unsigned : Signed;

  //setup local gather operator
  if (kind!=Halo) {
    gatherLocal = std::make_shared<ogsOperator_t>(platform,
                                                  knd,
                                                  NlocalP,
                                                  NlocalT,
                                                  N,
                                                  Nlocal,
                                                  localBaseIds,
                                                  localRowIds,
                                                  localColIds);
  }

  //setup local gather operator
  gatherHalo = std::make_shared<ogsOperator_t>(platform,
                                               knd,
                                               NhaloP,
                                               NhaloT,
                                               N,
                                               Nhalo,
                                               haloBaseIds,
                                               haloRowIds,
                                               haloColIds);
}

void ogsBase_t::ConstructSharedNodes(const memory<hlong> haloBaseIds,
                                     dlong &Nshared,
                                     memory<int>&   sharedRemoteRanks,
                                     memory<dlong>& sharedLocalRows,
                                     memory<dlong>& sharedRemoteRows,
                                     memory<hlong>& sharedLocalBaseIds,
                                     memory<hlong>& sharedRemoteBaseIds) {

  int rank = comm.rank();
  int size = comm.size();

  memory<dlong> rows(NhaloT);

  /*Create list of destination ranks*/
  memory<int> destRanks(NhaloT);
  memory<dlong> destRankSortIds(NhaloT);

  prim::range(NhaloT, 0, 1, rows);

  #pragma omp parallel for
  for (dlong n=0;n<NhaloT;++n) { destRanks[n] = std::abs(haloBaseIds[n]) % size; }

  /*Sort by destination to get the send ordering*/
  prim::sort(NhaloT, destRanks, destRankSortIds);

  memory<int> sendCounts(size);
  memory<int> recvCounts(size);
  memory<int> sendOffsets(size+1);
  memory<int> recvOffsets(size+1);

  /*Get length and offsets of groups of destinations*/
  prim::runLengthEncodeConsecutive(NhaloT, destRanks, size, sendOffsets);
  prim::adjacentDifference(size, sendOffsets+1, sendCounts);

  comm.Alltoall(sendCounts, recvCounts);
  recvOffsets[0] = 0;
  prim::inclusiveScan(size, recvCounts, recvOffsets+1);
  dlong recvN = recvOffsets[size]; //total ids to recv

  destRanks.free();

  //Arrange the list of haloBaseIds to send order
  memory<hlong> sendBaseIds(NhaloT);
  prim::transformGather(NhaloT, destRankSortIds, haloBaseIds, sendBaseIds);

  //Send the haloBaseIds to their destination rank.
  memory<hlong> recvBaseIds(recvN);
  comm.Alltoallv(sendBaseIds, sendCounts, sendOffsets,
                 recvBaseIds, recvCounts, recvOffsets);

  //Arrange the list of rows to send order
  memory<dlong> sendRows(NhaloT);
  prim::transformGather(NhaloT, destRankSortIds, rows, sendRows);

  //Send the rows to their destination rank.
  memory<dlong> recvRows(recvN);
  comm.Alltoallv(sendRows, sendCounts, sendOffsets,
                 recvRows, recvCounts, recvOffsets);

  // Send the source rank of each node to destination
  memory<int> sendRanks(NhaloT, rank); //initialize entries to 'rank'
  memory<int> recvRanks(recvN);
  comm.Alltoallv(sendRanks, sendCounts, sendOffsets,
                 recvRanks, recvCounts, recvOffsets);

  //free up some space
  destRankSortIds.free();
  rows.free();
  sendBaseIds.free();
  sendRows.free();
  sendRanks.free();

  /*Group recieved nodes together by the abs(baseId)*/
  memory<hlong> recvAbsBaseIds(recvN);
  memory<dlong> baseSortIds(recvN);

  /*Sort by abs(baseId) and record ordering*/
  prim::abs(recvN, recvBaseIds, recvAbsBaseIds);
  prim::sort(recvN, recvAbsBaseIds, baseSortIds);

  /*Get offsets to groups of baseIds*/
  dlong NbaseIdGroups = 0;
  memory<dlong> baseIdGroupOffsets;
  prim::runLengthEncode(recvN, recvAbsBaseIds, NbaseIdGroups, baseIdGroupOffsets);

  recvAbsBaseIds.free();

  // We now have a collection of nodes associated with some subset of all global Ids
  // Each of these baseId groups are shared between at least 2 MPI ranks.
  // We want to send the full sharing information back to every rank that is participating
  // in the group so that, for a given baseId group, each rank participating knows:
  // 1) All ranks sharing this baseId 2) the +/- sign of the baseId on each rank
  // 3) The gathered index of the baseId on each rank

  /*
  Each node in the list is going to record info for the whole baseId group,
  so grow the array and get offsets to where each node writes
  */
  memory<dlong> groupSize(recvN);
  memory<dlong> groupStart(recvN);
  memory<dlong> groupEnd(recvN);

  #pragma omp parallel for
  for (dlong n=0;n<NbaseIdGroups;++n) {
    const dlong start = baseIdGroupOffsets[n];
    const dlong end = baseIdGroupOffsets[n+1];

    //write the group size at every node location (excluding the rank's own node)
    for (dlong i=start;i<end;++i) {
      groupStart[baseSortIds[i]] = start;
      groupEnd[baseSortIds[i]] = end;
      groupSize[baseSortIds[i]] = end-start-1;
    }
  }
  memory<dlong> groupOffsets(recvN+1);
  groupOffsets[0] = 0;
  prim::inclusiveScan(recvN, groupSize, groupOffsets+1);

  /*Now fill the send data for each node*/
  dlong sendN = groupOffsets[recvN];
  sendRanks.malloc(sendN);
  sendBaseIds.malloc(sendN);
  memory<dlong> sendLocalRows(sendN);
  memory<dlong> sendRemoteRows(sendN);

  #pragma omp parallel for
  for (dlong n=0;n<recvN;++n) {
    const dlong start = groupOffsets[n];
    const dlong gStart = groupStart[n];
    const dlong gEnd = groupEnd[n];

    const dlong row = recvRows[n];

    //write all the group info
    dlong loc = start;
    for (dlong i=gStart;i<gEnd;++i) {
      const dlong gloc = baseSortIds[i];
      if (gloc==n) continue; //dont bother sending this rank's own node
      sendRanks[loc] = recvRanks[gloc];
      sendBaseIds[loc] = recvBaseIds[gloc];
      sendRemoteRows[loc] = recvRows[gloc];
      sendLocalRows[loc] = row;
      ++loc;
    }
  }
  baseSortIds.free();
  recvRows.free();
  recvRanks.free();
  recvBaseIds.free();

  /*Update the MPI offset data for the new payload size*/
  prim::transformGather(size+1, recvOffsets, groupOffsets, sendOffsets);
  prim::adjacentDifference(size, sendOffsets+1, sendCounts);

  comm.Alltoall(sendCounts, recvCounts);
  recvOffsets[0] = 0;
  prim::inclusiveScan(size, recvCounts, recvOffsets+1);
  Nshared = recvOffsets[size]; //total ids to recv

  recvRanks.malloc(Nshared);
  recvBaseIds.malloc(Nshared);
  memory<dlong> recvLocalRows(Nshared);
  memory<dlong> recvRemoteRows(Nshared);

  comm.Alltoallv(sendRanks, sendCounts, sendOffsets,
                 recvRanks, recvCounts, recvOffsets);
  comm.Alltoallv(sendBaseIds, sendCounts, sendOffsets,
                 recvBaseIds, recvCounts, recvOffsets);
  comm.Alltoallv(sendLocalRows, sendCounts, sendOffsets,
                 recvLocalRows, recvCounts, recvOffsets);
  comm.Alltoallv(sendRemoteRows, sendCounts, sendOffsets,
                 recvRemoteRows, recvCounts, recvOffsets);

  memory<dlong> localSortIds(Nshared);

  /*Sort list of shared nodes by their local row ordering*/
  prim::sort(Nshared, recvLocalRows, localSortIds);

  sharedLocalRows = recvLocalRows;

  sharedRemoteRanks.malloc(Nshared);
  prim::transformGather(Nshared, localSortIds, recvRanks, sharedRemoteRanks);

  sharedRemoteBaseIds.malloc(Nshared);
  prim::transformGather(Nshared, localSortIds, recvBaseIds, sharedRemoteBaseIds);

  sharedRemoteRows.malloc(Nshared);
  prim::transformGather(Nshared, localSortIds, recvRemoteRows, sharedRemoteRows);

  sharedLocalBaseIds.malloc(Nshared);

  #pragma omp parallel for
  for (dlong n=0;n<Nshared;++n) {
    const hlong baseId = std::abs(sharedRemoteBaseIds[n]);
    sharedLocalBaseIds[n] = (sharedLocalRows[n]<NhaloP) ? baseId : -baseId;
  }
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

  prim::range(NlocalT+NhaloT, 0, 1, ids);
  prim::set(N, -1, GlobalToLocal);

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
