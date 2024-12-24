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
#include "ogs/ogsExchange.hpp"
#include "primitives.hpp"

namespace libp {

namespace ogs {

/**********************************
* Host exchange
***********************************/
template<typename T>
inline void ogsCrystalRouter_t::HostStart(const int k, const Op op, const Transpose trans){}

template<typename T>
inline void ogsCrystalRouter_t::HostFinish(const int k, const Op op, const Transpose trans){

  memory<crLevel>& L = data[trans].levels;
  pinnedMemory<T> workBuf = h_workspace;

  for (int l=0;l<Nlevels;l++) {

    pinnedMemory<T> sendBuf = h_sendspace;
    pinnedMemory<T> recvBuf = h_recvspace;

    //post recvs
    if (L[l].Nmsg>0) {
      comm.Irecv(sendBuf + L[l].Nids*k,
                 L[l].partner,
                 k*L[l].Nrecv0,
                 L[l].partner,
                 requests[1]);
    }
    if (L[l].Nmsg==2) {
      comm.Irecv(sendBuf + L[l].Nids*k + L[l].Nrecv0*k,
                rank-1,
                k*L[l].Nrecv1,
                rank-1,
                requests[2]);
    }

    //assemble send buffer
    extract(L[l].Nsend, k, L[l].sendIds, sendBuf, workBuf);

    //post send
    comm.Isend(workBuf,
               L[l].partner,
               k*L[l].Nsend,
               rank,
               requests[0]);

    comm.Waitall(L[l].Nmsg+1, requests);

    //Gather the recv'd values into the haloBuffer
    L[l].gather.Gather(recvBuf, sendBuf, k, op, Trans);

    if (l==Nlevels-1) break;

    //swap buffers
    h_sendspace.swap(h_recvspace);
  }
}

void ogsCrystalRouter_t::HostStart(const Type type, const int k, const Op op, const Transpose trans) {
  switch (type) {
    case Int32:  HostStart<int>(k, op, trans); break;
    case Int64:  HostStart<long long int>(k, op, trans); break;
    case Float:  HostStart<float>(k, op, trans); break;
    case Double: HostStart<double>(k, op, trans); break;
  }
}
void ogsCrystalRouter_t::HostFinish(const Type type, const int k, const Op op, const Transpose trans) {
  switch (type) {
    case Int32:  HostFinish<int>(k, op, trans); break;
    case Int64:  HostFinish<long long int>(k, op, trans); break;
    case Float:  HostFinish<float>(k, op, trans); break;
    case Double: HostFinish<double>(k, op, trans); break;
  }
}

/**********************************
* GPU-aware exchange
***********************************/
template<typename T>
inline void ogsCrystalRouter_t::DeviceStart(const int k, const Op op, const Transpose trans){
  //wait for kernel to finish on default stream
  device_t &device = platform.device;
  device.finish();
}

template<typename T>
inline void ogsCrystalRouter_t::DeviceFinish(const int k, const Op op, const Transpose trans){

  device_t &device = platform.device;

  //get current stream
  stream_t currentStream = device.getStream();

  //the intermediate kernels are always overlapped with the default stream
  device.setStream(dataStream);

  memory<crLevel>& L = data[trans].levels;
  deviceMemory<T> o_workBuf = o_workspace;

  for (int l=0;l<Nlevels;l++) {
    deviceMemory<T> o_sendBuf = o_sendspace;
    deviceMemory<T> o_recvBuf = o_recvspace;

    //post recvs
    if (L[l].Nmsg>0) {
      comm.Irecv(o_sendBuf + L[l].Nids*k,
                 L[l].partner,
                 k*L[l].Nrecv0,
                 L[l].partner,
                 requests[1]);
    }
    if (L[l].Nmsg==2) {
      comm.Irecv(o_sendBuf + L[l].Nids*k + L[l].Nrecv0*k,
                rank-1,
                k*L[l].Nrecv1,
                rank-1,
                requests[2]);
    }

    //assemble send buffer
    if (L[l].Nsend) {
      extractKernel[ogsType<T>::get()](L[l].Nsend, k,
                                       L[l].o_sendIds,
                                       o_sendBuf, o_workBuf);
      device.finish();
    }

    //post send
    comm.Isend(o_workBuf,
               L[l].partner,
               k*L[l].Nsend,
               rank,
               requests[0]);

    comm.Waitall(L[l].Nmsg+1, requests);

    //Gather the recv'd values into the haloBuffer
    L[l].gather.Gather(o_recvBuf, o_sendBuf, k, op, Trans);

    if (l==Nlevels-1) break;

    //swap buffers
    o_sendspace.swap(o_recvspace);
  }

  device.finish();
  device.setStream(currentStream);
}

void ogsCrystalRouter_t::DeviceStart(const Type type, const int k, const Op op, const Transpose trans) {
  switch (type) {
    case Int32:  DeviceStart<int>(k, op, trans); break;
    case Int64:  DeviceStart<long long int>(k, op, trans); break;
    case Float:  DeviceStart<float>(k, op, trans); break;
    case Double: DeviceStart<double>(k, op, trans); break;
  }
}
void ogsCrystalRouter_t::DeviceFinish(const Type type, const int k, const Op op, const Transpose trans) {
  switch (type) {
    case Int32:  DeviceFinish<int>(k, op, trans); break;
    case Int64:  DeviceFinish<long long int>(k, op, trans); break;
    case Float:  DeviceFinish<float>(k, op, trans); break;
    case Double: DeviceFinish<double>(k, op, trans); break;
  }
}

/*
 *Crystal Router performs the needed MPI communcation via recursive
 * folding of a hypercube. Consider a set of NP ranks. We select a
 * pivot point n_half=(NP+1)/2, and pair all ranks r<n_half (called
 * lo half) the with ranks r>=n_half (called the hi half), as follows
 *
 *                0 <--> NP-1
 *                1 <--> NP-2
 *                2 <--> NP-3
 *                  * * *
 *         n_half-2 <--> NP-n_half+1
 *         n_half-1 <--> NP-n_half
 *
 * The communication can then be summarized thusly: if a rank in the lo
 * half has data needed by *any* rank in the hi half, it sends this data
 * to its hi partner, and analogously for ranks in the hi half. Each rank
 * therefore sends/receives a single message to/from its partner.
 *
 * The communication then proceeds recursively, applying the same folding
 * proceedure to the lo and hi halves seperately, and stopping when the size
 * of the local NP reaches 1.
 *
 * In the case where NP is odd, n_half-1 == NP-n_half and rank n_half-1 has
 * no partner to communicate with. In this case, we assign rank r to the
 * lo half of ranks, and rank n_half-1 sends its data to rank n_half (and
 * receives no message, as rank n_half-2 is receiving all rank n_half's data).

 * To perform the Crystal Router exchange, each rank gathers its halo nodes to
 * a coalesced buffer. At each step in the crystal router, a send buffer is
 * gathered from this buffer and sent to the rank's partner. Simultaneously, a
 * buffer is received from the rank's partner. This receive buffer is scattered
 * and added into the coalesced halo buffer. After all commincation is complete
 * the halo nodes are scattered back to the output array.
 */

static void EnumerateGroups(const dlong N,
                            const memory<hlong> baseIds,
                            dlong& Ngroups,
                            memory<dlong> gids,
                            memory<hlong>& newBaseIds) {

  /*
  Enumerate each baseId groups according to a gathered ordering
  */
  memory<hlong> bIds(N);
  memory<dlong> sortIds(N);

  /*Copy the baseIds*/
  bIds.copyFrom(baseIds);

  /*Sort to group baseIds*/
  prim::stableSort(N, bIds, sortIds);

  /*Compute how many baseId groups we have, and get offsets to groups of baseIds*/
  Ngroups = 0;
  memory<dlong> groupOffsets;
  prim::runLengthEncode(N, bIds, Ngroups, groupOffsets);

  prim::set(N, 0, gids);

  /* Mark the first appearance of each baseId group*/
  #pragma omp parallel for
  for (dlong n=0;n<Ngroups;++n) {
    const dlong start = groupOffsets[n];
    /*Since the sort was stable, 'start' is the first appearance of this baseId group*/
    gids[sortIds[start]] = 1;
  }

  Ngroups = prim::count(N, gids, 1);

  /* Get the ids of the first appearances of each baseId group, in their original ordering */
  memory<dlong> gatherIds(Ngroups);
  prim::select(N, gids,  1, gatherIds);

  /*Enumerate the first entry of the groups*/
  #pragma omp parallel for
  for (dlong n=0;n<Ngroups;n++) {
    gids[gatherIds[n]] = n;
  }

  newBaseIds.malloc(Ngroups);

  /* Propagate numbering to whole group */
  #pragma omp parallel for
  for (dlong n=0;n<Ngroups;++n) {
    const dlong start = groupOffsets[n];
    const dlong end = groupOffsets[n+1];

    const dlong gid = gids[sortIds[start]];
    for (dlong i=start+1;i<end;i++) {
      gids[sortIds[i]] = gid;
    }
    newBaseIds[gid] = bIds[start]; //Record the baseId of this group
  }
}

void ogsCrystalRouter_t::data_t::setupExchange(const int Nlevels,
                                               dlong NhaloP,
                                               dlong NhaloT,
                                               memory<hlong> haloBaseIds,
                                               dlong Nshared,
                                               memory<int>   remoteRanks,
                                               memory<hlong> remoteBaseIds,
                                               memory<dlong> localIds,
                                               comm_t comm,
                                               platform_t& platform) {

  if (Nlevels<1) return;

  int rank = comm.rank();
  int size = comm.size();

  levels.malloc(Nlevels);

  dlong Nids = NhaloP;
  memory<hlong> baseIds(Nids);
  baseIds.copyFrom(haloBaseIds);

  dlong Ncols=0, nnz=0;
  memory<dlong> cols;

  // The list of shared nodes to send is already sorted by baseId, so get group offsets
  dlong Ngroups = 0;
  memory<dlong> groupOffsets;
  prim::runLengthEncode(Nshared, remoteBaseIds, Ngroups, groupOffsets);

  memory<dlong> haloIds(Ngroups);
  prim::transformGather(Ngroups, groupOffsets, localIds, haloIds);

  comm_t::request_t requests[3];

  //Now build the levels
  int lvl = 0;
  int np = size;
  int np_offset=0;

  NsendMax = 0;
  NrecvMax = 0;

  while (np>1) {
    int np_half = (np+1)/2;
    int r_half = np_half + np_offset;

    int is_lo = (rank<r_half) ? 1 : 0;

    int partner = np-1-(rank-np_offset)+np_offset;
    int Nmsg=1;
    if (partner==rank) {
      partner=r_half;
      Nmsg=0;
    }
    if (np&1 && rank==r_half) {
      Nmsg=2;
    }
    levels[lvl].partner = partner;
    levels[lvl].Nmsg = Nmsg;

    levels[lvl].Nids = Nids;

    //flag lo/hi nodes and groups
    memory<int> hiFlags(Nids); //Flag if this group contributes to the hi part
    memory<int> loFlags(Nids); //Flag if this group contributes to the lo part
    memory<int> nodeFlags(Nshared); //Flag if individual shared node is in hi part

    if (Nids != Ngroups) {
      // This can happen if some ids in the halo dont actually get sent anywhere,
      // For example, in the Trans mode when a baseId group has no positive id anywhere
      prim::set(Nids, -1, hiFlags);
      prim::set(Nids, -1, loFlags);
    }

    #pragma omp parallel for
    for (dlong n=0;n<Ngroups;++n) {
      const dlong start = groupOffsets[n];
      const dlong end = groupOffsets[n+1];

      const dlong lid = haloIds[n]; //location of this baseId in the halo

      int hiFlag = 0;
      int loFlag = 0;
      for (dlong i=start;i<end;++i) {
        const int flag = (remoteRanks[i]<r_half) ? 0 : 1;
        hiFlag |= flag;
        loFlag |= (flag) ? 0 : 1;
        nodeFlags[i] = flag;
      }
      hiFlags[lid] = hiFlag;
      loFlags[lid] = loFlag;
    }

    dlong NidsLo = prim::count(Nids, loFlags, 1);
    dlong NidsHi = prim::count(Nids, hiFlags, 1);

    // Make the list of sendIds
    int Nsend = (is_lo) ? NidsHi : NidsLo;
    int Nkeep = (is_lo) ? NidsLo : NidsHi;
    memory<dlong> sendFlags = (is_lo) ? hiFlags : loFlags;
    memory<dlong> keepFlags = (is_lo) ? loFlags : hiFlags;

    levels[lvl].Nsend = Nsend;
    levels[lvl].sendIds.malloc(Nsend);
    prim::select(Nids, sendFlags, 1, levels[lvl].sendIds);
    levels[lvl].o_sendIds = platform.malloc(levels[lvl].sendIds);

    // Communicate how many baseIds we're sending
    comm.Isend(Nsend, partner, rank, requests[0]);

    int Nrecv0=0, Nrecv1=0;
    if (Nmsg>0)
      comm.Irecv(Nrecv0, partner, partner, requests[1]);
    if (Nmsg==2)
      comm.Irecv(Nrecv1, r_half-1, r_half-1, requests[2]);

    comm.Waitall(Nmsg+1, requests);

    levels[lvl].Nrecv0 = Nrecv0;
    levels[lvl].Nrecv1 = Nrecv1;

    int Nrecv = Nrecv0+Nrecv1;

    //Size of recv buffer
    Ncols = Nids + Nrecv;

    NsendMax = std::max(NsendMax, Nsend);
    NrecvMax = std::max(NrecvMax, Ncols);

    nnz = Nkeep + Nrecv;
    cols.malloc(nnz);
    memory<hlong> newBaseIds(nnz);

    memory<hlong> sendBaseIds(Nsend);
    prim::transformGather(Nsend, levels[lvl].sendIds, baseIds, sendBaseIds);

    // Send the list of baseIds to our partner
    comm.Isend(sendBaseIds, partner, Nsend, rank, requests[0]);

    if (Nmsg>0)
      comm.Irecv(newBaseIds+Nkeep, partner, Nrecv0, partner, requests[1]);
    if (Nmsg==2)
      comm.Irecv(newBaseIds+Nkeep+Nrecv0, r_half-1, Nrecv1, r_half-1, requests[2]);

    prim::select(Nids, keepFlags, 1, cols);
    prim::transformGather(Nkeep, cols, baseIds, newBaseIds);

    //Column ids of the recieved baseIds
    prim::range(Nrecv, Nids, 1, cols+Nkeep);

    comm.Waitall(Nmsg+1, requests);

    sendBaseIds.free();
    keepFlags.free();
    sendFlags.free();
    hiFlags.free();
    loFlags.free();
    localIds.free();
    groupOffsets.free();

    //Shrink the size of the hypercube
    if (is_lo) {
      np = np_half;
    } else {
      np -= np_half;
      np_offset = r_half;
    }

    //The last gather must be built with the desired row ordering, so exit early here
    if (np<=1) {
      baseIds = newBaseIds;
      break;
    }

    //We now have the list of baseIds on this rank after the comms.
    // Build the gather to compress the list to a unique set of baseIds

    memory<dlong> rows(nnz);
    EnumerateGroups(nnz, newBaseIds, Nids, rows, baseIds);

    /*Sort groups by their row*/
    memory<dlong> rowSortIds(nnz);
    prim::stableSort(nnz, rows, rowSortIds);

    memory<dlong> sortedCols(nnz);
    prim::transformGather(nnz, rowSortIds, cols, sortedCols);
    cols = sortedCols;

    /*Build the gather op to assemble the recieved data from MPI*/
    levels[lvl].gather = ogsOperator_t(platform,
                                       Unsigned,
                                       Nids,
                                       Nids,
                                       Ncols,
                                       nnz,
                                       memory<hlong>(),
                                       rows,
                                       cols);

    // To construct the next level, we need to forward the sharing info to our partner

    dlong NsharedLo = prim::count(Nshared, nodeFlags, 0);
    dlong NsharedHi = Nshared - NsharedLo;

    Nsend = (is_lo) ? NsharedHi : NsharedLo;
    Nkeep = (is_lo) ? NsharedLo : NsharedHi;


    // Communicate how many baseIds we're sending
    comm.Isend(Nsend, partner, rank, requests[0]);

    Nrecv0=0;
    Nrecv1=0;
    if (Nmsg>0)
      comm.Irecv(Nrecv0, partner, partner, requests[1]);
    if (Nmsg==2)
      comm.Irecv(Nrecv1, r_half-1, r_half-1, requests[2]);

    comm.Waitall(Nmsg+1, requests);

    Nrecv = Nrecv0+Nrecv1;

    dlong NsharedNew = Nkeep + Nrecv;

    memory<dlong> sendIds(Nsend);
    memory<dlong> keepIds(Nkeep);
    prim::select(Nshared, nodeFlags, ((is_lo) ? 1 : 0), sendIds);
    prim::select(Nshared, nodeFlags, ((is_lo) ? 0 : 1), keepIds);

    memory<int> sendRemoteRanks(Nsend);
    prim::transformGather(Nsend, sendIds, remoteRanks, sendRemoteRanks);

    memory<int> newRemoteRanks(NsharedNew);

    // Send the list of ranks to our partner
    comm.Isend(sendRemoteRanks, partner, Nsend, rank, requests[0]);

    if (Nmsg>0)
      comm.Irecv(newRemoteRanks+Nkeep, partner, Nrecv0, partner, requests[1]);
    if (Nmsg==2)
      comm.Irecv(newRemoteRanks+Nkeep+Nrecv0, r_half-1, Nrecv1, r_half-1, requests[2]);

    prim::transformGather(Nkeep, keepIds, remoteRanks, newRemoteRanks);
    comm.Waitall(Nmsg+1, requests);

    memory<hlong> sendRemoteBaseIds(Nsend);
    prim::transformGather(Nsend, sendIds, remoteBaseIds, sendRemoteBaseIds);

    memory<hlong> newRemoteBaseIds(NsharedNew);

    // Send the list of baseIds to our partner
    comm.Isend(sendRemoteBaseIds, partner, Nsend, rank, requests[0]);

    if (Nmsg>0)
      comm.Irecv(newRemoteBaseIds+Nkeep, partner, Nrecv0, partner, requests[1]);
    if (Nmsg==2)
      comm.Irecv(newRemoteBaseIds+Nkeep+Nrecv0, r_half-1, Nrecv1, r_half-1, requests[2]);

    prim::transformGather(Nkeep, keepIds, remoteBaseIds, newRemoteBaseIds);
    comm.Waitall(Nmsg+1, requests);

    Nshared = NsharedNew;
    remoteBaseIds.free();
    remoteRanks.free();
    remoteBaseIds = newRemoteBaseIds;

    // Finally, sort the new list to group baseIds together

    memory<dlong> sortIds(Nshared);
    prim::stableSort(Nshared, newRemoteBaseIds, sortIds);

    remoteRanks.malloc(Nshared);
    prim::transformGather(Nshared, sortIds, newRemoteRanks, remoteRanks);
    newRemoteRanks.free();

    // Get group offsets
    prim::runLengthEncode(Nshared, remoteBaseIds, Ngroups, groupOffsets);

    // Sort the baseIds that comprise our current buffer to get a baseId->haloId mapping
    memory<hlong> bIds(Nids);
    bIds.copyFrom(baseIds);
    haloIds.malloc(Nids);

    prim::sort(Ngroups, bIds, haloIds);
    bIds.free();

    lvl++;
  }

  // Build the last gather to assemble the recieved buffer to the desired halo ordering

  memory<hlong> hBaseIds(NhaloT);
  hBaseIds.copyFrom(haloBaseIds);
  haloIds.malloc(NhaloT);

  prim::sort(NhaloT, hBaseIds, haloIds);

  memory<dlong> sortIds(nnz);
  prim::sort(nnz, baseIds, sortIds);
  prim::runLengthEncode(nnz, baseIds, Ngroups, groupOffsets);

  memory<dlong> rows(nnz);

  #pragma omp parallel for
  for (dlong n=0;n<Ngroups;++n) {
    const dlong start = groupOffsets[n];
    const dlong end = groupOffsets[n+1];

    const dlong lid = haloIds[n]; //location of this baseId in the halo

    for (dlong i=start;i<end;++i) {
      rows[sortIds[i]] = lid;
    }
  }

  /*Sort groups by their row*/
  memory<dlong> rowSortIds(nnz);
  prim::stableSort(nnz, rows, rowSortIds);

  memory<dlong> sortedCols(nnz);
  prim::transformGather(nnz, rowSortIds, cols, sortedCols);
  cols = sortedCols;

  /*Build the final gather op to assemble the recieved data*/
  levels[lvl].gather = ogsOperator_t(platform,
                                     Unsigned,
                                     NhaloT,
                                     NhaloT,
                                     Ncols,
                                     nnz,
                                     memory<hlong>(),
                                     rows,
                                     cols);
}

ogsCrystalRouter_t::ogsCrystalRouter_t(Kind kind,
                                       const dlong Nshared,
                                       const memory<int>   sharedRanks,
                                       const memory<dlong> sharedLocalRows,
                                       const memory<dlong> sharedRemoteRows,
                                       const memory<hlong> sharedLocalBaseIds,
                                       const memory<hlong> sharedRemoteBaseIds,
                                       const memory<hlong> haloBaseIds,
                                       ogsOperator_t& gatherHalo,
                                       stream_t _dataStream,
                                       comm_t _comm,
                                       platform_t &_platform):
  ogsExchange_t(_platform,_comm,_dataStream) {

  NhaloP = gatherHalo.NrowsN;
  Nhalo  = gatherHalo.NrowsT;

  //first count how many levels we need
  Nlevels = 0;
  int np = size;
  int np_offset=0;
  while (np>1) {
    int np_half = (np+1)/2;
    int r_half = np_half + np_offset;

    int is_lo = (rank<r_half) ? 1 : 0;

    //Shrink the size of the hypercube
    if (is_lo) {
      np = np_half;
    } else {
      np -= np_half;
      np_offset = r_half;
    }
    Nlevels++;
  }

  // Expand the list of shared data to include the current halo nodes

  dlong N = Nshared + Nhalo;

  memory<dlong> localRows(N);
  prim::range(Nhalo, 0, 1, localRows);
  localRows.copyFrom(sharedLocalRows, Nshared, Nhalo);

  memory<hlong> localBaseIds(N);
  localBaseIds.copyFrom(haloBaseIds, Nhalo);
  localBaseIds.copyFrom(sharedLocalBaseIds, Nshared, Nhalo);

  memory<hlong> remoteBaseIds(N);
  remoteBaseIds.copyFrom(haloBaseIds, Nhalo);
  remoteBaseIds.copyFrom(sharedRemoteBaseIds, Nshared, Nhalo);

  memory<int> remoteRanks(N);
  prim::set(Nhalo, rank, remoteRanks);
  remoteRanks.copyFrom(sharedRanks, Nshared, Nhalo);


  memory<hlong> absRemoteBaseIds(N);
  memory<dlong> baseSortIds(N);

  /*Sort list of shared nodes into baseId groups*/
  prim::abs(N, remoteBaseIds, absRemoteBaseIds);
  prim::stableSort(N, absRemoteBaseIds, baseSortIds);

  memory<dlong> tmpLocalRows(N);
  prim::transformGather(N, baseSortIds, localRows, tmpLocalRows);
  localRows = tmpLocalRows;

  memory<hlong> tmpLocalBaseIds(N);
  prim::transformGather(N, baseSortIds, localBaseIds, tmpLocalBaseIds);
  localBaseIds = tmpLocalBaseIds;

  memory<hlong> tmpRemoteBaseIds(N);
  prim::transformGather(N, baseSortIds, remoteBaseIds, tmpRemoteBaseIds);
  remoteBaseIds = tmpRemoteBaseIds;

  memory<int> tmpRemoteRanks(N);
  prim::transformGather(N, baseSortIds, remoteRanks, tmpRemoteRanks);
  remoteRanks = tmpRemoteRanks;

  baseSortIds.free();

  memory<hlong> absHaloBaseIds(Nhalo);
  prim::abs(Nhalo, haloBaseIds, absHaloBaseIds);

  // Sym     mode - Send everything, gather to all Nhalo nodes
  // NoTrans mode - Only send positive baseIds, gather to all Nhalo nodes
  // Trans   mode - Only send to remote positive baseIds, gather to positive NhaloP nodes

  /*Build the symmetric exchange using all the shared data*/
  dlong NsendS = N;
  data[Sym].setupExchange(Nlevels,
                          Nhalo,
                          Nhalo,
                          absHaloBaseIds,
                          NsendS,
                          remoteRanks,
                          absRemoteBaseIds,
                          localRows,
                          comm,
                          platform);

  if (kind==Signed) {
    /*NoTrans: Get locations of shared nodes that have a local positive baseId*/
    memory<int> noTransFlags(N);

    #pragma omp parallel for
    for (dlong n=0; n<N;++n) {
      noTransFlags[n] = (localBaseIds[n]>0) ? 1 : 0;
    }

    dlong NsendN = prim::count(N, noTransFlags, 1);
    memory<dlong> noTransIds(NsendN);
    prim::select(N, noTransFlags, 1, noTransIds);

    /*Extract the subset of the shared node list for these nodes*/
    memory<int> remoteRanksN(NsendN);
    prim::transformGather(NsendN, noTransIds, remoteRanks, remoteRanksN);

    memory<hlong> absRemoteBaseIdsN(NsendN);
    prim::transformGather(NsendN, noTransIds, absRemoteBaseIds, absRemoteBaseIdsN);

    memory<dlong> localRowsN(NsendN);
    prim::transformGather(NsendN, noTransIds, localRows, localRowsN);

    /*Build the NoTrans exchange*/
    data[NoTrans].setupExchange(Nlevels,
                                NhaloP,
                                Nhalo,
                                absHaloBaseIds,
                                NsendN,
                                remoteRanksN,
                                absRemoteBaseIdsN,
                                localRowsN,
                                comm,
                                platform);

    /*Trans: Get locations of shared nodes that have a remote positive baseId*/
    memory<int> transFlags(N);

    #pragma omp parallel for
    for (dlong n=0; n<N;++n) {
      transFlags[n] = (remoteBaseIds[n]>0) ? 1 : 0;
    }

    dlong NsendT = prim::count(N, transFlags, 1);
    memory<dlong> transIds(NsendT);
    prim::select(N, transFlags,  1, transIds);

    /*Extract the subset of the shared node list for these nodes*/
    memory<int> remoteRanksT(NsendT);
    prim::transformGather(NsendT, transIds, remoteRanks, remoteRanksT);

    memory<hlong> absRemoteBaseIdsT(NsendT);
    prim::transformGather(NsendT, transIds, absRemoteBaseIds, absRemoteBaseIdsT);

    memory<dlong> localRowsT(NsendT);
    prim::transformGather(NsendT, transIds, localRows, localRowsT);

    /*Build the Trans exchange*/
    data[Trans].setupExchange(Nlevels,
                              Nhalo,
                              NhaloP,
                              absHaloBaseIds,
                              NsendT,
                              remoteRanksT,
                              absRemoteBaseIdsT,
                              localRowsT,
                              comm,
                              platform);
  } else {
    data[NoTrans] = data[Sym];
    data[Trans] = data[Sym];
  }

  //make scratch space
  AllocBuffer(sizeof(dfloat));
}

void ogsCrystalRouter_t::AllocBuffer(size_t Nbytes) {

  if (o_workspace.byte_size() < data[Sym].NsendMax*Nbytes) {
    h_workspace = platform.hostMalloc<char>(data[Sym].NsendMax*Nbytes);
    o_workspace = platform.malloc<char>(data[Sym].NsendMax*Nbytes);
  }
  if (o_sendspace.byte_size() < data[Sym].NrecvMax*Nbytes) {
    h_sendspace = platform.hostMalloc<char>(data[Sym].NrecvMax*Nbytes);
    h_recvspace = platform.hostMalloc<char>(data[Sym].NrecvMax*Nbytes);
    o_sendspace = platform.malloc<char>(data[Sym].NrecvMax*Nbytes);
    o_recvspace = platform.malloc<char>(data[Sym].NrecvMax*Nbytes);
  }
}

} //namespace ogs

} //namespace libp
