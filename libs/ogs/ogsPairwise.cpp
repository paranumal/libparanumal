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
inline void ogsPairwise_t::HostStart(const int k, const Op op, const Transpose trans){

  pinnedMemory<T> sendBuf = h_sendspace;
  pinnedMemory<T> workBuf = h_workspace;

  // Get the exchange data based on the mode
  data_t& d = data[trans];

  // extract the send buffer
  extract(d.Nsend, k, d.sendIds, sendBuf, workBuf);

  //post recvs
  for (int r=0;r<d.NranksRecv;r++) {
    comm.Irecv(sendBuf + d.NrowsP*k + d.recvOffsets[r]*k,
               d.recvRanks[r],
               k*d.recvCounts[r],
               d.recvRanks[r],
               requests[r]);
  }

  //post sends
  for (int r=0;r<d.NranksSend;r++) {
    comm.Isend(workBuf + d.sendOffsets[r]*k,
              d.sendRanks[r],
              k*d.sendCounts[r],
              rank,
              requests[d.NranksRecv+r]);
  }
}

template<typename T>
inline void ogsPairwise_t::HostFinish(const int k, const Op op, const Transpose trans){

  pinnedMemory<T> sendBuf = h_sendspace;
  pinnedMemory<T> recvBuf = h_recvspace;

  // Get the exchange data based on the mode
  data_t& d = data[trans];

  comm.Waitall(d.NranksRecv+d.NranksSend, requests);

  // gather the recieved nodes
  d.postmpi.Gather(recvBuf, sendBuf, k, op, Sym);
}

void ogsPairwise_t::HostStart(const Type type, const int k, const Op op, const Transpose trans) {
  switch (type) {
    case Int32:  HostStart<int>(k, op, trans); break;
    case Int64:  HostStart<long long int>(k, op, trans); break;
    case Float:  HostStart<float>(k, op, trans); break;
    case Double: HostStart<double>(k, op, trans); break;
  }
}
void ogsPairwise_t::HostFinish(const Type type, const int k, const Op op, const Transpose trans) {
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
void ogsPairwise_t::DeviceStart(const int k, const Op op, const Transpose trans){

  deviceMemory<T> o_sendBuf = o_sendspace;
  deviceMemory<T> o_workBuf = o_workspace;

  // Get the exchange data based on the mode
  data_t& d = data[trans];

  //  assemble the send buffer on device
  extractKernel[ogsType<T>::get()](d.Nsend, k, d.o_sendIds, o_sendBuf, o_workBuf);

  //wait for kernel to finish on default stream
  if (d.Nsend) {
    device_t &device = platform.device;
    device.finish();
  }
}

template<typename T>
void ogsPairwise_t::DeviceFinish(const int k, const Op op, const Transpose trans){

  deviceMemory<T> o_sendBuf = o_sendspace;
  deviceMemory<T> o_workBuf = o_workspace;
  deviceMemory<T> o_recvBuf = o_recvspace;

  // Get the exchange data based on the mode
  data_t& d = data[trans];

  //post recvs
  for (int r=0;r<d.NranksRecv;r++) {
    comm.Irecv(o_sendBuf + d.NrowsP*k + d.recvOffsets[r]*k,
              d.recvRanks[r],
              k*d.recvCounts[r],
              d.recvRanks[r],
              requests[r]);
  }

  //post sends
  for (int r=0;r<d.NranksSend;r++) {
    comm.Isend(o_workBuf + d.sendOffsets[r]*k,
              d.sendRanks[r],
              k*d.sendCounts[r],
              rank,
              requests[d.NranksRecv+r]);
  }

  comm.Waitall(d.NranksRecv+d.NranksSend, requests);

  // gather the recieved nodes on device
  d.postmpi.Gather(o_recvBuf, o_sendBuf, k, op, Sym);
}

void ogsPairwise_t::DeviceStart(const Type type, const int k, const Op op, const Transpose trans) {
  switch (type) {
    case Int32:  DeviceStart<int>(k, op, trans); break;
    case Int64:  DeviceStart<long long int>(k, op, trans); break;
    case Float:  DeviceStart<float>(k, op, trans); break;
    case Double: DeviceStart<double>(k, op, trans); break;
  }
}
void ogsPairwise_t::DeviceFinish(const Type type, const int k, const Op op, const Transpose trans) {
  switch (type) {
    case Int32:  DeviceFinish<int>(k, op, trans); break;
    case Int64:  DeviceFinish<long long int>(k, op, trans); break;
    case Float:  DeviceFinish<float>(k, op, trans); break;
    case Double: DeviceFinish<double>(k, op, trans); break;
  }
}


void ogsPairwise_t::data_t::setupExchange(const dlong Nsend_,
                                          const dlong NrowsP_,
                                          const dlong Nrows_,
                                          const memory<int>   destRanks,
                                          const memory<dlong> localRows,
                                          const memory<dlong> remoteRows,
                                          comm_t comm,
                                          platform_t& platform) {

  int size = comm.size();

  Nsend = Nsend_;
  NrowsP = NrowsP_;
  Nrows = Nrows_;

  memory<int> mpiSendCounts(size);
  memory<int> mpiRecvCounts(size);
  memory<int> mpiSendOffsets(size+1);
  memory<int> mpiRecvOffsets(size+1);

  /*Get length and offsets of groups of destinations*/
  prim::runLengthEncodeConsecutive(Nsend, destRanks, size, mpiSendOffsets);
  prim::adjacentDifference(size, mpiSendOffsets+1, mpiSendCounts);

  comm.Alltoall(mpiSendCounts, mpiRecvCounts);
  mpiRecvOffsets[0] = 0;
  prim::inclusiveScan(size, mpiRecvCounts, mpiRecvOffsets+1);
  dlong Nrecv = mpiRecvOffsets[size]; //total ids to recv

  sendIds.malloc(Nsend);
  sendIds.copyFrom(localRows);
  o_sendIds = platform.malloc(sendIds);

  //send the node lists so we know what we'll receive and in what order
  Ncols = NrowsP + Nrecv;
  memory<dlong> rows(Ncols);
  memory<dlong> cols(Ncols);

  prim::range(NrowsP, 0, 1, rows);

  //Send list of rows to each rank
  comm.Alltoallv(remoteRows,  mpiSendCounts, mpiSendOffsets,
                 rows+NrowsP, mpiRecvCounts, mpiRecvOffsets);

  /*Sort groups by their row*/
  prim::stableSort(Ncols, rows, cols);

  /*Build the gather op to assemble the recieved data from MPI*/
  postmpi = ogsOperator_t(platform,
                          Unsigned,
                          Nrows,
                          Nrows,
                          Ncols,
                          Ncols,
                          memory<hlong>(),
                          rows,
                          cols);

  //compress the send/recv counts to pairwise exchanges
  memory<int> sendFlag(size);
  memory<int> recvFlag(size);

  #pragma omp parallel for
  for (int r=0;r<size;++r) { sendFlag[r] = (mpiSendCounts[r]>0) ? 1 : 0; }

  #pragma omp parallel for
  for (int r=0;r<size;++r) { recvFlag[r] = (mpiRecvCounts[r]>0) ? 1 : 0; }

  NranksSend = prim::count(size, sendFlag, 1);
  NranksRecv = prim::count(size, recvFlag, 1);

  sendRanks.malloc(NranksSend);
  recvRanks.malloc(NranksRecv);
  sendCounts.malloc(NranksSend);
  recvCounts.malloc(NranksRecv);
  sendOffsets.malloc(NranksSend);
  recvOffsets.malloc(NranksRecv);

  prim::select(size, sendFlag, 1, sendRanks);
  prim::select(size, recvFlag, 1, recvRanks);

  prim::transformGather(NranksSend, sendRanks, mpiSendCounts,  sendCounts);
  prim::transformGather(NranksRecv, recvRanks, mpiRecvCounts,  recvCounts);
  prim::transformGather(NranksSend, sendRanks, mpiSendOffsets, sendOffsets);
  prim::transformGather(NranksRecv, recvRanks, mpiRecvOffsets, recvOffsets);
}

ogsPairwise_t::ogsPairwise_t(Kind kind,
                             const dlong Nshared,
                             const memory<int>   sharedRanks,
                             const memory<dlong> sharedLocalRows,
                             const memory<dlong> sharedRemoteRows,
                             const memory<hlong> sharedLocalBaseIds,
                             const memory<hlong> sharedRemoteBaseIds,
                             ogsOperator_t& gatherHalo,
                             stream_t _dataStream,
                             comm_t _comm,
                             platform_t &_platform):
  ogsExchange_t(_platform,_comm,_dataStream) {

  Nhalo  = gatherHalo.NrowsT;
  NhaloP = gatherHalo.NrowsN;

  memory<int>   destRanks(Nshared);
  memory<dlong> destRankSortIds(Nshared);
  destRanks.copyFrom(sharedRanks);

  /*Sort list of shared nodes by rank to get the ordering for MPI Alltoallv*/
  prim::stableSort(Nshared, destRanks, destRankSortIds);

  memory<dlong> localRows(Nshared);
  prim::transformGather(Nshared, destRankSortIds, sharedLocalRows, localRows);

  memory<dlong> remoteRows(Nshared);
  prim::transformGather(Nshared, destRankSortIds, sharedRemoteRows, remoteRows);

  // Sym     mode - Send everything, gather to all Nhalo nodes
  // NoTrans mode - Only send positive baseIds, gather to all Nhalo nodes
  // Trans   mode - Only send to remote positive baseIds, gather to positive NhaloP nodes

  /*Build the symmetric exchange using all the shared data*/
  dlong NsendS = Nshared;
  data[Sym].setupExchange(NsendS,
                          Nhalo,
                          Nhalo,
                          destRanks,
                          localRows,
                          remoteRows,
                          comm,
                          platform);

  if (kind==Signed) {
    /*Reorder baseId info*/
    memory<hlong> localBaseIds(Nshared);
    prim::transformGather(Nshared, destRankSortIds, sharedLocalBaseIds, localBaseIds);

    memory<hlong> remoteBaseIds(Nshared);
    prim::transformGather(Nshared, destRankSortIds, sharedRemoteBaseIds, remoteBaseIds);

    /*NoTrans: Get locations of shared nodes that have a local positive baseId*/
    memory<int> noTransFlags(Nshared);

    #pragma omp parallel for
    for (dlong n=0; n<Nshared;++n) {
      noTransFlags[n] = (localBaseIds[n]>0) ? 1 : 0;
    }

    dlong NsendN = prim::count(Nshared, noTransFlags, 1);
    memory<dlong> noTransIds(NsendN);
    prim::select(Nshared, noTransFlags, 1, noTransIds);

    /*Extract the subset of the shared node list for these nodes*/
    memory<int> destRanksN(NsendN);
    prim::transformGather(NsendN, noTransIds, destRanks, destRanksN);

    memory<dlong> localRowsN(NsendN);
    prim::transformGather(NsendN, noTransIds, localRows, localRowsN);

    memory<dlong> remoteRowsN(NsendN);
    prim::transformGather(NsendN, noTransIds, remoteRows, remoteRowsN);

    /*Build the NoTrans exchange*/
    data[NoTrans].setupExchange(NsendN,
                                NhaloP,
                                Nhalo,
                                destRanksN,
                                localRowsN,
                                remoteRowsN,
                                comm,
                                platform);

    /*Trans: Get locations of shared nodes that have a remote positive baseId*/
    memory<int> transFlags(Nshared);

    #pragma omp parallel for
    for (dlong n=0; n<Nshared;++n) {
      transFlags[n] = (remoteBaseIds[n]>0) ? 1 : 0;
    }

    dlong NsendT = prim::count(Nshared, transFlags, 1);
    memory<dlong> transIds(NsendT);
    prim::select(Nshared, transFlags,  1, transIds);

    /*Extract the subset of the shared node list for these nodes*/
    memory<int> destRanksT(NsendT);
    prim::transformGather(NsendT, transIds, destRanks, destRanksT);

    memory<dlong> localRowsT(NsendT);
    prim::transformGather(NsendT, transIds, localRows, localRowsT);

    memory<dlong> remoteRowsT(NsendT);
    prim::transformGather(NsendT, transIds, remoteRows, remoteRowsT);

    /*Build the Trans exchange*/
    data[Trans].setupExchange(NsendT,
                              NhaloP,
                              NhaloP,
                              destRanksT,
                              localRowsT,
                              remoteRowsT,
                              comm,
                              platform);
  } else {
    data[NoTrans] = data[Sym];
    data[Trans] = data[Sym];
  }

  requests.malloc(data[Sym].NranksSend+data[Sym].NranksRecv);

  //make scratch space
  AllocBuffer(sizeof(dfloat));
}

void ogsPairwise_t::AllocBuffer(size_t Nbytes) {
  if (o_sendspace.byte_size() < data[Sym].Ncols*Nbytes) {
    h_sendspace = platform.hostMalloc<char>(data[Sym].Ncols*Nbytes);
    o_sendspace = platform.malloc<char>(data[Sym].Ncols*Nbytes);
  }
  if (o_workspace.byte_size() < data[Sym].Nsend*Nbytes) {
    h_workspace = platform.hostMalloc<char>(data[Sym].Nsend*Nbytes);
    o_workspace = platform.malloc<char>(data[Sym].Nsend*Nbytes);
  }
  if (o_recvspace.byte_size() < data[Sym].Nrows*Nbytes) {
    h_recvspace = platform.hostMalloc<char>(data[Sym].Nrows*Nbytes);
    o_recvspace = platform.malloc<char>(data[Sym].Nrows*Nbytes);
  }
}

} //namespace ogs

} //namespace libp
