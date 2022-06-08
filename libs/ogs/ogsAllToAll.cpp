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

#ifdef GLIBCXX_PARALLEL
#include <parallel/algorithm>
using __gnu_parallel::sort;
#else
using std::sort;
#endif

namespace libp {

namespace ogs {

/**********************************
* Host exchange
***********************************/
template<typename T>
inline void ogsAllToAll_t::Start(pinnedMemory<T> &buf, const int k,
                          const Op op, const Transpose trans){

  pinnedMemory<T> sendBuf = h_sendspace;

  // extract the send buffer
  if (trans == NoTrans)
    extract(NsendN, k, sendIdsN, buf, sendBuf);
  else
    extract(NsendT, k, sendIdsT, buf, sendBuf);

  if (trans==NoTrans) {
    for (int r=0;r<size;++r) {
      sendCounts[r] = k*mpiSendCountsN[r];
      recvCounts[r] = k*mpiRecvCountsN[r];
      sendOffsets[r+1] = k*mpiSendOffsetsN[r+1];
      recvOffsets[r+1] = k*mpiRecvOffsetsN[r+1];
    }
  } else {
    for (int r=0;r<size;++r) {
      sendCounts[r] = k*mpiSendCountsT[r];
      recvCounts[r] = k*mpiRecvCountsT[r];
      sendOffsets[r+1] = k*mpiSendOffsetsT[r+1];
      recvOffsets[r+1] = k*mpiRecvOffsetsT[r+1];
    }
  }

  // collect everything needed with single MPI all to all
  comm.Ialltoallv(sendBuf,     sendCounts, sendOffsets,
                  buf+Nhalo*k, recvCounts, recvOffsets,
                  request);
}

template<typename T>
inline void ogsAllToAll_t::Finish(pinnedMemory<T> &buf, const int k,
                           const Op op, const Transpose trans){

  comm.Wait(request);

  //if we recvieved anything via MPI, gather the recv buffer and scatter
  // it back to to original vector
  dlong Nrecv = recvOffsets[size];
  if (Nrecv) {
    // gather the recieved nodes
    postmpi.Gather(buf, buf, k, op, trans);
  }
}

void ogsAllToAll_t::Start(pinnedMemory<float> &buf, const int k, const Op op, const Transpose trans) { Start<float>(buf, k, op, trans); }
void ogsAllToAll_t::Start(pinnedMemory<double> &buf, const int k, const Op op, const Transpose trans) { Start<double>(buf, k, op, trans); }
void ogsAllToAll_t::Start(pinnedMemory<int> &buf, const int k, const Op op, const Transpose trans) { Start<int>(buf, k, op, trans); }
void ogsAllToAll_t::Start(pinnedMemory<long long int> &buf, const int k, const Op op, const Transpose trans) { Start<long long int>(buf, k, op, trans); }
void ogsAllToAll_t::Finish(pinnedMemory<float> &buf, const int k, const Op op, const Transpose trans) { Finish<float>(buf, k, op, trans); }
void ogsAllToAll_t::Finish(pinnedMemory<double> &buf, const int k, const Op op, const Transpose trans) { Finish<double>(buf, k, op, trans); }
void ogsAllToAll_t::Finish(pinnedMemory<int> &buf, const int k, const Op op, const Transpose trans) { Finish<int>(buf, k, op, trans); }
void ogsAllToAll_t::Finish(pinnedMemory<long long int> &buf, const int k, const Op op, const Transpose trans) { Finish<long long int>(buf, k, op, trans); }

/**********************************
* GPU-aware exchange
***********************************/
template<typename T>
void ogsAllToAll_t::Start(deviceMemory<T> &o_buf,
                          const int k,
                          const Op op,
                          const Transpose trans){

  const dlong Nsend = (trans == NoTrans) ? NsendN : NsendT;

  if (Nsend) {
    deviceMemory<T> o_sendBuf = o_sendspace;

    // assemble the send buffer on device
    if (trans == NoTrans) {
      extractKernel[ogsType<T>::get()](NsendN, k, o_sendIdsN, o_buf, o_sendBuf);
    } else {
      extractKernel[ogsType<T>::get()](NsendT, k, o_sendIdsT, o_buf, o_sendBuf);
    }
    //wait for kernel to finish on default stream
    device_t &device = platform.device;
    device.finish();
  }
}

template<typename T>
void ogsAllToAll_t::Finish(deviceMemory<T> &o_buf,
                           const int k,
                           const Op op,
                           const Transpose trans){

  deviceMemory<T> o_sendBuf = o_sendspace;

  if (trans==NoTrans) {
    for (int r=0;r<size;++r) {
      sendCounts[r] = k*mpiSendCountsN[r];
      recvCounts[r] = k*mpiRecvCountsN[r];
      sendOffsets[r+1] = k*mpiSendOffsetsN[r+1];
      recvOffsets[r+1] = k*mpiRecvOffsetsN[r+1];
    }
  } else {
    for (int r=0;r<size;++r) {
      sendCounts[r] = k*mpiSendCountsT[r];
      recvCounts[r] = k*mpiRecvCountsT[r];
      sendOffsets[r+1] = k*mpiSendOffsetsT[r+1];
      recvOffsets[r+1] = k*mpiRecvOffsetsT[r+1];
    }
  }

  // collect everything needed with single MPI all to all
  comm.Alltoallv(o_sendBuf,     sendCounts, sendOffsets,
                 o_buf+Nhalo*k, recvCounts, recvOffsets);

  //if we recvieved anything via MPI, gather the recv buffer and scatter
  // it back to to original vector
  dlong Nrecv = recvOffsets[size];
  if (Nrecv) {
    // gather the recieved nodes on device
    postmpi.Gather(o_buf, o_buf, k, op, trans);
  }
}

void ogsAllToAll_t::Start(deviceMemory<float> &buf, const int k, const Op op, const Transpose trans) { Start<float>(buf, k, op, trans); }
void ogsAllToAll_t::Start(deviceMemory<double> &buf, const int k, const Op op, const Transpose trans) { Start<double>(buf, k, op, trans); }
void ogsAllToAll_t::Start(deviceMemory<int> &buf, const int k, const Op op, const Transpose trans) { Start<int>(buf, k, op, trans); }
void ogsAllToAll_t::Start(deviceMemory<long long int> &buf, const int k, const Op op, const Transpose trans) { Start<long long int>(buf, k, op, trans); }
void ogsAllToAll_t::Finish(deviceMemory<float> &buf, const int k, const Op op, const Transpose trans) { Finish<float>(buf, k, op, trans); }
void ogsAllToAll_t::Finish(deviceMemory<double> &buf, const int k, const Op op, const Transpose trans) { Finish<double>(buf, k, op, trans); }
void ogsAllToAll_t::Finish(deviceMemory<int> &buf, const int k, const Op op, const Transpose trans) { Finish<int>(buf, k, op, trans); }
void ogsAllToAll_t::Finish(deviceMemory<long long int> &buf, const int k, const Op op, const Transpose trans) { Finish<long long int>(buf, k, op, trans); }

ogsAllToAll_t::ogsAllToAll_t(dlong Nshared,
                             memory<parallelNode_t> &sharedNodes,
                             ogsOperator_t& gatherHalo,
                             stream_t _dataStream,
                             comm_t _comm,
                             platform_t &_platform):
  ogsExchange_t(_platform,_comm, _dataStream) {

  Nhalo  = gatherHalo.NrowsT;
  NhaloP = gatherHalo.NrowsN;

  // sort the list by rank to the order where they will be sent by MPI_Allgatherv
  sort(sharedNodes.ptr(), sharedNodes.ptr()+Nshared,
       [](const parallelNode_t& a, const parallelNode_t& b) {
         if(a.rank < b.rank) return true; //group by rank
         if(a.rank > b.rank) return false;

         return a.newId < b.newId; //then order by the localId relative to this rank
       });

  //make mpi allgatherv counts and offsets
  mpiSendCountsT.calloc(size);
  mpiSendCountsN.calloc(size);
  mpiRecvCountsT.malloc(size);
  mpiRecvCountsN.malloc(size);
  mpiSendOffsetsT.malloc(size+1);
  mpiSendOffsetsN.malloc(size+1);
  mpiRecvOffsetsN.malloc(size+1);
  mpiRecvOffsetsT.malloc(size+1);

  for (dlong n=0;n<Nshared;n++) { //loop through nodes we need to send
    const int r = sharedNodes[n].rank;
    if (sharedNodes[n].sign>0) mpiSendCountsN[r]++;
    mpiSendCountsT[r]++;
  }

  //shared counts
  comm.Alltoall(mpiSendCountsT, mpiRecvCountsT);
  comm.Alltoall(mpiSendCountsN, mpiRecvCountsN);

  //cumulative sum
  mpiSendOffsetsN[0] = 0;
  mpiSendOffsetsT[0] = 0;
  mpiRecvOffsetsN[0] = 0;
  mpiRecvOffsetsT[0] = 0;
  for (int r=0;r<size;r++) {
    mpiSendOffsetsN[r+1] = mpiSendOffsetsN[r]+mpiSendCountsN[r];
    mpiSendOffsetsT[r+1] = mpiSendOffsetsT[r]+mpiSendCountsT[r];
    mpiRecvOffsetsN[r+1] = mpiRecvOffsetsN[r]+mpiRecvCountsN[r];
    mpiRecvOffsetsT[r+1] = mpiRecvOffsetsT[r]+mpiRecvCountsT[r];
  }

  //make ops for scattering halo nodes before sending
  NsendN=mpiSendOffsetsN[size];
  NsendT=mpiSendOffsetsT[size];

  sendIdsN.malloc(NsendN);
  sendIdsT.malloc(NsendT);

  NsendN=0; //positive node count
  NsendT=0; //all node count

  for (dlong n=0;n<Nshared;n++) { //loop through nodes we need to send
    dlong id = sharedNodes[n].newId; //coalesced index for this baseId on this rank
    if (sharedNodes[n].sign==2) {
      sendIdsN[NsendN++] = id;
    }
    sendIdsT[NsendT++] = id;
  }
  o_sendIdsT = platform.malloc(sendIdsT);
  o_sendIdsN = platform.malloc(sendIdsN);

  //send the node lists so we know what we'll receive
  dlong Nrecv = mpiRecvOffsetsT[size];
  memory<parallelNode_t> recvNodes(Nrecv);

  //Send list of nodes to each rank
  comm.Alltoallv(sharedNodes, mpiSendCountsT, mpiSendOffsetsT,
                   recvNodes, mpiRecvCountsT, mpiRecvOffsetsT);

  //make ops for gathering halo nodes after an MPI_Allgatherv
  postmpi.platform = platform;
  postmpi.kind = Signed;

  postmpi.NrowsN = Nhalo;
  postmpi.NrowsT = Nhalo;
  postmpi.rowStartsN.malloc(Nhalo+1);
  postmpi.rowStartsT.malloc(Nhalo+1);

  //make array of counters
  memory<dlong> haloGatherTCounts(Nhalo);
  memory<dlong> haloGatherNCounts(Nhalo);

  //count the data that will already be in h_haloBuf.ptr()
  for (dlong n=0;n<Nhalo;n++) {
    haloGatherNCounts[n] = (n<NhaloP) ? 1 : 0;
    haloGatherTCounts[n] = 1;
  }

  for (dlong n=0;n<Nrecv;n++) { //loop through nodes needed for gathering halo nodes
    dlong id = recvNodes[n].localId; //coalesced index for this baseId on this rank
    if (recvNodes[n].sign==2) haloGatherNCounts[id]++;  //tally
    haloGatherTCounts[id]++;  //tally
  }

  postmpi.rowStartsN[0] = 0;
  postmpi.rowStartsT[0] = 0;
  for (dlong i=0;i<Nhalo;i++) {
    postmpi.rowStartsN[i+1] = postmpi.rowStartsN[i] + haloGatherNCounts[i];
    postmpi.rowStartsT[i+1] = postmpi.rowStartsT[i] + haloGatherTCounts[i];
    haloGatherNCounts[i] = 0;
    haloGatherTCounts[i] = 0;
  }
  postmpi.nnzN = postmpi.rowStartsN[Nhalo];
  postmpi.nnzT = postmpi.rowStartsT[Nhalo];
  postmpi.colIdsN.malloc(postmpi.nnzN);
  postmpi.colIdsT.malloc(postmpi.nnzT);

  for (dlong n=0;n<NhaloP;n++) {
    const dlong soffset = postmpi.rowStartsN[n];
    const int sindex  = haloGatherNCounts[n];
    postmpi.colIdsN[soffset+sindex] = n; //record id
    haloGatherNCounts[n]++;
  }
  for (dlong n=0;n<Nhalo;n++) {
    const dlong soffset = postmpi.rowStartsT[n];
    const int sindex  = haloGatherTCounts[n];
    postmpi.colIdsT[soffset+sindex] = n; //record id
    haloGatherTCounts[n]++;
  }

  dlong cnt=Nhalo; //positive node count
  for (dlong n=0;n<Nrecv;n++) { //loop through nodes we need to send
    dlong id = recvNodes[n].localId; //coalesced index for this baseId on this rank
    if (recvNodes[n].sign==2) {
      const dlong soffset = postmpi.rowStartsN[id];
      const int sindex  = haloGatherNCounts[id];
      postmpi.colIdsN[soffset+sindex] = cnt++; //record id
      haloGatherNCounts[id]++;
    }
    const dlong soffset = postmpi.rowStartsT[id];
    const int sindex  = haloGatherTCounts[id];
    postmpi.colIdsT[soffset+sindex] = n + Nhalo; //record id
    haloGatherTCounts[id]++;
  }

  postmpi.o_rowStartsN = platform.malloc(postmpi.rowStartsN);
  postmpi.o_rowStartsT = platform.malloc(postmpi.rowStartsT);
  postmpi.o_colIdsN = platform.malloc(postmpi.colIdsN);
  postmpi.o_colIdsT = platform.malloc(postmpi.colIdsT);

  //free up space
  recvNodes.free();
  haloGatherNCounts.free();
  haloGatherTCounts.free();

  postmpi.setupRowBlocks();

  sendCounts.malloc(size);
  recvCounts.malloc(size);
  sendOffsets.malloc(size+1);
  recvOffsets.malloc(size+1);

  sendOffsets[0]=0;
  recvOffsets[0]=0;

  //make scratch space
  AllocBuffer(sizeof(dfloat));
}

void ogsAllToAll_t::AllocBuffer(size_t Nbytes) {
  if (o_workspace.size() < postmpi.nnzT*Nbytes) {
    h_workspace = platform.hostMalloc<char>(postmpi.nnzT*Nbytes);
    o_workspace = platform.malloc<char>(postmpi.nnzT*Nbytes);
  }
  if (o_sendspace.size() < NsendT*Nbytes) {
    h_sendspace = platform.hostMalloc<char>(NsendT*Nbytes);
    o_sendspace = platform.malloc<char>(NsendT*Nbytes);
  }
}

} //namespace ogs

} //namespace libp
