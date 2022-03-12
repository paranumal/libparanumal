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
inline void ogsPairwise_t::Start(pinnedMemory<T> &buf, const int k,
                          const Op op, const Transpose trans){

  pinnedMemory<T> sendBuf = h_sendspace;

  const int NranksSend  = (trans==NoTrans) ? NranksSendN  : NranksSendT;
  const int NranksRecv  = (trans==NoTrans) ? NranksRecvN  : NranksRecvT;
  const int *sendRanks  = (trans==NoTrans) ? sendRanksN.ptr()   : sendRanksT.ptr();
  const int *recvRanks  = (trans==NoTrans) ? recvRanksN.ptr()   : recvRanksT.ptr();
  const int *sendCounts = (trans==NoTrans) ? sendCountsN.ptr()  : sendCountsT.ptr();
  const int *recvCounts = (trans==NoTrans) ? recvCountsN.ptr()  : recvCountsT.ptr();
  const int *sendOffsets= (trans==NoTrans) ? sendOffsetsN.ptr() : sendOffsetsT.ptr();
  const int *recvOffsets= (trans==NoTrans) ? recvOffsetsN.ptr() : recvOffsetsT.ptr();

  //post recvs
  for (int r=0;r<NranksRecv;r++) {
    comm.Irecv(buf + Nhalo*k + recvOffsets[r]*k,
               recvRanks[r],
               k*recvCounts[r],
               recvRanks[r],
               requests[r]);
  }

  // extract the send buffer
  if (trans == NoTrans)
    extract(NsendN, k, sendIdsN, buf, sendBuf);
  else
    extract(NsendT, k, sendIdsT, buf, sendBuf);

  //post sends
  for (int r=0;r<NranksSend;r++) {
    comm.Isend(sendBuf + sendOffsets[r]*k,
              sendRanks[r],
              k*sendCounts[r],
              rank,
              requests[NranksRecv+r]);
  }
}

template<typename T>
inline void ogsPairwise_t::Finish(pinnedMemory<T> &buf, const int k,
                           const Op op, const Transpose trans){

  const int NranksSend  = (trans==NoTrans) ? NranksSendN  : NranksSendT;
  const int NranksRecv  = (trans==NoTrans) ? NranksRecvN  : NranksRecvT;
  const int *recvOffsets= (trans==NoTrans) ? recvOffsetsN.ptr() : recvOffsetsT.ptr();

  comm.Waitall(NranksRecv+NranksSend, requests);

  //if we recvieved anything via MPI, gather the recv buffer and scatter
  // it back to to original vector
  dlong Nrecv = recvOffsets[NranksRecv];
  if (Nrecv) {
    // gather the recieved nodes
    postmpi.Gather(buf, buf, k, op, trans);
  }
}

void ogsPairwise_t::Start(pinnedMemory<float> &buf, const int k, const Op op, const Transpose trans) { Start<float>(buf, k, op, trans); }
void ogsPairwise_t::Start(pinnedMemory<double> &buf, const int k, const Op op, const Transpose trans) { Start<double>(buf, k, op, trans); }
void ogsPairwise_t::Start(pinnedMemory<int> &buf, const int k, const Op op, const Transpose trans) { Start<int>(buf, k, op, trans); }
void ogsPairwise_t::Start(pinnedMemory<long long int> &buf, const int k, const Op op, const Transpose trans) { Start<long long int>(buf, k, op, trans); }
void ogsPairwise_t::Finish(pinnedMemory<float> &buf, const int k, const Op op, const Transpose trans) { Finish<float>(buf, k, op, trans); }
void ogsPairwise_t::Finish(pinnedMemory<double> &buf, const int k, const Op op, const Transpose trans) { Finish<double>(buf, k, op, trans); }
void ogsPairwise_t::Finish(pinnedMemory<int> &buf, const int k, const Op op, const Transpose trans) { Finish<int>(buf, k, op, trans); }
void ogsPairwise_t::Finish(pinnedMemory<long long int> &buf, const int k, const Op op, const Transpose trans) { Finish<long long int>(buf, k, op, trans); }

/**********************************
* GPU-aware exchange
***********************************/
template<typename T>
void ogsPairwise_t::Start(deviceMemory<T> &o_buf,
                          const int k,
                          const Op op,
                          const Transpose trans){

  const dlong Nsend = (trans == NoTrans) ? NsendN : NsendT;

  if (Nsend) {
    deviceMemory<T> o_sendBuf = o_sendspace;

    //  assemble the send buffer on device
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
void ogsPairwise_t::Finish(deviceMemory<T> &o_buf,
                           const int k,
                           const Op op,
                           const Transpose trans){

  deviceMemory<T> o_sendBuf = o_sendspace;

  const int NranksSend  = (trans==NoTrans) ? NranksSendN  : NranksSendT;
  const int NranksRecv  = (trans==NoTrans) ? NranksRecvN  : NranksRecvT;
  const int *sendRanks  = (trans==NoTrans) ? sendRanksN.ptr()   : sendRanksT.ptr();
  const int *recvRanks  = (trans==NoTrans) ? recvRanksN.ptr()   : recvRanksT.ptr();
  const int *sendCounts = (trans==NoTrans) ? sendCountsN.ptr()  : sendCountsT.ptr();
  const int *recvCounts = (trans==NoTrans) ? recvCountsN.ptr()  : recvCountsT.ptr();
  const int *sendOffsets= (trans==NoTrans) ? sendOffsetsN.ptr() : sendOffsetsT.ptr();
  const int *recvOffsets= (trans==NoTrans) ? recvOffsetsN.ptr() : recvOffsetsT.ptr();

  //post recvs
  for (int r=0;r<NranksRecv;r++) {
    comm.Irecv(o_buf + Nhalo*k + recvOffsets[r]*k,
              recvRanks[r],
              k*recvCounts[r],
              recvRanks[r],
              requests[r]);
  }

  //post sends
  for (int r=0;r<NranksSend;r++) {
    comm.Isend(o_sendBuf + sendOffsets[r]*k,
              sendRanks[r],
              k*sendCounts[r],
              rank,
              requests[NranksRecv+r]);
  }

  comm.Waitall(NranksRecv+NranksSend, requests);

  //if we recvieved anything via MPI, gather the recv buffer and scatter
  // it back to to original vector
  dlong Nrecv = recvOffsets[NranksRecv];
  if (Nrecv) {
    // gather the recieved nodes on device
    postmpi.Gather(o_buf, o_buf, k, op, trans);
  }
}

void ogsPairwise_t::Start(deviceMemory<float> &buf, const int k, const Op op, const Transpose trans) { Start<float>(buf, k, op, trans); }
void ogsPairwise_t::Start(deviceMemory<double> &buf, const int k, const Op op, const Transpose trans) { Start<double>(buf, k, op, trans); }
void ogsPairwise_t::Start(deviceMemory<int> &buf, const int k, const Op op, const Transpose trans) { Start<int>(buf, k, op, trans); }
void ogsPairwise_t::Start(deviceMemory<long long int> &buf, const int k, const Op op, const Transpose trans) { Start<long long int>(buf, k, op, trans); }
void ogsPairwise_t::Finish(deviceMemory<float> &buf, const int k, const Op op, const Transpose trans) { Finish<float>(buf, k, op, trans); }
void ogsPairwise_t::Finish(deviceMemory<double> &buf, const int k, const Op op, const Transpose trans) { Finish<double>(buf, k, op, trans); }
void ogsPairwise_t::Finish(deviceMemory<int> &buf, const int k, const Op op, const Transpose trans) { Finish<int>(buf, k, op, trans); }
void ogsPairwise_t::Finish(deviceMemory<long long int> &buf, const int k, const Op op, const Transpose trans) { Finish<long long int>(buf, k, op, trans); }

ogsPairwise_t::ogsPairwise_t(dlong Nshared,
                             memory<parallelNode_t> &sharedNodes,
                             ogsOperator_t& gatherHalo,
                             stream_t _dataStream,
                             comm_t _comm,
                             platform_t &_platform):
  ogsExchange_t(_platform,_comm,_dataStream) {

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
  memory<int> mpiSendCountsT(size,0);
  memory<int> mpiSendCountsN(size,0);
  memory<int> mpiRecvCountsT(size);
  memory<int> mpiRecvCountsN(size);
  memory<int> mpiSendOffsetsT(size+1);
  memory<int> mpiSendOffsetsN(size+1);
  memory<int> mpiRecvOffsetsT(size+1);
  memory<int> mpiRecvOffsetsN(size+1);

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

  sendIdsN.calloc(NsendN);
  sendIdsT.calloc(NsendT);

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
  postmpi.rowStartsN.calloc(Nhalo+1);
  postmpi.rowStartsT.calloc(Nhalo+1);

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

  for (dlong i=0;i<Nhalo;i++) {
    postmpi.rowStartsN[i+1] = postmpi.rowStartsN[i] + haloGatherNCounts[i];
    postmpi.rowStartsT[i+1] = postmpi.rowStartsT[i] + haloGatherTCounts[i];
    haloGatherNCounts[i] = 0;
    haloGatherTCounts[i] = 0;
  }
  postmpi.nnzN = postmpi.rowStartsN[Nhalo];
  postmpi.nnzT = postmpi.rowStartsT[Nhalo];
  postmpi.colIdsN.calloc(postmpi.nnzN);
  postmpi.colIdsT.calloc(postmpi.nnzT);

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

  //compress the send/recv counts to pairwise exchanges
  NranksSendN=0;
  NranksSendT=0;
  NranksRecvN=0;
  NranksRecvT=0;
  for (int r=0;r<size;r++) {
    NranksSendN += (mpiSendCountsN[r]>0) ? 1 : 0;
    NranksSendT += (mpiSendCountsT[r]>0) ? 1 : 0;
    NranksRecvN += (mpiRecvCountsN[r]>0) ? 1 : 0;
    NranksRecvT += (mpiRecvCountsT[r]>0) ? 1 : 0;
  }

  sendRanksN.calloc(NranksSendN);
  sendRanksT.calloc(NranksSendT);
  recvRanksN.calloc(NranksRecvN);
  recvRanksT.calloc(NranksRecvT);
  sendCountsN.calloc(NranksSendN);
  sendCountsT.calloc(NranksSendT);
  recvCountsN.calloc(NranksRecvN);
  recvCountsT.calloc(NranksRecvT);
  sendOffsetsN.calloc(NranksSendN+1);
  sendOffsetsT.calloc(NranksSendT+1);
  recvOffsetsN.calloc(NranksRecvN+1);
  recvOffsetsT.calloc(NranksRecvT+1);

  //reset
  NranksSendN=0;
  NranksSendT=0;
  NranksRecvN=0;
  NranksRecvT=0;
  for (int r=0;r<size;r++) {
    if (mpiSendCountsN[r]>0) {
      sendRanksN[NranksSendN]  = r;
      sendCountsN[NranksSendN] = mpiSendCountsN[r];
      sendOffsetsN[NranksSendN] = mpiSendOffsetsN[r];
      NranksSendN++;
    }
    if (mpiSendCountsT[r]>0) {
      sendRanksT[NranksSendT]  = r;
      sendCountsT[NranksSendT] = mpiSendCountsT[r];
      sendOffsetsT[NranksSendT] = mpiSendOffsetsT[r];
      NranksSendT++;
    }
    if (mpiRecvCountsN[r]>0) {
      recvRanksN[NranksRecvN]   = r;
      recvCountsN[NranksRecvN]  = mpiRecvCountsN[r];
      recvOffsetsN[NranksRecvN] = mpiRecvOffsetsN[r];
      NranksRecvN++;
    }
    if (mpiRecvCountsT[r]>0) {
      recvRanksT[NranksRecvT]   = r;
      recvCountsT[NranksRecvT]  = mpiRecvCountsT[r];
      recvOffsetsT[NranksRecvT] = mpiRecvOffsetsT[r];
      NranksRecvT++;
    }
  }
  sendOffsetsN[NranksSendN] = mpiSendOffsetsN[size];
  sendOffsetsT[NranksSendT] = mpiSendOffsetsT[size];
  recvOffsetsN[NranksRecvN] = mpiRecvOffsetsN[size];
  recvOffsetsT[NranksRecvT] = mpiRecvOffsetsT[size];

  requests.malloc(NranksSendT+NranksRecvT);

  //make scratch space
  AllocBuffer(sizeof(dfloat));
}

void ogsPairwise_t::AllocBuffer(size_t Nbytes) {
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
