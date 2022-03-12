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

#ifndef OGS_EXCHANGE_HPP
#define OGS_EXCHANGE_HPP

#include "ogs.hpp"
#include "ogs/ogsOperator.hpp"

namespace libp {

namespace ogs {

//virtual base class to perform MPI exchange of gatherScatter
class ogsExchange_t {
public:
  platform_t platform;
  comm_t comm;
  int rank, size;

  dlong Nhalo, NhaloP;

  pinnedMemory<char> h_workspace, h_sendspace;
  deviceMemory<char> o_workspace, o_sendspace;

  stream_t dataStream;
  static kernel_t extractKernel[4];

#ifdef GPU_AWARE_MPI
  bool gpu_aware=true;
#else
  bool gpu_aware=false;
#endif

  ogsExchange_t(platform_t &_platform, comm_t _comm,
                stream_t _datastream):
    platform(_platform),
    comm(_comm),
    dataStream(_datastream) {
    rank = comm.rank();
    size = comm.size();
  }
  virtual ~ogsExchange_t() {}

  virtual void Start(pinnedMemory<float> &buf,const int k,const Op op,const Transpose trans)=0;
  virtual void Start(pinnedMemory<double> &buf,const int k,const Op op,const Transpose trans)=0;
  virtual void Start(pinnedMemory<int> &buf,const int k,const Op op,const Transpose trans)=0;
  virtual void Start(pinnedMemory<long long int> &buf,const int k,const Op op,const Transpose trans)=0;
  virtual void Finish(pinnedMemory<float> &buf,const int k,const Op op,const Transpose trans)=0;
  virtual void Finish(pinnedMemory<double> &buf,const int k,const Op op,const Transpose trans)=0;
  virtual void Finish(pinnedMemory<int> &buf,const int k,const Op op,const Transpose trans)=0;
  virtual void Finish(pinnedMemory<long long int> &buf,const int k,const Op op,const Transpose trans)=0;

  virtual void Start(deviceMemory<float> &buf,const int k,const Op op,const Transpose trans)=0;
  virtual void Start(deviceMemory<double> &buf,const int k,const Op op,const Transpose trans)=0;
  virtual void Start(deviceMemory<int> &buf,const int k,const Op op,const Transpose trans)=0;
  virtual void Start(deviceMemory<long long int> &buf,const int k,const Op op,const Transpose trans)=0;
  virtual void Finish(deviceMemory<float> &buf,const int k,const Op op,const Transpose trans)=0;
  virtual void Finish(deviceMemory<double> &buf,const int k,const Op op,const Transpose trans)=0;
  virtual void Finish(deviceMemory<int> &buf,const int k,const Op op,const Transpose trans)=0;
  virtual void Finish(deviceMemory<long long int> &buf,const int k,const Op op,const Transpose trans)=0;

  virtual void AllocBuffer(size_t Nbytes)=0;

  friend void InitializeKernels(platform_t& platform, const Type type, const Op op);
};

//MPI communcation via single MPI_Alltoallv call
class ogsAllToAll_t: public ogsExchange_t {
private:

  dlong NsendN=0, NsendT=0;
  memory<dlong> sendIdsN, sendIdsT;
  deviceMemory<dlong> o_sendIdsN, o_sendIdsT;

  ogsOperator_t postmpi;

  memory<int> mpiSendCountsN;
  memory<int> mpiSendCountsT;
  memory<int> mpiRecvCountsN;
  memory<int> mpiRecvCountsT;
  memory<int> mpiSendOffsetsN;
  memory<int> mpiSendOffsetsT;
  memory<int> mpiRecvOffsetsN;
  memory<int> mpiRecvOffsetsT;

  memory<int> sendCounts;
  memory<int> recvCounts;
  memory<int> sendOffsets;
  memory<int> recvOffsets;

  comm_t::request_t request;

public:
  ogsAllToAll_t(dlong Nshared,
               memory<parallelNode_t> &sharedNodes,
               ogsOperator_t &gatherHalo,
               stream_t _dataStream,
               comm_t _comm,
               platform_t &_platform);

  template<typename T>
  void Start(pinnedMemory<T> &buf,
                const int k,
                const Op op,
                const Transpose trans);

  template<typename T>
  void Finish(pinnedMemory<T> &buf,
                const int k,
                const Op op,
                const Transpose trans);

  virtual void Start(pinnedMemory<float> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(pinnedMemory<double> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(pinnedMemory<int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(pinnedMemory<long long int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(pinnedMemory<float> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(pinnedMemory<double> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(pinnedMemory<int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(pinnedMemory<long long int> &buf,const int k,const Op op,const Transpose trans);

  template<typename T>
  void Start(deviceMemory<T> &buf,
                const int k,
                const Op op,
                const Transpose trans);

  template<typename T>
  void Finish(deviceMemory<T> &buf,
                const int k,
                const Op op,
                const Transpose trans);

  virtual void Start(deviceMemory<float> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(deviceMemory<double> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(deviceMemory<int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(deviceMemory<long long int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(deviceMemory<float> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(deviceMemory<double> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(deviceMemory<int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(deviceMemory<long long int> &buf,const int k,const Op op,const Transpose trans);

  virtual void AllocBuffer(size_t Nbytes);

};

//MPI communcation via pairwise send/recvs
class ogsPairwise_t: public ogsExchange_t {
private:

  dlong NsendN=0, NsendT=0;
  memory<dlong> sendIdsN, sendIdsT;
  deviceMemory<dlong> o_sendIdsN, o_sendIdsT;

  ogsOperator_t postmpi;

  int NranksSendN=0, NranksRecvN=0;
  int NranksSendT=0, NranksRecvT=0;
  memory<int> sendRanksN;
  memory<int> sendRanksT;
  memory<int> recvRanksN;
  memory<int> recvRanksT;
  memory<int> sendCountsN;
  memory<int> sendCountsT;
  memory<int> recvCountsN;
  memory<int> recvCountsT;
  memory<int> sendOffsetsN;
  memory<int> sendOffsetsT;
  memory<int> recvOffsetsN;
  memory<int> recvOffsetsT;
  memory<comm_t::request_t> requests;

public:
  ogsPairwise_t(dlong Nshared,
               memory<parallelNode_t> &sharedNodes,
               ogsOperator_t &gatherHalo,
               stream_t _dataStream,
               comm_t _comm,
               platform_t &_platform);

  template<typename T>
  void Start(pinnedMemory<T> &buf,
                const int k,
                const Op op,
                const Transpose trans);

  template<typename T>
  void Finish(pinnedMemory<T> &buf,
                const int k,
                const Op op,
                const Transpose trans);

  virtual void Start(pinnedMemory<float> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(pinnedMemory<double> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(pinnedMemory<int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(pinnedMemory<long long int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(pinnedMemory<float> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(pinnedMemory<double> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(pinnedMemory<int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(pinnedMemory<long long int> &buf,const int k,const Op op,const Transpose trans);

  template<typename T>
  void Start(deviceMemory<T> &buf,
                const int k,
                const Op op,
                const Transpose trans);

  template<typename T>
  void Finish(deviceMemory<T> &buf,
                const int k,
                const Op op,
                const Transpose trans);

  virtual void Start(deviceMemory<float> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(deviceMemory<double> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(deviceMemory<int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(deviceMemory<long long int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(deviceMemory<float> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(deviceMemory<double> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(deviceMemory<int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(deviceMemory<long long int> &buf,const int k,const Op op,const Transpose trans);

  virtual void AllocBuffer(size_t Nbytes);
};

//MPI communcation via Crystal Router
class ogsCrystalRouter_t: public ogsExchange_t {
private:

  struct crLevel {
    int Nmsg;
    int partner;

    int Nsend, Nrecv0, Nrecv1;
    dlong recvOffset;

    memory<dlong> sendIds;
    deviceMemory<dlong> o_sendIds;

    ogsOperator_t gather;
  };

  int buf_id=0, hbuf_id=0;
  pinnedMemory<char> h_work[2];
  deviceMemory<char> o_work[2];

  memory<comm_t::request_t> request;

  int Nlevels=0;
  memory<crLevel> levelsN;
  memory<crLevel> levelsT;

  int NsendMax=0, NrecvMax=0;

public:
  ogsCrystalRouter_t(dlong Nshared,
                   memory<parallelNode_t> &sharedNodes,
                   ogsOperator_t &gatherHalo,
                   stream_t _dataStream,
                   comm_t _comm,
                   platform_t &_platform);

  template<typename T>
  void Start(pinnedMemory<T> &buf,
                const int k,
                const Op op,
                const Transpose trans);

  template<typename T>
  void Finish(pinnedMemory<T> &buf,
                const int k,
                const Op op,
                const Transpose trans);

  virtual void Start(pinnedMemory<float> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(pinnedMemory<double> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(pinnedMemory<int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(pinnedMemory<long long int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(pinnedMemory<float> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(pinnedMemory<double> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(pinnedMemory<int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(pinnedMemory<long long int> &buf,const int k,const Op op,const Transpose trans);

  template<typename T>
  void Start(deviceMemory<T> &buf,
                const int k,
                const Op op,
                const Transpose trans);

  template<typename T>
  void Finish(deviceMemory<T> &buf,
                const int k,
                const Op op,
                const Transpose trans);

  virtual void Start(deviceMemory<float> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(deviceMemory<double> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(deviceMemory<int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Start(deviceMemory<long long int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(deviceMemory<float> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(deviceMemory<double> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(deviceMemory<int> &buf,const int k,const Op op,const Transpose trans);
  virtual void Finish(deviceMemory<long long int> &buf,const int k,const Op op,const Transpose trans);

  virtual void AllocBuffer(size_t Nbytes);
};

} //namespace ogs

} //namespace libp

#endif
