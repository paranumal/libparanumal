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

  stream_t dataStream;
  static kernel_t extractKernel[4];

  bool gpu_aware=false;

  ogsExchange_t(platform_t &_platform, comm_t _comm,
                stream_t _datastream):
    platform(_platform),
    comm(_comm),
    dataStream(_datastream) {
    rank = comm.rank();
    size = comm.size();
    gpu_aware = comm.gpuAware();
  }
  virtual ~ogsExchange_t() {}

  pinnedMemory<char> getHostSendBuffer() { return h_sendspace; }
  pinnedMemory<char> getHostRecvBuffer() { return h_recvspace; }
  deviceMemory<char> getDeviceSendBuffer() { return o_sendspace; }
  deviceMemory<char> getDeviceRecvBuffer() { return o_recvspace; }

  virtual void HostStart(const Type type, const int k,const Op op,const Transpose trans)=0;
  virtual void HostFinish(const Type type, const int k,const Op op,const Transpose trans)=0;
  virtual void DeviceStart(const Type type, const int k,const Op op,const Transpose trans)=0;
  virtual void DeviceFinish(const Type type, const int k,const Op op,const Transpose trans)=0;

  virtual void AllocBuffer(size_t Nbytes)=0;

  friend void InitializeKernels(platform_t& platform, const Type type, const Op op);

protected:
  pinnedMemory<char> h_workspace, h_sendspace, h_recvspace;
  deviceMemory<char> o_workspace, o_sendspace, o_recvspace;
};

//MPI communcation via single MPI_Alltoallv call
class ogsAllToAll_t: public ogsExchange_t {
public:
  struct data_t {
    dlong Nsend=0;
    dlong NrowsP=0;
    dlong Nrows=0;
    dlong Ncols=0;

    memory<dlong> sendIds;
    deviceMemory<dlong> o_sendIds;

    memory<int> sendCounts;
    memory<int> recvCounts;
    memory<int> sendOffsets;
    memory<int> recvOffsets;

    ogsOperator_t postmpi;

    void setupExchange(const dlong Nsend,
                       const dlong NrowsP,
                       const dlong NrowsT,
                       const memory<int>   destRanks,
                       const memory<dlong> localRows,
                       const memory<dlong> remoteRows,
                       comm_t comm,
                       platform_t& platform);
  };

  ogsAllToAll_t(Kind kind,
                const dlong Nshared,
                const memory<int>   sharedRemoteRanks,
                const memory<dlong> sharedLocalRows,
                const memory<dlong> sharedRemoteRows,
                const memory<hlong> sharedLocalBaseIds,
                const memory<hlong> sharedRemoteBaseIds,
                ogsOperator_t &gatherHalo,
                stream_t _dataStream,
                comm_t _comm,
                platform_t &_platform);

  template<typename T>
  void HostStart(const int k,
                 const Op op,
                 const Transpose trans);

  template<typename T>
  void HostFinish(const int k,
                  const Op op,
                  const Transpose trans);

  void HostStart(const Type type, const int k,const Op op,const Transpose trans) override;
  void HostFinish(const Type type, const int k,const Op op,const Transpose trans) override;

  template<typename T>
  void DeviceStart(const int k,
                   const Op op,
                   const Transpose trans);

  template<typename T>
  void DeviceFinish(const int k,
                    const Op op,
                    const Transpose trans);

  void DeviceStart(const Type type, const int k,const Op op,const Transpose trans) override;
  void DeviceFinish(const Type type, const int k,const Op op,const Transpose trans) override;


  void AllocBuffer(size_t Nbytes) override;

private:
  data_t data[3];

  memory<int> sendCounts;
  memory<int> recvCounts;
  memory<int> sendOffsets;
  memory<int> recvOffsets;

  comm_t::request_t request;
};

//MPI communcation via pairwise send/recvs
class ogsPairwise_t: public ogsExchange_t {
public:
  struct data_t {
    dlong Nsend=0;
    dlong NrowsP=0;
    dlong Nrows=0;
    dlong Ncols=0;

    memory<dlong> sendIds;
    deviceMemory<dlong> o_sendIds;

    int NranksSend=0, NranksRecv=0;
    memory<int> sendRanks;
    memory<int> recvRanks;
    memory<int> sendCounts;
    memory<int> recvCounts;
    memory<int> sendOffsets;
    memory<int> recvOffsets;

    ogsOperator_t postmpi;

    void setupExchange(const dlong Nsend,
                       const dlong NrowsP,
                       const dlong NrowsT,
                       const memory<int>   destRanks,
                       const memory<dlong> localRows,
                       const memory<dlong> remoteRows,
                       comm_t comm,
                       platform_t& platform);
  };

  ogsPairwise_t(Kind kind,
                const dlong Nshared,
                const memory<int>   sharedRemoteRanks,
                const memory<dlong> sharedLocalRows,
                const memory<dlong> sharedRemoteRows,
                const memory<hlong> sharedLocalBaseIds,
                const memory<hlong> sharedRemoteBaseIds,
                ogsOperator_t &gatherHalo,
                stream_t _dataStream,
                comm_t _comm,
                platform_t &_platform);

  template<typename T>
  void HostStart(const int k,
                 const Op op,
                 const Transpose trans);

  template<typename T>
  void HostFinish(const int k,
                  const Op op,
                  const Transpose trans);

  void HostStart(const Type type, const int k,const Op op,const Transpose trans) override;
  void HostFinish(const Type type, const int k,const Op op,const Transpose trans) override;

  template<typename T>
  void DeviceStart(const int k,
                   const Op op,
                   const Transpose trans);

  template<typename T>
  void DeviceFinish(const int k,
                    const Op op,
                    const Transpose trans);

  void DeviceStart(const Type type, const int k,const Op op,const Transpose trans) override;
  void DeviceFinish(const Type type, const int k,const Op op,const Transpose trans) override;


  void AllocBuffer(size_t Nbytes) override;

private:
  data_t data[3];

  memory<comm_t::request_t> requests;
};

//MPI communcation via Crystal Router
class ogsCrystalRouter_t: public ogsExchange_t {

public:
  struct crLevel {
    int Nmsg;
    int partner;

    dlong Nids;
    int Nsend, Nrecv0, Nrecv1;

    memory<dlong> sendIds;
    deviceMemory<dlong> o_sendIds;

    ogsOperator_t gather;
  };

  struct data_t {
    memory<crLevel> levels;

    int NsendMax=0, NrecvMax=0;

    void setupExchange(const int Nlevels,
                       dlong NhaloP,
                       dlong NhaloT,
                       memory<hlong> haloBaseIds,
                       dlong Nshared,
                       memory<int>   remoteRanks,
                       memory<hlong> remoteBaseIds,
                       memory<dlong> localIds,
                       comm_t comm,
                       platform_t& platform);
  };

  ogsCrystalRouter_t(Kind kind,
                     const dlong Nshared,
                     const memory<int>   sharedRemoteRanks,
                     const memory<dlong> sharedLocalRows,
                     const memory<dlong> sharedRemoteRows,
                     const memory<hlong> sharedLocalBaseIds,
                     const memory<hlong> sharedRemoteBaseIds,
                     const memory<hlong> haloBaseIds,
                     ogsOperator_t& gatherHalo,
                     stream_t _dataStream,
                     comm_t _comm,
                     platform_t &_platform);

  template<typename T>
  void HostStart(const int k,
                 const Op op,
                 const Transpose trans);

  template<typename T>
  void HostFinish(const int k,
                  const Op op,
                  const Transpose trans);

  void HostStart(const Type type, const int k,const Op op,const Transpose trans) override;
  void HostFinish(const Type type, const int k,const Op op,const Transpose trans) override;

  template<typename T>
  void DeviceStart(const int k,
                   const Op op,
                   const Transpose trans);

  template<typename T>
  void DeviceFinish(const int k,
                    const Op op,
                    const Transpose trans);

  void DeviceStart(const Type type, const int k,const Op op,const Transpose trans) override;
  void DeviceFinish(const Type type, const int k,const Op op,const Transpose trans) override;


  void AllocBuffer(size_t Nbytes) override;

private:
  int Nlevels=0;

  data_t data[3];

  comm_t::request_t requests[3];
};

} //namespace ogs

} //namespace libp

#endif
