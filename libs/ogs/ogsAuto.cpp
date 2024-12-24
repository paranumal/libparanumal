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

namespace libp {

namespace ogs {

static void DeviceExchangeTest(ogsExchange_t* exchange, double time[3]) {
  const int Ncold = 10;
  const int Nhot  = 10;
  double localTime, sumTime, minTime, maxTime;

  comm_t& comm = exchange->comm;
  int size = comm.size();


  device_t &device = exchange->platform.device;

  //dry run
  for (int n=0;n<Ncold;++n) {
    if (exchange->gpu_aware) {
      /*GPU-aware exchange*/
      exchange->DeviceStart (Dfloat, 1, Add, Sym);
      exchange->DeviceFinish(Dfloat, 1, Add, Sym);
    } else {
      pinnedMemory<dfloat>   sendBuf = exchange->getHostSendBuffer();
      deviceMemory<dfloat> o_sendBuf = exchange->getDeviceSendBuffer();
      //if not using gpu-aware mpi move the halo buffer to the host
      sendBuf.copyFrom(o_sendBuf, exchange->Nhalo,
                       0, properties_t("async", true));
      device.finish();

      /*MPI exchange of host buffer*/
      exchange->HostStart (Dfloat, 1, Add, Sym);
      exchange->HostFinish(Dfloat, 1, Add, Sym);

      pinnedMemory<dfloat>   recvBuf = exchange->getHostRecvBuffer();
      deviceMemory<dfloat> o_recvBuf = exchange->getDeviceRecvBuffer();

      // copy recv back to device
      recvBuf.copyTo(o_recvBuf, exchange->Nhalo,
                     0, properties_t("async", true));
      device.finish(); //wait for transfer to finish
    }
  }

  //hot runs
  timePoint_t start = Time();
  for (int n=0;n<Nhot;++n) {
    if (exchange->gpu_aware) {
      /*GPU-aware exchange*/
      exchange->DeviceStart (Dfloat, 1, Add, Sym);
      exchange->DeviceFinish(Dfloat, 1, Add, Sym);
    } else {
      pinnedMemory<dfloat>   sendBuf = exchange->getHostSendBuffer();
      deviceMemory<dfloat> o_sendBuf = exchange->getDeviceSendBuffer();
      //if not using gpu-aware mpi move the halo buffer to the host
      sendBuf.copyFrom(o_sendBuf, exchange->Nhalo,
                       0, properties_t("async", true));
      device.finish();

      /*MPI exchange of host buffer*/
      exchange->HostStart (Dfloat, 1, Add, Sym);
      exchange->HostFinish(Dfloat, 1, Add, Sym);

      pinnedMemory<dfloat>   recvBuf = exchange->getHostRecvBuffer();
      deviceMemory<dfloat> o_recvBuf = exchange->getDeviceRecvBuffer();

      // copy recv back to device
      recvBuf.copyTo(o_recvBuf, exchange->Nhalo,
                     0, properties_t("async", true));
      device.finish(); //wait for transfer to finish
    }
  }
  timePoint_t end = Time();

  localTime = ElapsedTime(start,end)/Nhot;
  comm.Allreduce(localTime, sumTime, comm_t::Sum);
  comm.Allreduce(localTime, maxTime, comm_t::Max);
  comm.Allreduce(localTime, minTime, comm_t::Min);

  time[0] = minTime;      //min
  time[1] = sumTime/size; //avg
  time[2] = maxTime;      //max
}

static void HostExchangeTest(ogsExchange_t* exchange, double time[3]) {
  const int Ncold = 10;
  const int Nhot  = 10;
  double localTime, sumTime, minTime, maxTime;

  comm_t& comm = exchange->comm;
  int size = comm.size();

  //dry run
  for (int n=0;n<Ncold;++n) {
    exchange->HostStart (Dfloat, 1, Add, Sym);
    exchange->HostFinish(Dfloat, 1, Add, Sym);
  }

  //hot runs
  timePoint_t start = Time();
  for (int n=0;n<Nhot;++n) {
    exchange->HostStart (Dfloat, 1, Add, Sym);
    exchange->HostFinish(Dfloat, 1, Add, Sym);
  }
  timePoint_t end = Time();

  localTime = ElapsedTime(start,end)/Nhot;
  comm.Allreduce(localTime, sumTime, comm_t::Sum);
  comm.Allreduce(localTime, maxTime, comm_t::Max);
  comm.Allreduce(localTime, minTime, comm_t::Min);

  time[0] = minTime;      //min
  time[1] = sumTime/size; //avg
  time[2] = maxTime;      //max
}

ogsExchange_t* ogsBase_t::AutoSetup(const dlong Nshared,
                                    const memory<int>   sharedRemoteRanks,
                                    const memory<dlong> sharedLocalRows,
                                    const memory<dlong> sharedRemoteRows,
                                    const memory<hlong> sharedLocalBaseIds,
                                    const memory<hlong> sharedRemoteBaseIds,
                                    const memory<hlong> haloBaseIds,
                                    ogsOperator_t& _gatherHalo,
                                    comm_t _comm,
                                    platform_t &_platform,
                                    const int verbose) {

  int rank, size;
  rank = comm.rank();
  size = comm.size();

  Kind knd = (kind == Unsigned) ? Unsigned : Signed;

  if (size==1) return new ogsPairwise_t(knd,
                                        Nshared,
                                        sharedRemoteRanks,
                                        sharedLocalRows,
                                        sharedRemoteRows,
                                        sharedLocalBaseIds,
                                        sharedRemoteBaseIds,
                                        _gatherHalo,
                                        dataStream,
                                        comm,
                                        platform);

  ogsExchange_t* bestExchange;
  Method method;
  double bestTime;

  if (comm.gpuAware()) {
    if (rank==0 && verbose)
      printf("   Method         Device Exchange (min, avg, max)  Device Exchange (GPU-aware)      Host Exchange \n");
  } else {
    if (rank==0 && verbose)
      printf("   Method         Device Exchange (min, avg, max)  Host Exchange \n");
  }

  //Trigger JIT kernel builds
  InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  /********************************
   * Pairwise
   ********************************/
  ogsExchange_t* pairwise = new ogsPairwise_t(knd,
                                              Nshared,
                                              sharedRemoteRanks,
                                              sharedLocalRows,
                                              sharedRemoteRows,
                                              sharedLocalBaseIds,
                                              sharedRemoteBaseIds,
                                              _gatherHalo,
                                              dataStream,
                                              comm,
                                              platform);

  //standard copy to host - exchange - copy back to device
  pairwise->gpu_aware=false;

  double pairwiseTime[3];
  double pairwiseGATime[3];

  DeviceExchangeTest(pairwise, pairwiseTime);
  double pairwiseAvg = pairwiseTime[1];

  if (comm.gpuAware()) {
    //test GPU-aware exchange
    pairwise->gpu_aware=true;

    DeviceExchangeTest(pairwise, pairwiseGATime);

    if (pairwiseGATime[1] < pairwiseAvg)
      pairwiseAvg = pairwiseGATime[1];
    else
      pairwise->gpu_aware=false;
  }

  //test exchange from host memory (just for reporting)
  double pairwiseHostTime[3];
  HostExchangeTest(pairwise, pairwiseHostTime);

  bestExchange = pairwise;
  method = Pairwise;
  bestTime = pairwiseAvg;

  if (comm.gpuAware()) {
    if (rank==0 && verbose)
      printf("   Pairwise       %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e \n",
              pairwiseTime[0],     pairwiseTime[1],     pairwiseTime[2],
              pairwiseGATime[0],   pairwiseGATime[1],   pairwiseGATime[2],
              pairwiseHostTime[0], pairwiseHostTime[1], pairwiseHostTime[2]);
  } else {
    if (rank==0 && verbose)
      printf("   Pairwise       %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e \n",
              pairwiseTime[0],     pairwiseTime[1],     pairwiseTime[2],
              pairwiseHostTime[0], pairwiseHostTime[1], pairwiseHostTime[2]);
  }

  /********************************
   * All-to-All
   ********************************/
  ogsExchange_t* alltoall = new ogsAllToAll_t(knd,
                                              Nshared,
                                              sharedRemoteRanks,
                                              sharedLocalRows,
                                              sharedRemoteRows,
                                              sharedLocalBaseIds,
                                              sharedRemoteBaseIds,
                                              _gatherHalo,
                                              dataStream,
                                              comm,
                                              platform);
  //standard copy to host - exchange - copy back to device
  alltoall->gpu_aware=false;

  double alltoallTime[3];
  double alltoallGATime[3];

  DeviceExchangeTest(alltoall, alltoallTime);
  double alltoallAvg = alltoallTime[0];

  if (comm.gpuAware()) {
    //test GPU-aware exchange
    alltoall->gpu_aware=true;

    DeviceExchangeTest(alltoall, alltoallGATime);

    if (alltoallGATime[1] < alltoallAvg)
      alltoallAvg = alltoallGATime[1];
    else
      alltoall->gpu_aware=false;
  }

  //test exchange from host memory (just for reporting)
  double alltoallHostTime[3];
  HostExchangeTest(alltoall, alltoallHostTime);

  if (alltoallAvg < bestTime) {
    delete bestExchange;
    bestExchange = alltoall;
    method = AllToAll;
    bestTime = alltoallAvg;
  } else {
    delete alltoall;
  }

  if (comm.gpuAware()) {
    if (rank==0 && verbose)
      printf("   AllToAll       %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e \n",
              alltoallTime[0],     alltoallTime[1],     alltoallTime[2],
              alltoallGATime[0],   alltoallGATime[1],   alltoallGATime[2],
              alltoallHostTime[0], alltoallHostTime[1], alltoallHostTime[2]);
  } else {
    if (rank==0 && verbose)
      printf("   AllToAll       %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e \n",
              alltoallTime[0],     alltoallTime[1],     alltoallTime[2],
              alltoallHostTime[0], alltoallHostTime[1], alltoallHostTime[2]);
  }

  /********************************
   * Crystal Router
   ********************************/
  ogsExchange_t* crystal = new ogsCrystalRouter_t(knd,
                                                  Nshared,
                                                  sharedRemoteRanks,
                                                  sharedLocalRows,
                                                  sharedRemoteRows,
                                                  sharedLocalBaseIds,
                                                  sharedRemoteBaseIds,
                                                  haloBaseIds,
                                                  _gatherHalo,
                                                  dataStream,
                                                  comm,
                                                  platform);

  //standard copy to host - exchange - copy back to device
  crystal->gpu_aware=false;

  double crystalTime[3];
  double crystalGATime[3];

  DeviceExchangeTest(crystal, crystalTime);
  double crystalAvg = crystalTime[0];

  if (comm.gpuAware()) {
    //test GPU-aware exchange
    crystal->gpu_aware=true;

    DeviceExchangeTest(crystal, crystalGATime);

    if (crystalGATime[1] < crystalAvg)
      crystalAvg = crystalGATime[1];
    else
      crystal->gpu_aware=false;
  }

  //test exchange from host memory (just for reporting)
  double crystalHostTime[3];
  HostExchangeTest(crystal, crystalHostTime);

  if (crystalAvg < bestTime) {
    delete bestExchange;
    bestExchange = crystal;
    method = CrystalRouter;
    bestTime = crystalAvg;
  } else {
    delete crystal;
  }

  if (comm.gpuAware()) {
    if (rank==0 && verbose)
      printf("   CrystalRouter  %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e \n",
              crystalTime[0],     crystalTime[1],     crystalTime[2],
              crystalGATime[0],   crystalGATime[1],   crystalGATime[2],
              crystalHostTime[0], crystalHostTime[1], crystalHostTime[2]);
  } else {
    if (rank==0 && verbose)
      printf("   CrystalRouter  %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e \n",
              crystalTime[0],     crystalTime[1],     crystalTime[2],
              crystalHostTime[0], crystalHostTime[1], crystalHostTime[2]);
  }

  if (rank==0 && verbose) {
    switch (method) {
      case AllToAll:
        printf("   Exchange method selected: AllToAll"); break;
      case Pairwise:
        printf("   Exchange method selected: Pairwise"); break;
      case CrystalRouter:
        printf("   Exchange method selected: CrystalRouter"); break;
      default:
        break;
    }
    if (bestExchange->gpu_aware) printf(" (GPU-aware)");
    printf("\n");
  }

  return bestExchange;
}


} //namespace ogs

} //namespace libp
