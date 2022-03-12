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

  pinnedMemory<dfloat>   buf = exchange->h_workspace;
  deviceMemory<dfloat> o_buf = exchange->o_workspace;

  device_t &device = exchange->platform.device;

  //dry run
  for (int n=0;n<Ncold;++n) {
    if (exchange->gpu_aware) {
      /*GPU-aware exchange*/
      exchange->Start (o_buf, 1, Add, Sym);
      exchange->Finish(o_buf, 1, Add, Sym);
    } else {
      //if not using gpu-aware mpi move the halo buffer to the host
      o_buf.copyTo(buf, exchange->Nhalo,
                   0, "async: true");
      device.finish();

      /*MPI exchange of host buffer*/
      exchange->Start (buf, 1, Add, Sym);
      exchange->Finish(buf, 1, Add, Sym);

      // copy recv back to device
      o_buf.copyFrom(buf, exchange->Nhalo,
                     0, "async: true");
      device.finish(); //wait for transfer to finish
    }
  }

  //hot runs
  timePoint_t start = Time();
  for (int n=0;n<Nhot;++n) {
    if (exchange->gpu_aware) {
      /*GPU-aware exchange*/
      exchange->Start (o_buf, 1, Add, Sym);
      exchange->Finish(o_buf, 1, Add, Sym);
    } else {
      //if not using gpu-aware mpi move the halo buffer to the host
      o_buf.copyTo(buf, exchange->Nhalo,
                   0, "async: true");
      device.finish();

      /*MPI exchange of host buffer*/
      exchange->Start (buf, 1, Add, Sym);
      exchange->Finish(buf, 1, Add, Sym);

      // copy recv back to device
      o_buf.copyFrom(buf, exchange->Nhalo,
                     0, "async: true");
      device.finish(); //wait for transfer to finish
    }
  }
  timePoint_t end = Time();

  localTime = ElapsedTime(start,end)/Nhot;
  comm.Allreduce(localTime, sumTime, comm_t::Sum);
  comm.Allreduce(localTime, maxTime, comm_t::Max);
  comm.Allreduce(localTime, minTime, comm_t::Min);

  time[0] = sumTime/size; //avg
  time[1] = minTime;      //min
  time[2] = maxTime;      //max
}

static void HostExchangeTest(ogsExchange_t* exchange, double time[3]) {
  const int Ncold = 10;
  const int Nhot  = 10;
  double localTime, sumTime, minTime, maxTime;

  comm_t& comm = exchange->comm;
  int size = comm.size();

  pinnedMemory<dfloat> buf = exchange->h_workspace;

  //dry run
  for (int n=0;n<Ncold;++n) {
    exchange->Start (buf, 1, Add, Sym);
    exchange->Finish(buf, 1, Add, Sym);
  }

  //hot runs
  timePoint_t start = Time();
  for (int n=0;n<Nhot;++n) {
    exchange->Start (buf, 1, Add, Sym);
    exchange->Finish(buf, 1, Add, Sym);
  }
  timePoint_t end = Time();

  localTime = ElapsedTime(start,end)/Nhot;
  comm.Allreduce(localTime, sumTime, comm_t::Sum);
  comm.Allreduce(localTime, maxTime, comm_t::Max);
  comm.Allreduce(localTime, minTime, comm_t::Min);

  time[0] = sumTime/size; //avg
  time[1] = minTime;      //min
  time[2] = maxTime;      //max
}

ogsExchange_t* ogsBase_t::AutoSetup(dlong Nshared,
                                    memory<parallelNode_t> &sharedNodes,
                                    ogsOperator_t& _gatherHalo,
                                    comm_t _comm,
                                    platform_t &_platform,
                                    const int verbose) {

  int rank, size;
  rank = comm.rank();
  size = comm.size();

  if (size==1) return new ogsPairwise_t(Nshared, sharedNodes,
                                        _gatherHalo, dataStream,
                                        comm, platform);

  ogsExchange_t* bestExchange;
  Method method;
  double bestTime;

#ifdef GPU_AWARE_MPI
  if (rank==0 && verbose)
    printf("   Method         Device Exchange (avg, min, max)  Device Exchange (GPU-aware)      Host Exchange \n");
#else
  if (rank==0 && verbose)
    printf("   Method         Device Exchange (avg, min, max)  Host Exchange \n");
#endif

  //Trigger JIT kernel builds
  InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  /********************************
   * Pairwise
   ********************************/
  ogsExchange_t* pairwise = new ogsPairwise_t(Nshared, sharedNodes,
                                              _gatherHalo, dataStream,
                                              comm, platform);

  //standard copy to host - exchange - copy back to device
  pairwise->gpu_aware=false;

  double pairwiseTime[3];
  DeviceExchangeTest(pairwise, pairwiseTime);
  double pairwiseAvg = pairwiseTime[0];

#ifdef GPU_AWARE_MPI
  //test GPU-aware exchange
  pairwise->gpu_aware=true;

  double pairwiseGATime[3];
  DeviceExchangeTest(pairwise, pairwiseGATime);

  if (pairwiseGATime[0] < pairwiseAvg)
    pairwiseAvg = pairwiseGATime[0];
  else
    pairwise->gpu_aware=false;

#endif

  //test exchange from host memory (just for reporting)
  double pairwiseHostTime[3];
  HostExchangeTest(pairwise, pairwiseHostTime);

  bestExchange = pairwise;
  method = Pairwise;
  bestTime = pairwiseAvg;

#ifdef GPU_AWARE_MPI
  if (rank==0 && verbose)
    printf("   Pairwise       %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e \n",
            pairwiseTime[0],     pairwiseTime[1],     pairwiseTime[2],
            pairwiseGATime[0],   pairwiseGATime[1],   pairwiseGATime[2],
            pairwiseHostTime[0], pairwiseHostTime[1], pairwiseHostTime[2]);
#else
  if (rank==0 && verbose)
    printf("   Pairwise       %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e \n",
            pairwiseTime[0],     pairwiseTime[1],     pairwiseTime[2],
            pairwiseHostTime[0], pairwiseHostTime[1], pairwiseHostTime[2]);
#endif

  /********************************
   * All-to-All
   ********************************/
  ogsExchange_t* alltoall = new ogsAllToAll_t(Nshared, sharedNodes,
                                           _gatherHalo, dataStream,
                                           comm, platform);
  //standard copy to host - exchange - copy back to device
  alltoall->gpu_aware=false;

  double alltoallTime[3];
  DeviceExchangeTest(alltoall, alltoallTime);
  double alltoallAvg = alltoallTime[0];

#ifdef GPU_AWARE_MPI
  //test GPU-aware exchange
  alltoall->gpu_aware=true;

  double alltoallGATime[3];
  DeviceExchangeTest(alltoall, alltoallGATime);

  if (alltoallGATime[0] < alltoallAvg)
    alltoallAvg = alltoallGATime[0];
  else
    alltoall->gpu_aware=false;

#endif

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

#ifdef GPU_AWARE_MPI
  if (rank==0 && verbose)
    printf("   AllToAll       %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e \n",
            alltoallTime[0],     alltoallTime[1],     alltoallTime[2],
            alltoallGATime[0],   alltoallGATime[1],   alltoallGATime[2],
            alltoallHostTime[0], alltoallHostTime[1], alltoallHostTime[2]);
#else
  if (rank==0 && verbose)
    printf("   AllToAll       %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e \n",
            alltoallTime[0],     alltoallTime[1],     alltoallTime[2],
            alltoallHostTime[0], alltoallHostTime[1], alltoallHostTime[2]);
#endif

  /********************************
   * Crystal Router
   ********************************/
  ogsExchange_t* crystal = new ogsCrystalRouter_t(Nshared, sharedNodes,
                                                 _gatherHalo, dataStream,
                                                 comm, platform);

  //standard copy to host - exchange - copy back to device
  crystal->gpu_aware=false;

  double crystalTime[3];
  DeviceExchangeTest(crystal, crystalTime);
  double crystalAvg = crystalTime[0];

#ifdef GPU_AWARE_MPI
  //test GPU-aware exchange
  crystal->gpu_aware=true;

  double crystalGATime[3];
  DeviceExchangeTest(crystal, crystalGATime);

  if (crystalGATime[0] < crystalAvg)
    crystalAvg = crystalGATime[0];
  else
    crystal->gpu_aware=false;

#endif

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

#ifdef GPU_AWARE_MPI
  if (rank==0 && verbose)
    printf("   CrystalRouter  %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e \n",
            crystalTime[0],     crystalTime[1],     crystalTime[2],
            crystalGATime[0],   crystalGATime[1],   crystalGATime[2],
            crystalHostTime[0], crystalHostTime[1], crystalHostTime[2]);
#else
  if (rank==0 && verbose)
    printf("   CrystalRouter  %5.3e %5.3e %5.3e    %5.3e %5.3e %5.3e \n",
            crystalTime[0],     crystalTime[1],     crystalTime[2],
            crystalHostTime[0], crystalHostTime[1], crystalHostTime[2]);
#endif

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
