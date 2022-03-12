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

namespace libp {

namespace ogs {

/********************************
 * Device GatherScatter
 ********************************/
template<typename T>
void ogs_t::GatherScatter(deviceMemory<T> o_v,
                          const int k,
                          const Op op,
                          const Transpose trans){
  GatherScatterStart (o_v, k, op, trans);
  GatherScatterFinish(o_v, k, op, trans);
}

template<typename T>
void ogs_t::GatherScatterStart(deviceMemory<T> o_v,
                               const int k,
                               const Op op,
                               const Transpose trans){
  exchange->AllocBuffer(k*sizeof(T));

  deviceMemory<T> o_haloBuf = exchange->o_workspace;

  //collect halo buffer
  gatherHalo->Gather(o_haloBuf, o_v, k, op, trans);

  if (exchange->gpu_aware) {
    //prepare MPI exchange
    exchange->Start(o_haloBuf, k, op, trans);
  } else {
    //get current stream
    device_t &device = platform.device;
    stream_t currentStream = device.getStream();

    pinnedMemory<T> haloBuf = exchange->h_workspace;

    //if not using gpu-aware mpi move the halo buffer to the host
    const dlong Nhalo = (trans == NoTrans) ? NhaloP : NhaloT;

    //wait for o_haloBuf to be ready
    device.finish();

    //queue copy to host
    device.setStream(dataStream);
    haloBuf.copyFrom(o_haloBuf, Nhalo*k,
                     0, "async: true");
    device.setStream(currentStream);
  }
}

template<typename T>
void ogs_t::GatherScatterFinish(deviceMemory<T> o_v,
                                const int k,
                                const Op op,
                                const Transpose trans){

  //queue local gs operation
  gatherLocal->GatherScatter(o_v, k, op, trans);

  deviceMemory<T> o_haloBuf = exchange->o_workspace;

  if (exchange->gpu_aware) {
    //finish MPI exchange
    exchange->Finish(o_haloBuf, k, op, trans);
  } else {
    pinnedMemory<T> haloBuf = exchange->h_workspace;

    //get current stream
    device_t &device = platform.device;
    stream_t currentStream = device.getStream();

    //synchronize data stream to ensure the buffer is on the host
    device.setStream(dataStream);
    device.finish();

    /*MPI exchange of host buffer*/
    exchange->Start (haloBuf, k, op, trans);
    exchange->Finish(haloBuf, k, op, trans);

    // copy recv back to device
    const dlong Nhalo = (trans == Trans) ? NhaloP : NhaloT;
    haloBuf.copyTo(o_haloBuf, Nhalo*k,
                   0, "async: true");
    device.finish(); //wait for transfer to finish
    device.setStream(currentStream);
  }

  //write exchanged halo buffer back to vector
  gatherHalo->Scatter(o_v, o_haloBuf, k, trans);
}

template
void ogs_t::GatherScatter(deviceMemory<float> v, const int k,
                          const Op op, const Transpose trans);
template
void ogs_t::GatherScatter(deviceMemory<double> v, const int k,
                          const Op op, const Transpose trans);
template
void ogs_t::GatherScatter(deviceMemory<int> v, const int k,
                          const Op op, const Transpose trans);
template
void ogs_t::GatherScatter(deviceMemory<long long int> v, const int k,
                          const Op op, const Transpose trans);

/********************************
 * Host GatherScatter
 ********************************/
template<typename T>
void ogs_t::GatherScatter(memory<T> v,
                          const int k,
                          const Op op,
                          const Transpose trans){
  GatherScatterStart (v, k, op, trans);
  GatherScatterFinish(v, k, op, trans);
}

template<typename T>
void ogs_t::GatherScatterStart(memory<T> v,
                               const int k,
                               const Op op,
                               const Transpose trans){
  exchange->AllocBuffer(k*sizeof(T));

  /*Cast workspace to type T*/
  pinnedMemory<T> haloBuf = exchange->h_workspace;

  //collect halo buffer
  gatherHalo->Gather(haloBuf, v, k, op, trans);

  //prepare MPI exchange
  exchange->Start(haloBuf, k, op, trans);
}

template<typename T>
void ogs_t::GatherScatterFinish(memory<T> v,
                                const int k,
                                const Op op,
                                const Transpose trans){

  /*Cast workspace to type T*/
  pinnedMemory<T> haloBuf = exchange->h_workspace;

  //queue local gs operation
  gatherLocal->GatherScatter(v, k, op, trans);

  //finish MPI exchange
  exchange->Finish(haloBuf, k, op, trans);

  //write exchanged halo buffer back to vector
  gatherHalo->Scatter(v, haloBuf, k, trans);
}

template
void ogs_t::GatherScatter(memory<float> v, const int k,
                          const Op op, const Transpose trans);
template
void ogs_t::GatherScatter(memory<double> v, const int k,
                          const Op op, const Transpose trans);
template
void ogs_t::GatherScatter(memory<int> v, const int k,
                          const Op op, const Transpose trans);
template
void ogs_t::GatherScatter(memory<long long int> v, const int k,
                          const Op op, const Transpose trans);

/********************************
 * Device Gather
 ********************************/
template<typename T>
void ogs_t::Gather(deviceMemory<T> o_gv,
                   deviceMemory<T> o_v,
                   const int k,
                   const Op op,
                   const Transpose trans){
  GatherStart (o_gv, o_v, k, op, trans);
  GatherFinish(o_gv, o_v, k, op, trans);
}

template<typename T>
void ogs_t::GatherStart(deviceMemory<T> o_gv,
                        deviceMemory<T> o_v,
                        const int k,
                        const Op op,
                        const Transpose trans){
  AssertGatherDefined();

  deviceMemory<T> o_haloBuf = exchange->o_workspace;

  if (trans==Trans) { //if trans!=ogs::Trans theres no comms required
    exchange->AllocBuffer(k*sizeof(T));

    //collect halo buffer
    gatherHalo->Gather(o_haloBuf, o_v, k, op, Trans);

    if (exchange->gpu_aware) {
      //prepare MPI exchange
      exchange->Start(o_haloBuf, k, op, Trans);
    } else {
      //get current stream
      device_t &device = platform.device;
      stream_t currentStream = device.getStream();

      //if not using gpu-aware mpi move the halo buffer to the host
      pinnedMemory<T> haloBuf = exchange->h_workspace;

      //wait for o_haloBuf to be ready
      device.finish();

      //queue copy to host
      device.setStream(dataStream);
      haloBuf.copyFrom(o_haloBuf, NhaloT*k,
                       0, "async: true");
      device.setStream(currentStream);
    }
  } else {
    //gather halo
    gatherHalo->Gather(o_gv + k*NlocalT, o_v, k, op, trans);
  }
}

template<typename T>
void ogs_t::GatherFinish(deviceMemory<T> o_gv,
                         deviceMemory<T> o_v,
                         const int k,
                         const Op op,
                         const Transpose trans){
  AssertGatherDefined();

  deviceMemory<T> o_haloBuf = exchange->o_workspace;

  //queue local g operation
  gatherLocal->Gather(o_gv, o_v, k, op, trans);

  if (trans==Trans) { //if trans!=ogs::Trans theres no comms required
    if (exchange->gpu_aware) {
      //finish MPI exchange
      exchange->Finish(o_haloBuf, k, op, Trans);

      //put the result at the end of o_gv
      o_haloBuf.copyTo(o_gv + k*NlocalT,
                       k*NhaloP, 0, "async: true");
    } else {
      pinnedMemory<T> haloBuf = exchange->h_workspace;

      //get current stream
      device_t &device = platform.device;
      stream_t currentStream = device.getStream();

      //synchronize data stream to ensure the buffer is on the host
      device.setStream(dataStream);
      device.finish();

      /*MPI exchange of host buffer*/
      exchange->Start (haloBuf, k, op, trans);
      exchange->Finish(haloBuf, k, op, trans);

      // copy recv back to device
      //put the result at the end of o_gv
      haloBuf.copyTo(o_gv + k*NlocalT, k*NhaloP,
                     0, "async: true");
      device.finish(); //wait for transfer to finish
      device.setStream(currentStream);
    }
  }
}

template
void ogs_t::Gather(deviceMemory<float> v, const deviceMemory<float> gv,
                   const int k, const Op op, const Transpose trans);
template
void ogs_t::Gather(deviceMemory<double> v, const deviceMemory<double> gv,
                   const int k, const Op op, const Transpose trans);
template
void ogs_t::Gather(deviceMemory<int> v, const deviceMemory<int> gv,
                   const int k, const Op op, const Transpose trans);
template
void ogs_t::Gather(deviceMemory<long long int> v, const deviceMemory<long long int> gv,
                   const int k, const Op op, const Transpose trans);

/********************************
 * Host Gather
 ********************************/

//host versions
template<typename T>
void ogs_t::Gather(memory<T> gv,
                   const memory<T> v,
                   const int k,
                   const Op op,
                   const Transpose trans){
  GatherStart (gv, v, k, op, trans);
  GatherFinish(gv, v, k, op, trans);
}

template<typename T>
void ogs_t::GatherStart(memory<T> gv,
                        const memory<T> v,
                        const int k,
                        const Op op,
                        const Transpose trans){
  AssertGatherDefined();

  if (trans==Trans) { //if trans!=ogs::Trans theres no comms required
    exchange->AllocBuffer(k*sizeof(T));

    /*Cast workspace to type T*/
    pinnedMemory<T> haloBuf = exchange->h_workspace;

    //collect halo buffer
    gatherHalo->Gather(haloBuf, v, k, op, Trans);

    //prepare MPI exchange
    exchange->Start(haloBuf, k, op, Trans);
  } else {
    //gather halo
    gatherHalo->Gather(gv + k*NlocalT, v, k, op, trans);
  }
}

template<typename T>
void ogs_t::GatherFinish(memory<T> gv,
                         const memory<T> v,
                         const int k,
                         const Op op,
                         const Transpose trans){
  AssertGatherDefined();

  //queue local g operation
  gatherLocal->Gather(gv, v, k, op, trans);

  if (trans==Trans) { //if trans!=ogs::Trans theres no comms required
    /*Cast workspace to type T*/
    pinnedMemory<T> haloBuf = exchange->h_workspace;

    //finish MPI exchange
    exchange->Finish(haloBuf, k, op, Trans);

    //put the result at the end of o_gv
    haloBuf.copyTo(gv+k*NlocalT, k*NhaloP);
  }
}

template
void ogs_t::Gather(memory<float> v, const memory<float> gv,
                   const int k, const Op op, const Transpose trans);
template
void ogs_t::Gather(memory<double> v, const memory<double> gv,
                   const int k, const Op op, const Transpose trans);
template
void ogs_t::Gather(memory<int> v, const memory<int> gv,
                   const int k, const Op op, const Transpose trans);
template
void ogs_t::Gather(memory<long long int> v, const memory<long long int> gv,
                   const int k, const Op op, const Transpose trans);

/********************************
 * Device Scatter
 ********************************/
template<typename T>
void ogs_t::Scatter(deviceMemory<T> o_v,
                    deviceMemory<T> o_gv,
                    const int k,
                    const Transpose trans){
  ScatterStart (o_v, o_gv, k, trans);
  ScatterFinish(o_v, o_gv, k, trans);
}

template<typename T>
void ogs_t::ScatterStart(deviceMemory<T> o_v,
                         deviceMemory<T> o_gv,
                         const int k,
                         const Transpose trans){
  AssertGatherDefined();

  deviceMemory<T> o_haloBuf = exchange->o_workspace;

  if (trans==NoTrans) { //if trans!=ogs::NoTrans theres no comms required
    exchange->AllocBuffer(k*sizeof(T));

    device_t &device = platform.device;

    if (exchange->gpu_aware) {
      //collect halo buffer
      o_haloBuf.copyFrom(o_gv + k*NlocalT,
                         k*NhaloP, 0, "async: true");

      //wait for o_haloBuf to be ready
      device.finish();

      //prepare MPI exchange
      exchange->Start(o_haloBuf, k, Add, NoTrans);
    } else {
      //get current stream
      stream_t currentStream = device.getStream();

      //if not using gpu-aware mpi move the halo buffer to the host
      pinnedMemory<T> haloBuf = exchange->h_workspace;

      //wait for o_gv to be ready
      device.finish();

      //queue copy to host
      device.setStream(dataStream);
      haloBuf.copyFrom(o_gv + k*NlocalT, NhaloP*k,
                       0, "async: true");
      device.setStream(currentStream);
    }
  }
}

template<typename T>
void ogs_t::ScatterFinish(deviceMemory<T> o_v,
                          deviceMemory<T> o_gv,
                          const int k,
                          const Transpose trans){
  AssertGatherDefined();

  deviceMemory<T> o_haloBuf = exchange->o_workspace;

  //queue local s operation
  gatherLocal->Scatter(o_v, o_gv, k, trans);

  if (trans==NoTrans) { //if trans!=ogs::NoTrans theres no comms required
    if (exchange->gpu_aware) {
      //finish MPI exchange
      exchange->Finish(o_haloBuf, k, Add, NoTrans);
    } else {
      pinnedMemory<T> haloBuf = exchange->h_workspace;

      //get current stream
      device_t &device = platform.device;
      stream_t currentStream = device.getStream();

      //synchronize data stream to ensure the buffer is on the host
      device.setStream(dataStream);
      device.finish();

      /*MPI exchange of host buffer*/
      exchange->Start (haloBuf, k, Add, NoTrans);
      exchange->Finish(haloBuf, k, Add, NoTrans);

      // copy recv back to device
      haloBuf.copyTo(o_haloBuf, NhaloT*k,
                     0, "async: true");
      device.finish(); //wait for transfer to finish
      device.setStream(currentStream);
    }

    //scatter halo buffer
    gatherHalo->Scatter(o_v, o_haloBuf, k, NoTrans);
  } else {
    //scatter halo
    gatherHalo->Scatter(o_v, o_gv + k*NlocalT, k, trans);
  }
}

template
void ogs_t::Scatter(deviceMemory<float> v, const deviceMemory<float> gv,
                    const int k, const Transpose trans);
template
void ogs_t::Scatter(deviceMemory<double> v, const deviceMemory<double> gv,
                    const int k, const Transpose trans);
template
void ogs_t::Scatter(deviceMemory<int> v, const deviceMemory<int> gv,
                    const int k, const Transpose trans);
template
void ogs_t::Scatter(deviceMemory<long long int> v, const deviceMemory<long long int> gv,
                    const int k, const Transpose trans);

/********************************
 * Host Scatter
 ********************************/

//host versions
template<typename T>
void ogs_t::Scatter(memory<T> v,
                    const memory<T> gv,
                    const int k,
                    const Transpose trans){
  ScatterStart (v, gv, k, trans);
  ScatterFinish(v, gv, k, trans);
}

template<typename T>
void ogs_t::ScatterStart(memory<T> v,
                         const memory<T> gv,
                         const int k,
                         const Transpose trans){
  AssertGatherDefined();

  if (trans==NoTrans) { //if trans!=ogs::NoTrans theres no comms required
    exchange->AllocBuffer(k*sizeof(T));

    /*Cast workspace to type T*/
    pinnedMemory<T> haloBuf = exchange->h_workspace;

    //collect halo buffer
    haloBuf.copyFrom(gv + k*NlocalT, k*NhaloP);

    //prepare MPI exchange
    exchange->Start(haloBuf, k, Add, NoTrans);
  }
}

template<typename T>
void ogs_t::ScatterFinish(memory<T> v,
                          const memory<T> gv,
                          const int k,
                          const Transpose trans){
  AssertGatherDefined();

  //queue local s operation
  gatherLocal->Scatter(v, gv, k, trans);

  if (trans==NoTrans) { //if trans!=ogs::NoTrans theres no comms required
    /*Cast workspace to type T*/
    pinnedMemory<T> haloBuf = exchange->h_workspace;

    //finish MPI exchange (and put the result at the end of o_gv)
    exchange->Finish(haloBuf, k, Add, NoTrans);

    //scatter halo buffer
    gatherHalo->Scatter(v, haloBuf, k, NoTrans);
  } else {
    //scatter halo
    gatherHalo->Scatter(v, gv + k*NlocalT, k, trans);
  }
}

template
void ogs_t::Scatter(memory<float> v, const memory<float> gv,
                    const int k, const Transpose trans);
template
void ogs_t::Scatter(memory<double> v, const memory<double> gv,
                    const int k, const Transpose trans);
template
void ogs_t::Scatter(memory<int> v, const memory<int> gv,
                    const int k, const Transpose trans);
template
void ogs_t::Scatter(memory<long long int> v, const memory<long long int> gv,
                    const int k, const Transpose trans);
} //namespace ogs

} //namespace libp
