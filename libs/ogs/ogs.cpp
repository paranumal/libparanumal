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

  deviceMemory<T> o_sendBuf = exchange->getDeviceSendBuffer();

  //collect halo buffer
  gatherHalo->Gather(o_sendBuf, o_v, k, op, trans);

  if (exchange->gpu_aware) {
    //prepare MPI exchange
    exchange->DeviceStart(ogsType<T>::get(), k, op, trans);
  } else {
    //get current stream
    device_t &device = platform.device;
    stream_t currentStream = device.getStream();

    pinnedMemory<T> sendBuf = exchange->getHostSendBuffer();

    //if not using gpu-aware mpi move the halo buffer to the host
    const dlong Nhalo = (trans == NoTrans) ? NhaloP : NhaloT;

    //wait for o_sendBuf to be ready
    device.finish();

    //queue copy to host
    device.setStream(dataStream);
    sendBuf.copyFrom(o_sendBuf, Nhalo*k,
                     0, properties_t("async", true));
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

  if (exchange->gpu_aware) {
    //finish MPI exchange
    exchange->DeviceFinish(ogsType<T>::get(),k, op, trans);
  } else {

    //get current stream
    device_t &device = platform.device;
    stream_t currentStream = device.getStream();

    //synchronize data stream to ensure the buffer is on the host
    device.setStream(dataStream);
    device.finish();

    /*MPI exchange of host buffer*/
    exchange->HostStart (ogsType<T>::get(), k, op, trans);
    exchange->HostFinish(ogsType<T>::get(), k, op, trans);

    pinnedMemory<T> recvBuf = exchange->getHostRecvBuffer();
    deviceMemory<T> o_recvBuf = exchange->getDeviceRecvBuffer();

    // copy recv back to device
    const dlong Nhalo = (trans == Trans) ? NhaloP : NhaloT;
    recvBuf.copyTo(o_recvBuf, Nhalo*k,
                   0, properties_t("async", true));
    device.finish(); //wait for transfer to finish
    device.setStream(currentStream);
  }

  deviceMemory<T> o_recvBuf = exchange->getDeviceRecvBuffer();

  //write exchanged halo buffer back to vector
  gatherHalo->Scatter(o_v, o_recvBuf, k, trans);
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
  pinnedMemory<T> sendBuf = exchange->getHostSendBuffer();

  //collect halo buffer
  gatherHalo->Gather(sendBuf, v, k, op, trans);

  //prepare MPI exchange
  exchange->HostStart(ogsType<T>::get(), k, op, trans);
}

template<typename T>
void ogs_t::GatherScatterFinish(memory<T> v,
                                const int k,
                                const Op op,
                                const Transpose trans){

  /*Cast workspace to type T*/

  //queue local gs operation
  gatherLocal->GatherScatter(v, k, op, trans);

  //finish MPI exchange
  exchange->HostFinish(ogsType<T>::get(), k, op, trans);

  pinnedMemory<T> recvBuf = exchange->getHostRecvBuffer();

  //write exchanged halo buffer back to vector
  gatherHalo->Scatter(v, recvBuf, k, trans);
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

  if (trans==Trans) { //if trans!=ogs::Trans theres no comms required
    exchange->AllocBuffer(k*sizeof(T));

    deviceMemory<T> o_sendBuf = exchange->getDeviceSendBuffer();

    //collect halo buffer
    gatherHalo->Gather(o_sendBuf, o_v, k, op, Trans);

    if (exchange->gpu_aware) {
      //prepare MPI exchange
      exchange->DeviceStart(ogsType<T>::get(), k, op, Trans);
    } else {
      //get current stream
      device_t &device = platform.device;
      stream_t currentStream = device.getStream();

      //if not using gpu-aware mpi move the halo buffer to the host
      pinnedMemory<T> sendBuf = exchange->getHostSendBuffer();

      //wait for o_sendBuf to be ready
      device.finish();

      //queue copy to host
      device.setStream(dataStream);
      sendBuf.copyFrom(o_sendBuf, NhaloT*k,
                       0, properties_t("async", true));
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

  //queue local g operation
  gatherLocal->Gather(o_gv, o_v, k, op, trans);

  if (trans==Trans) { //if trans!=ogs::Trans theres no comms required
    if (exchange->gpu_aware) {
      //finish MPI exchange
      exchange->DeviceFinish(ogsType<T>::get(),k, op, Trans);

      deviceMemory<T> o_recvBuf = exchange->getDeviceRecvBuffer();

      //put the result at the end of o_gv
      o_recvBuf.copyTo(o_gv + k*NlocalT,
                       k*NhaloP, 0, properties_t("async", true));
    } else {
      //get current stream
      device_t &device = platform.device;
      stream_t currentStream = device.getStream();

      //synchronize data stream to ensure the buffer is on the host
      device.setStream(dataStream);
      device.finish();

      /*MPI exchange of host buffer*/
      exchange->HostStart (ogsType<T>::get(), k, op, trans);
      exchange->HostFinish(ogsType<T>::get(), k, op, trans);

      pinnedMemory<T> recvBuf = exchange->getHostRecvBuffer();

      // copy recv back to device
      //put the result at the end of o_gv
      recvBuf.copyTo(o_gv + k*NlocalT, k*NhaloP,
                     0, properties_t("async", true));
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
    pinnedMemory<T> sendBuf = exchange->getHostSendBuffer();

    //collect halo buffer
    gatherHalo->Gather(sendBuf, v, k, op, Trans);

    //prepare MPI exchange
    exchange->HostStart(ogsType<T>::get(), k, op, Trans);
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
    //finish MPI exchange
    exchange->HostFinish(ogsType<T>::get(), k, op, Trans);

    pinnedMemory<T> recvBuf = exchange->getHostRecvBuffer();

    //put the result at the end of o_gv
    recvBuf.copyTo(gv+k*NlocalT, k*NhaloP);
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

  if (trans==NoTrans) { //if trans!=ogs::NoTrans theres no comms required
    exchange->AllocBuffer(k*sizeof(T));

    deviceMemory<T> o_sendBuf = exchange->getDeviceSendBuffer();

    device_t &device = platform.device;

    if (exchange->gpu_aware) {
      //collect halo buffer
      o_sendBuf.copyFrom(o_gv + k*NlocalT,
                         k*NhaloP, 0, properties_t("async", true));

      //wait for o_sendBuf to be ready
      device.finish();

      //prepare MPI exchange
      exchange->DeviceStart(ogsType<T>::get(), k, Add, NoTrans);
    } else {
      //get current stream
      stream_t currentStream = device.getStream();

      //if not using gpu-aware mpi move the halo buffer to the host
      pinnedMemory<T> sendBuf = exchange->getHostSendBuffer();

      //wait for o_gv to be ready
      device.finish();

      //queue copy to host
      device.setStream(dataStream);
      sendBuf.copyFrom(o_gv + k*NlocalT, NhaloP*k,
                       0, properties_t("async", true));
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

  //queue local s operation
  gatherLocal->Scatter(o_v, o_gv, k, trans);

  if (trans==NoTrans) { //if trans!=ogs::NoTrans theres no comms required
    if (exchange->gpu_aware) {
      //finish MPI exchange
      exchange->DeviceFinish(ogsType<T>::get(),k, Add, NoTrans);
    } else {

      //get current stream
      device_t &device = platform.device;
      stream_t currentStream = device.getStream();

      //synchronize data stream to ensure the buffer is on the host
      device.setStream(dataStream);
      device.finish();

      /*MPI exchange of host buffer*/
      exchange->HostStart (ogsType<T>::get(), k, Add, NoTrans);
      exchange->HostFinish(ogsType<T>::get(), k, Add, NoTrans);

      pinnedMemory<T> recvBuf = exchange->getHostRecvBuffer();
      deviceMemory<T> o_recvBuf = exchange->getDeviceRecvBuffer();

      // copy recv back to device
      recvBuf.copyTo(o_recvBuf, NhaloT*k,
                     0, properties_t("async", true));
      device.finish(); //wait for transfer to finish
      device.setStream(currentStream);
    }

    deviceMemory<T> o_recvBuf = exchange->getDeviceRecvBuffer();

    //scatter halo buffer
    gatherHalo->Scatter(o_v, o_recvBuf, k, NoTrans);
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
    pinnedMemory<T> sendBuf = exchange->getHostSendBuffer();

    //collect halo buffer
    sendBuf.copyFrom(gv + k*NlocalT, k*NhaloP);

    //prepare MPI exchange
    exchange->HostStart(ogsType<T>::get(), k, Add, NoTrans);
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
    //finish MPI exchange (and put the result at the end of o_gv)
    exchange->HostFinish(ogsType<T>::get(), k, Add, NoTrans);

    pinnedMemory<T> recvBuf = exchange->getHostRecvBuffer();

    //scatter halo buffer
    gatherHalo->Scatter(v, recvBuf, k, NoTrans);
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
