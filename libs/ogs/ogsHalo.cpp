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
 * Device Exchange
 ********************************/
template<typename T>
void halo_t::Exchange(deviceMemory<T> o_v, const int k) {
  ExchangeStart (o_v, k);
  ExchangeFinish(o_v, k);
}

template<typename T>
void halo_t::ExchangeStart(deviceMemory<T> o_v, const int k){
  exchange->AllocBuffer(k*sizeof(T));

  deviceMemory<T> o_sendBuf = exchange->getDeviceSendBuffer();

  if (exchange->gpu_aware) {
    if (gathered_halo) {
      //if this halo was build from a gathered ogs the halo nodes are at the end
      o_sendBuf.copyFrom(o_v + k*NlocalT, k*NhaloP,
                         0, properties_t("async", true));
    } else {
      //collect halo buffer
      gatherHalo->Gather(o_sendBuf, o_v, k, Add, NoTrans);
    }

    //prepare MPI exchange
    exchange->DeviceStart(ogsType<T>::get(), k, Add, NoTrans);

  } else {
    //get current stream
    device_t &device = platform.device;
    stream_t currentStream = device.getStream();

    //if not using gpu-aware mpi move the halo buffer to the host
    pinnedMemory<T> sendBuf = exchange->getHostSendBuffer();

    if (gathered_halo) {
      //wait for o_v to be ready
      device.finish();

      //queue copy to host
      device.setStream(dataStream);
      sendBuf.copyFrom(o_v + k*NlocalT, NhaloP*k,
                       0, properties_t("async", true));
      device.setStream(currentStream);
    } else {
      //collect halo buffer
      gatherHalo->Gather(o_sendBuf, o_v, k, Add, NoTrans);

      //wait for o_haloBuf to be ready
      device.finish();

      //queue copy to host
      device.setStream(dataStream);
      sendBuf.copyFrom(o_sendBuf, NhaloP*k,
                       0, properties_t("async", true));
      device.setStream(currentStream);
    }
  }
}

template<typename T>
void halo_t::ExchangeFinish(deviceMemory<T> o_v, const int k){


  //write exchanged halo buffer back to vector
  if (exchange->gpu_aware) {
    //finish MPI exchange
    exchange->DeviceFinish(ogsType<T>::get(), k, Add, NoTrans);

    deviceMemory<T> o_recvBuf = exchange->getDeviceRecvBuffer();

    if (gathered_halo) {
      o_recvBuf.copyTo(o_v + k*(NlocalT+NhaloP), k*Nhalo,
                       k*NhaloP, properties_t("async", true));
    } else {
      gatherHalo->Scatter(o_v, o_recvBuf, k, NoTrans);
    }
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
    if (gathered_halo) {
      recvBuf.copyTo(o_v + k*(NlocalT+NhaloP), k*Nhalo,
                     k*NhaloP, properties_t("async", true));
      device.finish(); //wait for transfer to finish
      device.setStream(currentStream);
    } else {
      recvBuf.copyTo(o_recvBuf, k*NhaloT,
                     0, properties_t("async", true));
      device.finish(); //wait for transfer to finish
      device.setStream(currentStream);

      gatherHalo->Scatter(o_v, o_recvBuf, k, NoTrans);
    }
  }
}

template void halo_t::Exchange(deviceMemory<float> o_v, const int k);
template void halo_t::Exchange(deviceMemory<double> o_v, const int k);
template void halo_t::Exchange(deviceMemory<int> o_v, const int k);
template void halo_t::Exchange(deviceMemory<long long int> o_v, const int k);

//host version
template<typename T>
void halo_t::Exchange(memory<T> v, const int k) {
  ExchangeStart (v, k);
  ExchangeFinish(v, k);
}

template<typename T>
void halo_t::ExchangeStart(memory<T> v, const int k) {
  exchange->AllocBuffer(k*sizeof(T));

  pinnedMemory<T> sendBuf = exchange->getHostSendBuffer();

  //collect halo buffer
  if (gathered_halo) {
    //if this halo was build from a gathered ogs the halo nodes are at the end
    sendBuf.copyFrom(v + k*NlocalT, k*NhaloP);
  } else {
    gatherHalo->Gather(sendBuf, v, k, Add, NoTrans);
  }

  //Prepare MPI exchange
  exchange->HostStart(ogsType<T>::get(), k, Add, NoTrans);
}

template<typename T>
void halo_t::ExchangeFinish(memory<T> v, const int k) {

  //finish MPI exchange
  exchange->HostFinish(ogsType<T>::get(), k, Add, NoTrans);

  pinnedMemory<T> recvBuf = exchange->getHostRecvBuffer();

  //write exchanged halo buffer back to vector
  if (gathered_halo) {
    //if this halo was build from a gathered ogs the halo nodes are at the end
    recvBuf.copyTo(v + k*(NlocalT+NhaloP),
                   k*Nhalo,
                   k*NhaloP);
  } else {
    gatherHalo->Scatter(v, recvBuf, k, NoTrans);
  }
}

template void halo_t::Exchange(memory<float> v, const int k);
template void halo_t::Exchange(memory<double> v, const int k);
template void halo_t::Exchange(memory<int> v, const int k);
template void halo_t::Exchange(memory<long long int> v, const int k);

/********************************
 * Combine
 ********************************/
template<typename T>
void halo_t::Combine(deviceMemory<T> o_v, const int k) {
  CombineStart (o_v, k);
  CombineFinish(o_v, k);
}

template<typename T>
void halo_t::CombineStart(deviceMemory<T> o_v, const int k){
  exchange->AllocBuffer(k*sizeof(T));

  deviceMemory<T> o_sendBuf = exchange->getDeviceSendBuffer();

  if (exchange->gpu_aware) {
    if (gathered_halo) {
      //if this halo was build from a gathered ogs the halo nodes are at the end
      o_sendBuf.copyFrom(o_v + k*NlocalT, k*NhaloT,
                         0, properties_t("async", true));
    } else {
      //collect halo buffer
      gatherHalo->Gather(o_sendBuf, o_v, k, Add, Trans);
    }

    //prepare MPI exchange
    exchange->DeviceStart(ogsType<T>::get(), k, Add, Trans);
  } else {
    //get current stream
    device_t &device = platform.device;
    stream_t currentStream = device.getStream();

    //if not using gpu-aware mpi move the halo buffer to the host
    pinnedMemory<T> sendBuf = exchange->getHostSendBuffer();

    if (gathered_halo) {
      //wait for o_v to be ready
      device.finish();

      //queue copy to host
      device.setStream(dataStream);
      sendBuf.copyFrom(o_v + k*NlocalT, NhaloT*k,
                       0, properties_t("async", true));
      device.setStream(currentStream);
    } else {
      //collect halo buffer
      gatherHalo->Gather(o_sendBuf, o_v, k, Add, Trans);

      //wait for o_sendBuf to be ready
      device.finish();

      //queue copy to host
      device.setStream(dataStream);
      sendBuf.copyFrom(o_sendBuf, NhaloT*k,
                       0, properties_t("async", true));
      device.setStream(currentStream);
    }
  }
}

template<typename T>
void halo_t::CombineFinish(deviceMemory<T> o_v, const int k){


  //write exchanged halo buffer back to vector
  if (exchange->gpu_aware) {
    //finish MPI exchange
    exchange->DeviceFinish(ogsType<T>::get(), k, Add, Trans);

    deviceMemory<T> o_recvBuf = exchange->getDeviceRecvBuffer();

    if (gathered_halo) {
      //if this halo was build from a gathered ogs the halo nodes are at the end
      o_recvBuf.copyTo(o_v + k*NlocalT, k*NhaloP,
                       0, properties_t("async", true));
    } else {
      gatherHalo->Scatter(o_v, o_recvBuf, k, Trans);
    }
  } else {

    //get current stream
    device_t &device = platform.device;
    stream_t currentStream = device.getStream();

    //synchronize data stream to ensure the buffer is on the host
    device.setStream(dataStream);
    device.finish();

    /*MPI exchange of host buffer*/
    exchange->HostStart (ogsType<T>::get(), k, Add, Trans);
    exchange->HostFinish(ogsType<T>::get(), k, Add, Trans);

    pinnedMemory<T> recvBuf = exchange->getHostRecvBuffer();
    deviceMemory<T> o_recvBuf = exchange->getDeviceRecvBuffer();

    if (gathered_halo) {
      // copy recv back to device
      recvBuf.copyTo(o_v + k*NlocalT, NhaloP*k,
                     0, properties_t("async", true));
      device.finish(); //wait for transfer to finish
      device.setStream(currentStream);
    } else {
      recvBuf.copyTo(o_recvBuf, NhaloP*k,
                     0, properties_t("async", true));
      device.finish(); //wait for transfer to finish
      device.setStream(currentStream);

      gatherHalo->Scatter(o_v, o_recvBuf, k, Trans);
    }
  }
}

template void halo_t::Combine(deviceMemory<float> o_v, const int k);
template void halo_t::Combine(deviceMemory<double> o_v, const int k);
template void halo_t::Combine(deviceMemory<int> o_v, const int k);
template void halo_t::Combine(deviceMemory<long long int> o_v, const int k);

//host version
template<typename T>
void halo_t::Combine(memory<T> v, const int k) {
  CombineStart (v, k);
  CombineFinish(v, k);
}

template<typename T>
void halo_t::CombineStart(memory<T> v, const int k) {
  exchange->AllocBuffer(k*sizeof(T));

  pinnedMemory<T> sendBuf = exchange->getHostSendBuffer();

  //collect halo buffer
  if (gathered_halo) {
    //if this halo was build from a gathered ogs the halo nodes are at the end
    sendBuf.copyFrom(v + k*NlocalT, k*NhaloT);
  } else {
    gatherHalo->Gather(sendBuf, v, k, Add, Trans);
  }

  //Prepare MPI exchange
  exchange->HostStart(ogsType<T>::get(), k, Add, Trans);
}


template<typename T>
void halo_t::CombineFinish(memory<T> v, const int k) {

  //finish MPI exchange
  exchange->HostFinish(ogsType<T>::get(), k, Add, Trans);

  pinnedMemory<T> recvBuf = exchange->getHostRecvBuffer();

  //write exchanged halo buffer back to vector
  if (gathered_halo) {
    //if this halo was build from a gathered ogs the halo nodes are at the end
    recvBuf.copyTo(v + k*NlocalT, k*NhaloP);
  } else {
    gatherHalo->Scatter(v, recvBuf, k, Trans);
  }
}

template void halo_t::Combine(memory<float> v, const int k);
template void halo_t::Combine(memory<double> v, const int k);
template void halo_t::Combine(memory<int> v, const int k);
template void halo_t::Combine(memory<long long int> v, const int k);

} //namespace ogs

} //namespace libp
