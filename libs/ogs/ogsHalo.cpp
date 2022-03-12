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

  deviceMemory<T> o_haloBuf = exchange->o_workspace;

  if (exchange->gpu_aware) {
    if (gathered_halo) {
      //if this halo was build from a gathered ogs the halo nodes are at the end
      o_haloBuf.copyFrom(o_v + k*NlocalT, k*NhaloP,
                         0, "async: true");
    } else {
      //collect halo buffer
      gatherHalo->Gather(o_haloBuf, o_v, k, Add, NoTrans);
    }

    //prepare MPI exchange
    exchange->Start(o_haloBuf, k, Add, NoTrans);

  } else {
    //get current stream
    device_t &device = platform.device;
    stream_t currentStream = device.getStream();

    //if not using gpu-aware mpi move the halo buffer to the host
    pinnedMemory<T> haloBuf = exchange->h_workspace;

    if (gathered_halo) {
      //wait for o_v to be ready
      device.finish();

      //queue copy to host
      device.setStream(dataStream);
      haloBuf.copyFrom(o_v + k*NlocalT, NhaloP*k,
                       0, "async: true");
      device.setStream(currentStream);
    } else {
      //collect halo buffer
      gatherHalo->Gather(o_haloBuf, o_v, k, Add, NoTrans);

      //wait for o_haloBuf to be ready
      device.finish();

      //queue copy to host
      device.setStream(dataStream);
      haloBuf.copyFrom(o_haloBuf, NhaloP*k,
                       0, "async: true");
      device.setStream(currentStream);
    }
  }
}

template<typename T>
void halo_t::ExchangeFinish(deviceMemory<T> o_v, const int k){

  deviceMemory<T> o_haloBuf = exchange->o_workspace;

  //write exchanged halo buffer back to vector
  if (exchange->gpu_aware) {
    //finish MPI exchange
    exchange->Finish(o_haloBuf, k, Add, NoTrans);

    if (gathered_halo) {
      o_haloBuf.copyTo(o_v + k*(NlocalT+NhaloP), k*Nhalo,
                       k*NhaloP, "async: true");
    } else {
      gatherHalo->Scatter(o_v, o_haloBuf, k, NoTrans);
    }
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
    if (gathered_halo) {
      haloBuf.copyTo(o_v + k*(NlocalT+NhaloP), k*Nhalo,
                     k*NhaloP, "async: true");
      device.finish(); //wait for transfer to finish
      device.setStream(currentStream);
    } else {
      haloBuf.copyTo(o_haloBuf+k*NhaloP, k*Nhalo,
                     k*NhaloP, "async: true");
      device.finish(); //wait for transfer to finish
      device.setStream(currentStream);

      gatherHalo->Scatter(o_v, o_haloBuf, k, NoTrans);
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

  pinnedMemory<T> haloBuf = exchange->h_workspace;

  //collect halo buffer
  if (gathered_halo) {
    //if this halo was build from a gathered ogs the halo nodes are at the end
    haloBuf.copyFrom(v + k*NlocalT, k*NhaloP);
  } else {
    gatherHalo->Gather(haloBuf, v, k, Add, NoTrans);
  }

  //Prepare MPI exchange
  exchange->Start(haloBuf, k, Add, NoTrans);
}

template<typename T>
void halo_t::ExchangeFinish(memory<T> v, const int k) {

  pinnedMemory<T> haloBuf = exchange->h_workspace;

  //finish MPI exchange
  exchange->Finish(haloBuf, k, Add, NoTrans);

  //write exchanged halo buffer back to vector
  if (gathered_halo) {
    //if this halo was build from a gathered ogs the halo nodes are at the end
    haloBuf.copyTo(v + k*(NlocalT+NhaloP),
                   k*Nhalo,
                   k*NhaloP);
  } else {
    gatherHalo->Scatter(v, haloBuf, k, NoTrans);
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

  deviceMemory<T> o_haloBuf = exchange->o_workspace;

  if (exchange->gpu_aware) {
    if (gathered_halo) {
      //if this halo was build from a gathered ogs the halo nodes are at the end
      o_haloBuf.copyFrom(o_v + k*NlocalT, k*NhaloT,
                         0, "async: true");
    } else {
      //collect halo buffer
      gatherHalo->Gather(o_haloBuf, o_v, k, Add, Trans);
    }

    //prepare MPI exchange
    exchange->Start(o_haloBuf, k, Add, Trans);
  } else {
    //get current stream
    device_t &device = platform.device;
    stream_t currentStream = device.getStream();

    //if not using gpu-aware mpi move the halo buffer to the host
    pinnedMemory<T> haloBuf = exchange->h_workspace;

    if (gathered_halo) {
      //wait for o_v to be ready
      device.finish();

      //queue copy to host
      device.setStream(dataStream);
      haloBuf.copyFrom(o_v + k*NlocalT, NhaloT*k,
                       0, "async: true");
      device.setStream(currentStream);
    } else {
      //collect halo buffer
      gatherHalo->Gather(o_haloBuf, o_v, k, Add, Trans);

      //wait for o_haloBuf to be ready
      device.finish();

      //queue copy to host
      device.setStream(dataStream);
      haloBuf.copyFrom(o_haloBuf, NhaloT*k,
                       0, "async: true");
      device.setStream(currentStream);
    }
  }
}

template<typename T>
void halo_t::CombineFinish(deviceMemory<T> o_v, const int k){

  deviceMemory<T> o_haloBuf = exchange->o_workspace;

  //write exchanged halo buffer back to vector
  if (exchange->gpu_aware) {
    //finish MPI exchange
    exchange->Finish(o_haloBuf, k, Add, Trans);

    if (gathered_halo) {
      //if this halo was build from a gathered ogs the halo nodes are at the end
      o_haloBuf.copyTo(o_v + k*NlocalT, k*NhaloP,
                       0, "async: true");
    } else {
      gatherHalo->Scatter(o_v, o_haloBuf, k, Trans);
    }
  } else {
    pinnedMemory<T> haloBuf = exchange->h_workspace;

    //get current stream
    device_t &device = platform.device;
    stream_t currentStream = device.getStream();

    //synchronize data stream to ensure the buffer is on the host
    device.setStream(dataStream);
    device.finish();

    /*MPI exchange of host buffer*/
    exchange->Start (haloBuf, k, Add, Trans);
    exchange->Finish(haloBuf, k, Add, Trans);

    if (gathered_halo) {
      // copy recv back to device
      haloBuf.copyTo(o_v + k*NlocalT, NhaloP*k,
                     0, "async: true");
      device.finish(); //wait for transfer to finish
      device.setStream(currentStream);
    } else {
      haloBuf.copyTo(o_haloBuf, NhaloP*k,
                     0, "async: true");
      device.finish(); //wait for transfer to finish
      device.setStream(currentStream);

      gatherHalo->Scatter(o_v, o_haloBuf, k, Trans);
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

  pinnedMemory<T> haloBuf = exchange->h_workspace;

  //collect halo buffer
  if (gathered_halo) {
    //if this halo was build from a gathered ogs the halo nodes are at the end
    haloBuf.copyFrom(v + k*NlocalT, k*NhaloT);
  } else {
    gatherHalo->Gather(haloBuf, v, k, Add, Trans);
  }

  //Prepare MPI exchange
  exchange->Start(haloBuf, k, Add, Trans);
}


template<typename T>
void halo_t::CombineFinish(memory<T> v, const int k) {

  pinnedMemory<T> haloBuf = exchange->h_workspace;

  //finish MPI exchange
  exchange->Finish(haloBuf, k, Add, Trans);

  //write exchanged halo buffer back to vector
  if (gathered_halo) {
    //if this halo was build from a gathered ogs the halo nodes are at the end
    haloBuf.copyTo(v + k*NlocalT, k*NhaloP);
  } else {
    gatherHalo->Scatter(v, haloBuf, k, Trans);
  }
}

template void halo_t::Combine(memory<float> v, const int k);
template void halo_t::Combine(memory<double> v, const int k);
template void halo_t::Combine(memory<int> v, const int k);
template void halo_t::Combine(memory<long long int> v, const int k);

} //namespace ogs

} //namespace libp
