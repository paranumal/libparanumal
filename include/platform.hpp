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

#ifndef PLATFORM_HPP
#define PLATFORM_HPP

#define LIBP_MAJOR_VERSION 0
#define LIBP_MINOR_VERSION 5
#define LIBP_PATCH_VERSION 0
#define LIBP_VERSION       00500
#define LIBP_VERSION_STR   "0.5.0"

#include "core.hpp"
#include "memory.hpp"
#include "comm.hpp"
#include "settings.hpp"
#include "linAlg.hpp"

namespace libp {

class platformSettings_t: public settings_t {
public:
  platformSettings_t(comm_t _comm);
  void report();
};

namespace internal {

class iplatform_t {
public:
  platformSettings_t settings;
  properties_t props;

  iplatform_t(platformSettings_t& _settings):
    settings(_settings) {
  }
};

} //namespace internal

class platform_t {
 public:
  comm_t comm;
  device_t device;

 private:
  std::shared_ptr<internal::iplatform_t> iplatform;
  std::shared_ptr<linAlg_t<dfloat>> ilinAlg;

  memPool_t deviceMemPool;
  memPool_t pinnedMemPool;

public:
  platform_t()=default;

  platform_t(platformSettings_t& _settings) {

    iplatform = std::make_shared<internal::iplatform_t>(_settings);

    comm = settings().comm;

    if (comm.rank()==0) {
      std::cout << "\n";
      std::cout << "\033[1m";
      std::cout << " _ _ _     ____                                             _ \n";
      std::cout << "| (_) |__ |  _ \\ __ _ _ __ __ _ _ __  _   _ _ __ ___   __ _| |\n";
      std::cout << "| | | '_ \\| |_) / _` | '__/ _` | '_ \\| | | | '_ ` _ \\ / _` | |\n";
      std::cout << "| | | |_) |  __/ (_| | | | (_| | | | | |_| | | | | | | (_| | |\n";
      std::cout << "|_|_|_.__/|_|   \\__,_|_|  \\__,_|_| |_|\\__,_|_| |_| |_|\\__,_|_|\n";
      std::cout << "\033[0m";
      std::cout << "\n";
      std::cout << "Version: " LIBP_VERSION_STR " \n";
      std::cout << "Contributing developers: Noel Chalmers, Ali Karakus, Kasia Swirydowicz,\n";
      std::cout << "                         Anthony Austin, & Tim Warburton\n";
      std::cout << "\n";
    }

    DeviceConfig();
    DeviceProperties();

    properties_t props;
#if defined(LIBP_DEBUG)
    props["verbose"] = true;
#endif

    deviceMemPool = device.createMemoryPool(props);

    props["host"] = true;
    pinnedMemPool = device.createMemoryPool(props);

    ilinAlg = std::make_shared<linAlg_t<dfloat>>(this);
  }

  platform_t(const platform_t &other)=default;
  platform_t& operator = (const platform_t &other)=default;

  bool isInitialized() {
    return (iplatform!=nullptr);
  }

  void assertInitialized() {
    LIBP_ABORT("Platform not initialized.",
                  !isInitialized());
  }

  kernel_t buildKernel(std::string fileName, std::string kernelName,
                       properties_t& kernelInfo);

  template <typename T>
  deviceMemory<T> malloc(const size_t count,
                         const properties_t &prop = properties_t()) {
    assertInitialized();
    if (occa::dtype::get<T>() == occa::dtype::none) {
      return deviceMemory<T>(device.malloc(count*sizeof(T), occa::dtype::byte, prop));
    } else {
      return deviceMemory<T>(device.malloc<T>(count, prop));
    }
  }

  template <typename T>
  deviceMemory<T> malloc(const size_t count,
                         const memory<T> src,
                         const properties_t &prop = properties_t()) {
    assertInitialized();
    if (occa::dtype::get<T>() == occa::dtype::none) {
      return deviceMemory<T>(device.malloc(count*sizeof(T), occa::dtype::byte, src.ptr(), prop));
    } else {
      return deviceMemory<T>(device.malloc<T>(count, src.ptr(), prop));
    }
  }

  template <typename T>
  deviceMemory<T> malloc(const memory<T> src,
                         const properties_t &prop = properties_t()) {
    assertInitialized();
    if (occa::dtype::get<T>() == occa::dtype::none) {
      return deviceMemory<T>(device.malloc(src.size(), occa::dtype::byte, src.ptr(), prop));
    } else {
      return deviceMemory<T>(device.malloc<T>(src.length(), src.ptr(), prop));
    }
  }

  template <typename T>
  pinnedMemory<T> hostMalloc(const size_t count){
    assertInitialized();
    properties_t hostProp("host", true);
    if (occa::dtype::get<T>() == occa::dtype::none) {
      return pinnedMemory<T>(device.malloc(count*sizeof(T), occa::dtype::byte, nullptr, hostProp));
    } else {
      return pinnedMemory<T>(device.malloc<T>(count, nullptr, hostProp));
    }
  }

  template <typename T>
  pinnedMemory<T> hostMalloc(const size_t count,
                             const memory<T> src){
    assertInitialized();
    properties_t hostProp("host", true);
    if (occa::dtype::get<T>() == occa::dtype::none) {
      return pinnedMemory<T>(device.malloc(count*sizeof(T), occa::dtype::byte, src.ptr(), hostProp));
    } else {
      return pinnedMemory<T>(device.malloc<T>(count, src.ptr(), hostProp));
    }
  }

  template <typename T>
  pinnedMemory<T> hostMalloc(const memory<T> src){
    assertInitialized();
    properties_t hostProp("host", true);
    if (occa::dtype::get<T>() == occa::dtype::none) {
      return pinnedMemory<T>(device.malloc(src.size(), occa::dtype::byte, src.ptr(), hostProp));
    } else {
      return pinnedMemory<T>(device.malloc<T>(src.length(), src.ptr(), hostProp));
    }
  }

  template <typename T>
  deviceMemory<T> reserve(const size_t count) {
    assertInitialized();
    if (occa::dtype::get<T>() == occa::dtype::none) {
      return deviceMemory<T>(deviceMemPool.reserve(count*sizeof(T), occa::dtype::byte));
    } else {
      return deviceMemory<T>(deviceMemPool.reserve<T>(count));
    }
  }

  template <typename T>
  pinnedMemory<T> hostReserve(const size_t count) {
    assertInitialized();
    if (occa::dtype::get<T>() == occa::dtype::none) {
      return pinnedMemory<T>(pinnedMemPool.reserve(count*sizeof(T), occa::dtype::byte));
    } else {
      return pinnedMemory<T>(pinnedMemPool.reserve<T>(count));
    }
  }

  /*Return the alignment of the memory pool as count of type Ts*/
  template<typename T = std::byte>
  size_t memPoolAlignment() {
    return (deviceMemPool.alignment() + sizeof(T)-1)/sizeof(T);
  }

  linAlg_t<dfloat>& linAlg() {
    assertInitialized();
    return *ilinAlg;
  }

  settings_t& settings() {
    assertInitialized();
    return iplatform->settings;
  }

  properties_t& props() {
    assertInitialized();
    return iplatform->props;
  }

  void finish() {
    device.finish();
  }

  stream_t getStream() {
    return device.getStream();
  }

  void setStream(stream_t stream) {
    device.setStream(stream);
  }

  const int rank() const {
    return comm.rank();
  }

  const int size() const {
    return comm.size();
  }

  int getDeviceCount(const std::string mode) {
    return occa::getDeviceCount(mode);
  }

  void setCacheDir(const std::string cacheDir) {
    occa::env::setOccaCacheDir(cacheDir);
  }

 private:
  void DeviceConfig();
  void DeviceProperties();
};

} //namespace libp

#endif
