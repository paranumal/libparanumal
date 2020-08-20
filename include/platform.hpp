/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "core.hpp"
#include "settings.hpp"
#include "linAlg.hpp"

class platformSettings_t: public settings_t {
public:
  platformSettings_t(MPI_Comm _comm);
  void report();
};

class platform_t {
public:
  const MPI_Comm& comm;
  platformSettings_t& settings;
  occa::properties props;

  occa::device device;
  linAlg_t linAlg;

  int rank, size;

  platform_t(platformSettings_t& _settings):
    comm(_settings.comm),
    settings(_settings) {

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank==0) {
      std::cout << "\n";
      std::cout << "\033[1m";
      std::cout << " _ _ _     ____                                             _ \n";
      std::cout << "| (_) |__ |  _ \\ __ _ _ __ __ _ _ __  _   _ _ __ ___   __ _| |\n";
      std::cout << "| | | '_ \\| |_) / _` | '__/ _` | '_ \\| | | | '_ ` _ \\ / _` | |\n";
      std::cout << "| | | |_) |  __/ (_| | | | (_| | | | | |_| | | | | | | (_| | |\n";
      std::cout << "|_|_|_.__/|_|   \\__,_|_|  \\__,_|_| |_|\\__,_|_| |_| |_|\\__,_|_|\n";
      std::cout << "\033[0m";
      std::cout << "\n";
      std::cout << "Version: 0.3.1\n";
      std::cout << "Contributing developers: Noel Chalmers, Ali Karakus, Kasia Swirydowicz,\n";
      std::cout << "                         Anthony Austin, & Tim Warburton\n";
      std::cout << "\n";
    }

    DeviceConfig();
    DeviceProperties();

    linAlg.Setup(this);
  }

  ~platform_t(){}

  occa::kernel buildKernel(std::string fileName, std::string kernelName,
                           occa::properties& kernelInfo);

  occa::memory malloc(const size_t bytes,
                      const void *src = NULL,
                      const occa::properties &prop = occa::properties()) {
    return device.malloc(bytes, src, prop);
  }

  occa::memory malloc(const size_t bytes,
                      const occa::memory &src,
                      const occa::properties &prop = occa::properties()) {
    return device.malloc(bytes, src, prop);
  }

  occa::memory malloc(const size_t bytes,
                      const occa::properties &prop) {
    return device.malloc(bytes, prop);
  }

  void *hostMalloc(const size_t bytes,
                   const void *src,
                   occa::memory &h_mem){
    occa::properties prop;
    prop["mapped"] = true;
    h_mem = device.malloc(bytes, prop);
    return h_mem.ptr(prop);
  }

private:
  void DeviceConfig();
  void DeviceProperties();

};

#endif