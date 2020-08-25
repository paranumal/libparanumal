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

#include "platform.hpp"
// #include "omp.h"

//hack the hook to ask OCCA to return a device count
namespace occa {
#if OCCA_CUDA_ENABLED
  namespace cuda {
    int getDeviceCount();
  }
#endif
#if OCCA_HIP_ENABLED
  namespace hip {
    int getDeviceCount();
  }
#endif
#if OCCA_OPENCL_ENABLED
  namespace opencl {
    namespace info {
      static const int CPU     = (1 << 0);
      static const int GPU     = (1 << 1);
      static const int FPGA    = (1 << 3);
      static const int XeonPhi = (1 << 2);
      static const int anyType = (CPU | GPU | FPGA | XeonPhi);

      static const int Intel     = (1 << 4);
      static const int AMD       = (1 << 5);
      static const int Altera    = (1 << 6);
      static const int NVIDIA    = (1 << 7);
      static const int anyVendor = (Intel | AMD | Altera | NVIDIA);

      static const int any = (anyType | anyVendor);

      std::string deviceType(int type);
      std::string vendor(int type);
    }

    int getDeviceCountInPlatform(int pID, int type = info::any);
  }
#endif
}

// OCCA build stuff
void platform_t::DeviceConfig(){

  int plat=0;
  int device_id=0;

  //for testing a single device, run with 1 rank and specify DEVICE NUMBER
  if (size==1) {
    settings.getSetting("DEVICE NUMBER",device_id);
  } else {
    //find out how many ranks and devices are on this system
    char* hostnames = (char *) ::malloc(size*sizeof(char)*MPI_MAX_PROCESSOR_NAME);
    char* hostname = hostnames+rank*MPI_MAX_PROCESSOR_NAME;

    int namelen;
    MPI_Get_processor_name(hostname,&namelen);

    MPI_Allgather(hostname , MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
                  hostnames, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MPI_COMM_WORLD);

    int localRank = 0;
    int localSize = 0;
    for (int n=0; n<rank; n++){
      if (!strcmp(hostname, hostnames+n*MPI_MAX_PROCESSOR_NAME)) localRank++;
    }
    for (int n=0; n<size; n++){
      if (!strcmp(hostname, hostnames+n*MPI_MAX_PROCESSOR_NAME)) localSize++;
    }

    device_id = localRank;

    //check for over-subscribing devices
    if(settings.compareSetting("THREAD MODEL", "CUDA")){
#if OCCA_CUDA_ENABLED
      int deviceCount = occa::cuda::getDeviceCount();
      if (deviceCount>0 && localRank>=deviceCount) {
        stringstream ss;
        ss << "Rank " << rank << " oversubscribing CUDA device " << device_id%deviceCount << " on node \"" << hostname<< "\"";
        LIBP_WARNING(ss.str());
        device_id = device_id%deviceCount;
      }
#endif
    }
    else if(settings.compareSetting("THREAD MODEL", "HIP")){
#if OCCA_HIP_ENABLED
      int deviceCount = occa::hip::getDeviceCount();
      if (deviceCount>0 && localRank>=deviceCount) {
        stringstream ss;
        ss << "Rank " << rank << " oversubscribing HIP device " << device_id%deviceCount << " on node \"" << hostname<< "\"";
        LIBP_WARNING(ss.str());
        device_id = device_id%deviceCount;
      }
#endif
    }
    else if(settings.compareSetting("THREAD MODEL", "OpenCL")){
#if OCCA_OPENCL_ENABLED
      settings.getSetting("PLATFORM NUMBER", plat);
      int deviceCount = occa::opencl::getDeviceCountInPlatform(plat);
      if (deviceCount>0 && localRank>=deviceCount) {
        stringstream ss;
        ss << "Rank " << rank << " oversubscribing OpenCL device " << device_id%deviceCount << " on node \"" << hostname<< "\"";
        LIBP_WARNING(ss.str());
        device_id = device_id%deviceCount;
      }
#endif
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(hostnames);
  }

  // read thread model/device/platform from settings
  std::string mode;

  if(settings.compareSetting("THREAD MODEL", "CUDA")){
    mode = "mode: 'CUDA', device_id: " + std::to_string(device_id);
  }
  else if(settings.compareSetting("THREAD MODEL", "HIP")){
    mode = "mode: 'HIP', device_id: " + std::to_string(device_id);
  }
  else if(settings.compareSetting("THREAD MODEL", "OpenCL")){
    mode = "mode: 'OpenCL', platform_id : " + std::to_string(plat)
                          + ", device_id: " + std::to_string(device_id);
  }
  else if(settings.compareSetting("THREAD MODEL", "OpenMP")){
    mode = "mode: 'OpenMP'";
  }
  else{
    mode = "mode: 'Serial'";
  }

  //set number of omp threads to use
  //int Ncores = sysconf(_SC_NPROCESSORS_ONLN);
  //int Nthreads = Ncores/localSize;
  // Nthreads = mymax(1,Nthreads/2);
  // omp_set_num_threads(Nthreads);

  // if (settings.compareSetting("VERBOSE","TRUE"))
  //   printf("Rank %d: Ncores = %d, Nthreads = %d, device_id = %d \n", rank, Ncores, Nthreads, device_id);

  device.setup(mode);

#ifdef USE_OCCA_MEM_BYTE_ALIGN
  // change OCCA MEM BYTE ALIGNMENT
  occa::env::OCCA_MEM_BYTE_ALIGN = USE_OCCA_MEM_BYTE_ALIGN;
#endif
}
