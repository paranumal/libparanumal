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

#include <unistd.h>
#include <mpi.h>
#include "occa.hpp"
#include "types.h"
#include "utils.hpp"
#include "core.hpp"
#include "settings.hpp"
#include "omp.h"

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

void occaDeviceConfig(occa::device &device, MPI_Comm comm,
                      settings_t& settings, occa::properties& props){

  // OCCA build stuff
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  char deviceConfig[BUFSIZ];

  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);
  long int hostId = gethostid();

  long int* hostIds = (long int*) calloc(size,sizeof(long int));
  MPI_Allgather(&hostId,1,MPI_LONG,hostIds,1,MPI_LONG,comm);

  int device_id = 0;
  int totalDevices = 0;
  for (int r=0;r<rank;r++) {
    if (hostIds[r]==hostId) device_id++;
  }
  for (int r=0;r<size;r++) {
    if (hostIds[r]==hostId) totalDevices++;
  }

  //for testing a single device, run with 1 rank and specify DEVICE NUMBER
  if (size==1) settings.getSetting("DEVICE NUMBER",device_id);

  //check for over-subscribing devices
  if(settings.compareSetting("THREAD MODEL", "CUDA")){
#if OCCA_CUDA_ENABLED
    int deviceCount = occa::cuda::getDeviceCount();
    if (device_id>=deviceCount) {
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
    if (device_id>=deviceCount) {
      stringstream ss;
      ss << "Rank " << rank << " oversubscribing HIP device " << device_id%deviceCount << " on node \"" << hostname<< "\"";
      LIBP_WARNING(ss.str());
      device_id = device_id%deviceCount;
    }
#endif

  }
  else if(settings.compareSetting("THREAD MODEL", "OpenCL")){
#if OCCA_OPENCL_ENABLED
    int plat;
    settings.getSetting("PLATFORM NUMBER", plat);
    int deviceCount = occa::opencl::getDeviceCountInPlatform(plat);
    if (device_id>=deviceCount) {
      stringstream ss;
      ss << "Rank " << rank << " oversubscribing OpenCL device " << device_id%deviceCount << " on node \"" << hostname<< "\"";
      LIBP_WARNING(ss.str());
      device_id = device_id%deviceCount;
    }
#endif
  }

  // read thread model/device/platform from settings
  if(settings.compareSetting("THREAD MODEL", "CUDA")){
    sprintf(deviceConfig, "mode: 'CUDA', device_id: %d", device_id);
  }
  else if(settings.compareSetting("THREAD MODEL", "HIP")){
    sprintf(deviceConfig, "mode: 'HIP', device_id: %d", device_id);
  }
  else if(settings.compareSetting("THREAD MODEL", "OpenCL")){
    int plat;
    settings.getSetting("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "mode: 'OpenCL', device_id: %d, platform_id: %d", device_id, plat);
  }
  else if(settings.compareSetting("THREAD MODEL", "OpenMP")){
    sprintf(deviceConfig, "mode: 'OpenMP' ");
  }
  else{
    sprintf(deviceConfig, "mode: 'Serial' ");
  }

  //set number of omp threads to use
  int Ncores = sysconf(_SC_NPROCESSORS_ONLN);
  int Nthreads = Ncores/totalDevices;

  //  Nthreads = mymax(1,Nthreads/2);
  Nthreads = mymax(1,Nthreads/2);
  omp_set_num_threads(Nthreads);

  // if (settings.compareSetting("VERBOSE","TRUE"))
  //   printf("Rank %d: Ncores = %d, Nthreads = %d, device_id = %d \n", rank, Ncores, Nthreads, device_id);

  //  mesh->device.setup( (std::string) deviceConfig); // deviceProps);
  device.setup( (std::string)deviceConfig);

#ifdef USE_OCCA_MEM_BYTE_ALIGN
  // change OCCA MEM BYTE ALIGNMENT
  occa::env::OCCA_MEM_BYTE_ALIGN = USE_OCCA_MEM_BYTE_ALIGN;
#endif

  occaDeviceProperties(device, props);
}
