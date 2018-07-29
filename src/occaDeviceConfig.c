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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include <unistd.h>
#include  "mpi.h"
#include "mesh.h"

void occaDeviceConfig(mesh_t *mesh, setupAide &options){

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  rank = mesh->rank;
  size = mesh->size;

  long int hostId = gethostid();

  long int* hostIds = (long int*) calloc(size,sizeof(long int));
  MPI_Allgather(&hostId,1,MPI_LONG,hostIds,1,MPI_LONG,mesh->comm);

  int device_id = 0;
  int totalDevices = 0;
  for (int r=0;r<rank;r++) {
    if (hostIds[r]==hostId) device_id++;
  }
  for (int r=0;r<size;r++) {
    if (hostIds[r]==hostId) totalDevices++;
  }

  if (size==1) options.getArgs("DEVICE NUMBER" ,device_id);

#ifdef OCCA_VERSION_1_0
  // read thread model/device/platform from options
  if(options.compareArgs("THREAD MODEL", "CUDA")){
    sprintf(deviceConfig, "mode: 'CUDA', device_id: %d",device_id);
  }
  else if(options.compareArgs("THREAD MODEL", "HIP")){
    sprintf(deviceConfig, "mode: 'HIP', device_id: %d",device_id);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenCL")){
    int plat;
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "mode: 'OpenCL', device_id: %d, platform_id: %d", device_id, plat);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenMP")){
    sprintf(deviceConfig, "mode: 'OpenMP' ");
  }
  else{
    sprintf(deviceConfig, "mode: 'Serial' ");
  }
#else
  // read thread model/device/platform from options
  if(options.compareArgs("THREAD MODEL", "CUDA")){
    sprintf(deviceConfig, "mode=CUDA, deviceID=%d",device_id);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenCL")){
    int plat;
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "mode=OpenCL, deviceID=%d, platformID=%d", device_id, plat);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenMP")){
    sprintf(deviceConfig, "mode=OpenMP ");
  }
  else{
    sprintf(deviceConfig, "mode=Serial");
  }
#endif

  //set number of omp threads to use
  int Ncores = sysconf(_SC_NPROCESSORS_ONLN);
  int Nthreads = Ncores/totalDevices;

  //  Nthreads = mymax(1,Nthreads/2);
  Nthreads = mymax(1,Nthreads/2);
  omp_set_num_threads(Nthreads);

  if (rank==0 && options.compareArgs("VERBOSE","TRUE"))
    printf("Rank %d: Ncores = %d, Nthreads = %d\n", rank, Ncores, Nthreads);
  
  mesh->device.setup(deviceConfig);

  occa::initTimer(mesh->device);
}
