
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
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  long int hostId = gethostid();

  long int* hostIds = (long int*) calloc(size,sizeof(long int));
  MPI_Allgather(&hostId,1,MPI_LONG,hostIds,1,MPI_LONG,MPI_COMM_WORLD);

  int deviceID = 0;
  int totalDevices = 0;
  for (int r=0;r<rank;r++) {
    if (hostIds[r]==hostId) deviceID++;
  }
  for (int r=0;r<size;r++) {
    if (hostIds[r]==hostId) totalDevices++;
  }

  if (size==1) options.getArgs("DEVICE NUMBER" ,deviceID);

  // read thread model/device/platform from options
  if(options.compareArgs("THREAD MODEL", "CUDA")){
    sprintf(deviceConfig, "mode = CUDA, deviceID = %d",deviceID);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenCL")){
    int plat;
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "mode = OpenCL, deviceID = %d, platformID = %d", deviceID, plat);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenMP")){
    sprintf(deviceConfig, "mode = OpenMP");
  }
  else{
    sprintf(deviceConfig, "mode = Serial");
  }

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
