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

#include "acoustics.hpp"

acoustics_t& acoustics_t::Setup(mesh_t& mesh, linAlg_t& linAlg){

  acoustics_t* acoustics = new acoustics_t(mesh, linAlg);

  settings_t& settings = acoustics->settings;

  acoustics->Nfields = (mesh.dim==3) ? 4:3;

  //setup timeStepper
  dlong Nlocal = mesh.Nelements*mesh.Np*acoustics->Nfields;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np*acoustics->Nfields;
  acoustics->timeStepper = timeStepper_t::Setup(Nlocal, Nhalo, *acoustics);
  acoustics->timeStepper->Init();

  // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat cfl = 0.5; // depends on the stability region size

  dfloat dt = cfl*hmin/((mesh.N+1.)*(mesh.N+1.));
  acoustics->timeStepper->SetTimeStep(dt);

  //setup linear algebra module
  acoustics->linAlg.InitKernels({"innerProd"}, mesh.comm);

  // set penalty parameter
  dfloat Lambda2 = 0.5;

  /*setup trace halo exchange */

  //get number of elements in each rank for global offset
  hlong *globalNelements = (hlong*) calloc(mesh.size+1,sizeof(hlong));
  MPI_Allgather(&(mesh.Nelements), 1, MPI_HLONG, globalNelements+1, 1, MPI_HLONG, mesh.comm);
  for(int r=0;r<mesh.size;++r) globalNelements[r+1] = globalNelements[r]+globalNelements[r+1];

  hlong globalOffset = globalNelements[mesh.rank];
  free(globalNelements);

  // populate a global numbering system
  hlong *globalIds = (hlong *) calloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np,sizeof(hlong));

  for (dlong e=0;e<mesh.Nelements;e++) {
    for (int n=0;n<mesh.Np;n++) {
      dlong id = e*mesh.Np + n;
      globalIds[id] = (e+globalOffset)*mesh.Np + n + 1;
    }
  }
  //populate the trace ids using the trace halo from mesh
  mesh.traceHalo->Exchange(globalIds, 1, ogs_hlong);

  //populate another global numbering system which has the Nfields stride
  hlong *traceGlobalIds = (hlong *) calloc((mesh.Nelements+mesh.totalHaloPairs)
                                           *mesh.Np*acoustics->Nfields,sizeof(hlong));
  for (dlong e=0;e<mesh.Nelements;e++) {
    for (int k=0;k<acoustics->Nfields;k++) {
      for (int n=0;n<mesh.Np;n++) {
        dlong id = e*mesh.Np*acoustics->Nfields + k*mesh.Np + n;
        traceGlobalIds[id] = (e+globalOffset)*mesh.Np*acoustics->Nfields + k*mesh.Np + n + 1;
      }
    }
  }

  //translate the halo ids from globalIds to traceGlobalIds
  for (dlong e=mesh.Nelements;e<mesh.Nelements+mesh.totalHaloPairs;e++) {
    for (int n=0;n<mesh.Np;n++) {
      if (globalIds[e*mesh.Np + n]!=0) {
        hlong gid = globalIds[e*mesh.Np + n] - 1;
        hlong eP = gid/mesh.Np;
        int vidP = gid%mesh.Np;

        dlong id = e*mesh.Np*acoustics->Nfields + n;
        for (int k=0;k<acoustics->Nfields;k++) {
          traceGlobalIds[id + k*mesh.Np] = -(eP*mesh.Np*acoustics->Nfields + k*mesh.Np + vidP + 1);
        }
      }
    }
  }
  free(globalIds);

  int verbose = 0;
  acoustics->traceHalo = halo_t::Setup((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*acoustics->Nfields,
                                        traceGlobalIds, mesh.comm, verbose, mesh.device);
  free(traceGlobalIds);

  // compute samples of q at interpolation nodes
  acoustics->q = (dfloat*) calloc(Nlocal+Nhalo, sizeof(dfloat));
  acoustics->o_q = mesh.device.malloc((Nlocal+Nhalo)*sizeof(dfloat), acoustics->q);

  //storage for M*q during reporting
  acoustics->o_Mq = mesh.device.malloc((Nlocal+Nhalo)*sizeof(dfloat), acoustics->q);

  // OCCA build stuff
  occa::properties kernelInfo = acoustics->props; //copy base occa properties

  //add boundary data to kernel info
  string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;


  kernelInfo["defines/" "p_Nfields"]= acoustics->Nfields;

  const dfloat p_half = 1./2.;
  kernelInfo["defines/" "p_half"]= p_half;

  int maxNodes = mymax(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = 1024/mesh.Np; // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  int cubMaxNodes = mymax(mesh.Np, (mesh.intNfp*mesh.Nfaces));
  kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
  int cubMaxNodes1 = mymax(mesh.Np, (mesh.intNfp));
  kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

  kernelInfo["defines/" "p_Lambda2"]= Lambda2;

  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  // set kernel name suffix
  char *suffix;
  if(mesh.elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(mesh.elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(mesh.elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(mesh.elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  // kernels from volume file
  sprintf(fileName, DACOUSTICS "/okl/acousticsVolume%s.okl", suffix);
  sprintf(kernelName, "acousticsVolume%s", suffix);

  acoustics->volumeKernel =  buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);
  // kernels from surface file
  sprintf(fileName, DACOUSTICS "/okl/acousticsSurface%s.okl", suffix);
  sprintf(kernelName, "acousticsSurface%s", suffix);

  acoustics->surfaceKernel = buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);

  // mass matrix operator
  sprintf(fileName, LIBP_DIR "/core/okl/MassMatrixOperator%s.okl", suffix);
  sprintf(kernelName, "MassMatrixOperator%s", suffix);

  acoustics->MassMatrixKernel = buildKernel(mesh.device, fileName, kernelName,
                                            kernelInfo, mesh.comm);


  if (mesh.dim==2) {
    sprintf(fileName, DACOUSTICS "/okl/acousticsInitialCondition2D.okl");
    sprintf(kernelName, "acousticsInitialCondition2D");
  } else {
    sprintf(fileName, DACOUSTICS "/okl/acousticsInitialCondition3D.okl");
    sprintf(kernelName, "acousticsInitialCondition3D");
  }

  acoustics->initialConditionKernel = buildKernel(mesh.device, fileName, kernelName,
                                                  kernelInfo, mesh.comm);

  return *acoustics;
}
