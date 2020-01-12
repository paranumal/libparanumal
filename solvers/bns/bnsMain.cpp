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

#include "bns.hpp"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  if(argc!=2)
    LIBP_ABORT(string("Usage: ./bnsMain setupfile"));

  //create default settings
  occaSettings_t occaSettings(comm);
  meshSettings_t meshSettings(comm);
  bnsSettings_t bnsSettings(comm);

  //load settings from file
  bnsSettings.parseFromFile(occaSettings, meshSettings,
                            argv[1]);

  occaSettings.report();
  meshSettings.report();
  bnsSettings.report();

  // set up occa device
  occa::device device;
  occa::properties props;
  occaDeviceConfig(device, comm, occaSettings, props);

  // set up mesh
  mesh_t& mesh = mesh_t::Setup(device, comm, meshSettings, props);

  // set up linear algebra module
  linAlg_t& linAlg = linAlg_t::Setup(device, props);

  // set up bns solver
  bns_t& bns = bns_t::Setup(mesh, linAlg, bnsSettings);

  // run
  bns.Run();

  // clean up
  delete &bns;
  delete &linAlg;
  delete &mesh;
  device.free();

  // close down MPI
  MPI_Finalize();
  return LIBP_SUCCESS;
}
