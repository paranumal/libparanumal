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

#include "wave.hpp"

int main(int argc, char **argv){

  // start up MPI
  Comm::Init(argc, argv);

  LIBP_ABORT("Usage: ./waveMain setupfile", argc!=2);

  { /*Scope so everything is destructed before MPI_Finalize */
    comm_t comm(Comm::World().Dup());

    //create default settings
    platformSettings_t platformSettings(comm);
    meshSettings_t meshSettings(comm);
    waveSettings_t waveSettings(comm);

    // set up platform
    platform_t platform(platformSettings);

    platformSettings.report();
    meshSettings.report();

    // set up mesh
    mesh_t mesh(platform, meshSettings, comm);

    // Boundary Type translation. Just defaults.
    int NBCTypes = 3;
    memory<int> BCType(3);
    BCType[0] = 0;
    BCType[1] = 1;
    BCType[2] = 2;

    wave_t wave(platform, mesh, waveSettings);

    // run
    wave.Run();
  }

  // close down MPI
  Comm::Finalize();
  return LIBP_SUCCESS;
}
