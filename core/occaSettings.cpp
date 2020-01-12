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

#include "core.hpp"

occaSettings_t::occaSettings_t(MPI_Comm& _comm):
  settings_t(_comm) {

  //settings format
  newSetting("FORMAT",
             "2.0",
             "Setting file version number",
             {"2.0"});

  newSetting("THREAD MODEL",
             "CUDA",
             "OCCA's Parallel execution platform",
             {"Serial", "OpenMP", "CUDA", "HIP", "OpenCL"});

  newSetting("PLATFORM NUMBER",
             "0",
             "Parallel platform number (used in OpenCL mode)");

  newSetting("DEVICE NUMBER",
             "0"
             "Parallel device number");
}

void occaSettings_t::report() {

  int rank;
  MPI_Comm_rank(comm, &rank);

  if (rank==0) {
    std::cout << "OCCA Settings:\n\n";

    reportSetting("THREAD MODEL");

    if (compareSetting("THREAD MODEL","OpenCL"))
      reportSetting("PLATFORM NUMBER");

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if ((size==1)
        &&(compareSetting("THREAD MODEL","CUDA")
        ||compareSetting("THREAD MODEL","HIP")
        ||compareSetting("THREAD MODEL","OpenCL") ))
      reportSetting("DEVICE NUMBER");
  }
}