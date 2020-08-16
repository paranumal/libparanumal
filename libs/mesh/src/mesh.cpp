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

#include "mesh.hpp"
#include "mesh2D.hpp"
#include "mesh3D.hpp"

//makeing a mesh object requires it to be bound to a device and communicator
mesh_t::mesh_t(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  platform(_platform), device(platform.device), settings(_settings), comm(_comm) {
  props = platform.props;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
}

mesh2D::mesh2D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh_t(_platform, _settings, _comm) {}

mesh3D::mesh3D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh_t(_platform, _settings, _comm) {}

meshTri2D::meshTri2D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh2D(_platform, _settings, _comm) {}

meshQuad2D::meshQuad2D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh2D(_platform, _settings, _comm) {}

meshTri3D::meshTri3D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh3D(_platform, _settings, _comm) {}

meshQuad3D::meshQuad3D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh3D(_platform, _settings, _comm) {}

meshTet3D::meshTet3D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh3D(_platform, _settings, _comm) {}

meshHex3D::meshHex3D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh3D(_platform, _settings, _comm) {}

mesh_t::~mesh_t() {
  if (halo) halo->Free();
  if (ringHalo) ringHalo->Free();
  if (ogs) ogs->Free();
}