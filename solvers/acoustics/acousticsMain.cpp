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

#include "acoustics.hpp"

void meshReport(mesh_t &mesh, const char *s){

  std::cout << "mesh: " << std::string(s) << std::endl;
  std::cout << "mesh:etype " << mesh.elementType << std::endl;
  std::cout << "mesh:Nel   " << mesh.Nelements << std::endl;
  std::cout << "mesh:Np    " << mesh.Np << std::endl;
}


int main(int argc, char **argv){

  // start up MPI
  comm_t::Init(argc, argv);

  LIBP_ABORT("Usage: ./acousticsMain setupfile", argc!=2);

  { /*Scope so everything is destructed before MPI_Finalize */
    comm_t comm(comm_t::world().Dup());

    //create default settings
    platformSettings_t platformSettings(comm);
    meshSettings_t meshSettings(comm);
    acousticsSettings_t acousticsSettings(comm);

    //load settings from file
    acousticsSettings.parseFromFile(platformSettings, meshSettings,
                              argv[1]);

    // set up platform
    platform_t platform(platformSettings);

    platformSettings.report();
    meshSettings.report();
    acousticsSettings.report();

    // set up mesh
    meshSettings.changeSetting("ELEMENT TYPE", std::to_string(Mesh::TRIANGLES));
    mesh_t mesh(platform, meshSettings, comm);
    meshReport(mesh, "read 1");
    
    meshSettings.changeSetting("ELEMENT TYPE", std::to_string(Mesh::TRIANGLES));
    mesh_t meshTri(platform, meshSettings, comm);
    meshReport(meshTri, "Tri");

    meshSettings.changeSetting("ELEMENT TYPE", std::to_string(Mesh::QUADRILATERALS));
    mesh_t meshQuad(platform, meshSettings, comm);
    meshReport(meshTri, "Quad");

    // set up acoustics solver
    acoustics_t acoustics(platform, mesh, acousticsSettings);

    // run
    acoustics.Run();
  }

  // close down MPI
  comm_t::Finalize();
  return LIBP_SUCCESS;
}
