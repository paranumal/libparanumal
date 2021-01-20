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

#include "elliptic.hpp"
#include "mesh/mesh2D.hpp"
#include "mesh/mesh3D.hpp"

void plotGeometry(elliptic_t &elliptic, int num){
  mesh_t &mesh = elliptic.mesh;

  dlong nx, ny, nz;
  dfloat *q = (dfloat*) calloc(elliptic.Nfields*mesh.Np*mesh.Nelements, sizeof(dfloat));

  mesh.settings.getSetting("BOX NX", nx);
  mesh.settings.getSetting("BOX NY", ny);

  // plot mesh
  if(mesh.dim==3){
    mesh.settings.getSetting("BOX NZ", nz);
    
    for(dlong e=0;e<mesh.Nelements;++e){
      dlong ex = e%nx;
      dlong ey = (e/nx)%ny;
      dlong ez = (e/(nx*ny));
      dlong flag = ((ex+ey+ez)%2)==0;
      //    printf("e=%d,flag=%d\n", e, flag);
      for(dlong n=0;n<mesh.Np;++n){
	dlong id = e*mesh.Np*elliptic.Nfields + n;
	q[id] = flag;
      }
    }
  }
  else{
    for(dlong e=0;e<mesh.Nelements;++e){
      dlong ex = e%nx;
      dlong ey = (e/nx);
      dlong flag = ((ex+ey)%2)==0;
      for(dlong n=0;n<mesh.Np;++n){
	dlong id = e*mesh.Np*elliptic.Nfields + n;
	q[id] = flag;
      }
    }
  }
  
  char fileName[BUFSIZ];
  sprintf(fileName, "geometry%05d.vtu", num);
  elliptic.PlotFields(q, fileName);
}



int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  if(argc!=2)
    LIBP_ABORT(string("Usage: ./ellipticMain setupfile"));

  //create default settings
  platformSettings_t platformSettings(comm);
  meshSettings_t meshSettings(comm);
  ellipticSettings_t ellipticSettings(comm);
  ellipticAddRunSettings(ellipticSettings);

  //load settings from file
  ellipticSettings.parseFromFile(platformSettings, meshSettings,
                                 argv[1]);

  // set up platform
  platform_t platform(platformSettings);

  platformSettings.report();
  meshSettings.report();
  ellipticSettings.report();

  // set up mesh
  mesh_t& mesh = mesh_t::Setup(platform, meshSettings, comm);
  
  // map coords
  if(mesh.elementType==HEXAHEDRA)
    ((meshHex3D&)mesh).CoordinateTransform(mesh.N+1, "GL");
  if(mesh.elementType==QUADRILATERALS)
    ((meshQuad2D&)mesh).CoordinateTransform(mesh.N+1, "GL");
  
  dfloat lambda = 0.0;
  ellipticSettings.getSetting("LAMBDA", lambda);

  // Boundary Type translation. Just defaults.
  int NBCTypes = 3;
  int BCType[NBCTypes] = {0,1,2};

  // set up elliptic solver
  elliptic_t& elliptic = elliptic_t::Setup(platform, mesh, ellipticSettings,
                                           lambda, NBCTypes, BCType);

  
  //  plotGeometry(elliptic, 100);
  
  // run
  elliptic.Run();

  
  // close down MPI
  MPI_Finalize();
  return LIBP_SUCCESS;
}
