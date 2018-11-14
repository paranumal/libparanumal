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

#include "advection.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2){
    printf("usage2: ./advectionMain setupfile\n");
    exit(-1);
  }

  // if argv > 2 then should load input data from argv
  setupAide newOptions(argv[1]);
  
  // set up mesh stuff
  string fileName;
  int N, dim, elementType;

  newOptions.getArgs("MESH FILE", fileName);
  newOptions.getArgs("POLYNOMIAL DEGREE", N);
  newOptions.getArgs("ELEMENT TYPE", elementType);
  newOptions.getArgs("MESH DIMENSION", dim);
  
  // set up mesh
  mesh_t *mesh;
  switch(elementType){
  case TRIANGLES:
  case TETRAHEDRA:
    printf("Triangles and tetrahedra are not currently supported for this code, exiting ...\n");
    exit(-1);
  case QUADRILATERALS:
    mesh = meshSetupQuad2D((char*)fileName.c_str(), N); break;
  case HEXAHEDRA:
    mesh = meshSetupHex3D((char*)fileName.c_str(), N); break;
  }

  if(elementType==HEXAHEDRA){
    
    /* rescale to unit box and transform */
    hlong allNelements = mesh->Nelements+mesh->totalHaloPairs;
    for(int n=0;n<allNelements*mesh->Np;++n){
      mesh->x[n] = 0.5*(mesh->x[n]+1);
      mesh->y[n] = 0.5*(mesh->y[n]+1);
      mesh->z[n] = 0.5*(mesh->z[n]+1);
    }
    
    // compute geometric factors
    meshGeometricFactorsHex3D(mesh);
    meshSurfaceGeometricFactorsHex3D(mesh);
  }
  
  char *boundaryHeaderFileName; // could sprintf
  if(dim==2)
    boundaryHeaderFileName = strdup(DADVECTION "/advectionBox2D.h"); // default
  if(dim==3)
    boundaryHeaderFileName = strdup(DADVECTION "/advectionBox3D.h"); // default

  // set up advection stuff
  advection_t *advection = advectionSetup(mesh, newOptions, boundaryHeaderFileName);

  // run
  advectionRun(advection, newOptions);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
