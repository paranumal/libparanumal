#include "acoustics.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2){
    printf("usage2: ./acousticsMain setupfile\n");
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
    mesh = meshSetupTri2D((char*)fileName.c_str(), N); break;
  case QUADRILATERALS:
    mesh = meshSetupQuad2D((char*)fileName.c_str(), N); break;
  case TETRAHEDRA:
    mesh = meshSetupTet3D((char*)fileName.c_str(), N); break;
  case HEXAHEDRA:
    mesh = meshSetupHex3D((char*)fileName.c_str(), N); break;
  }

  char *boundaryHeaderFileName; // could sprintf
  if(dim==2)
    boundaryHeaderFileName = strdup(DACOUSTICS "/acousticsUniform2D.h"); // default
  if(dim==3)
    boundaryHeaderFileName = strdup(DACOUSTICS "/acousticsUniform3D.h"); // default

  // set up acoustics stuff
  acoustics_t *acoustics = acousticsSetup(mesh, newOptions, boundaryHeaderFileName);

  // run
  acousticsRun(acoustics, newOptions);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
