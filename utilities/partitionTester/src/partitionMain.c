#include "partition.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2){
    // to run cavity test case with degree N elements
    printf("usage: ./main setupTri2D.rc \n");
    exit(-1);
  }
  
  // if argv > 2 then should load input data from argv
  setupAide options(argv[1]);
  
  // set up mesh stuff
  string fileName;
  int N, dim, elementType;
  
  options.getArgs("MESH FILE", fileName);
  options.getArgs("POLYNOMIAL DEGREE", N);
  options.getArgs("ELEMENT TYPE", elementType);
  options.getArgs("MESH DIMENSION", dim);
  
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

  // set up
  partitionSetup(mesh);
  
  // plot mesh file
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  char fname[BUFSIZ];
  sprintf(fname, "foo_%05d_00000.vtu", rank);
  
  meshPlotVTU2D(mesh, fname, 0);
  
  // close down MPI
  MPI_Finalize();

  return 0;
}
