#include "mppf.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2){
    printf("usage: ./insMain setupfile\n");
    MPI_Finalize();
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

  mppf_t *mppf = mppfSetup(mesh,options);

  // if(ins->readRestartFile){
  //   printf("Reading restart file..."); 
  //   insRestartRead(ins, ins->options); 
  //   printf("done\n");   
  //  }  
  
  // if (ins->options.compareArgs("TIME INTEGRATOR", "ARK"))  insRunARK(ins);
  // if (ins->options.compareArgs("TIME INTEGRATOR", "EXTBDF"))  insRunEXTBDF(ins);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
