#include "bns.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  // Check input
  if(argc!=2){
    printf("usage2: ./bnsMain setupfile\n");
    exit(-1);
  }
  
  // load input data
  setupAide options(argv[1]);
  
  // setup mesh structure
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

  

   bns_t *bns = bnsSetup(mesh,options);
   if(bns->readRestartFile){
    printf("Reading restart file \n"); 
    bnsRestartRead(bns, options);   
   }  


   bnsRun(bns,options);
   
  // close down MPI
  MPI_Finalize();
  exit(0);
  return 0;
}
