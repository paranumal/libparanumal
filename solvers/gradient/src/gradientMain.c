#include "gradient.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2){
    printf("usage: ./gradientMain setupfile\n");
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

  // set up gradient stuff
  gradient_t *gradient = gradientSetup(mesh, options);

  // populate q
  for(int n=0;n<mesh->Nelements*mesh->Np;++n){
    dfloat x = mesh->x[n];
    dfloat y = mesh->y[n];
    dfloat z = mesh->z[n];

    gradient->q[n] = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
  }

  gradient->o_q.copyFrom(gradient->q);
  
  // extract isosurface
  gradientReport(gradient, 0.0,  options);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
