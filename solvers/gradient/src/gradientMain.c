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

    gradient->q[n] = cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*z);
  }

  gradient->o_q.copyFrom(gradient->q);

  // call gradient kernel
  gradient->gradientKernel(mesh->Nelements,
			   mesh->o_vgeo,
			   mesh->o_Dmatrices,
			   gradient->o_q,
			   gradient->o_gradientq);

  // copy gradient back to host
  gradient->o_gradientq.copyTo(gradient->gradientq);

#if 0
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      int id = e*mesh->Np*mesh->dim + n;
      if(mesh->dim==3)
	printf("(%g,%g,%g)\n ",
	       gradient->gradientq[id],
	       gradient->gradientq[id+mesh->Np],
	       gradient->gradientq[id+2*mesh->Np]);
      else
	printf("(%g,%g)\n ",
	       gradient->gradientq[id],
	       gradient->gradientq[id+mesh->Np]);
    }
  }
#endif
  
  //extract isosurface
  if(mesh->dim==3)
   gradientReport(gradient, 0.0,  options);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
