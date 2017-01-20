#include "ellipticHex3D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=3){
    // to run cavity test case with degree N elements
    printf("usage: ./main meshes/cavityH005.msh N\n");
    exit(-1);
  }

  // int specify polynomial degree 
  int N = atoi(argv[2]);

  // set up mesh stuff
  mesh3D *meshSetupHex3D(char *, iint);
  mesh3D *mesh = meshSetupHex3D(argv[1], N);

  // set up elliptic stuff
  void ellipticSetupHex3D(mesh3D *mesh);
  ellipticSetupHex3D(mesh);

  // at this point gather-scatter is available
  dfloat lambda = 10;
  mesh->AxKernel(mesh->Nelements, mesh->o_ggeo, mesh->o_D, lambda, mesh->o_q, mesh->o_rhsq); // store A*q in rhsq

  for(iint n=0;n<mesh->Nelements*mesh->Np;++n)
    mesh->rhsq[n] = 1;

  mesh->o_rhsq.copyFrom(mesh->rhsq);
  
  // do parallel gather scatter
  void meshParallelGatherScatter3D(mesh3D *mesh, occa::memory &o_v, occa::memory &o_gsv, const char *type);
  meshParallelGatherScatter3D(mesh, mesh->o_rhsq, mesh->o_rhsq, "float");

  mesh->o_rhsq.copyTo(mesh->rhsq);

  for(iint n=0;n<mesh->Nelements*mesh->Np;++n)
    printf("deg[%d] = %g\n", n, mesh->rhsq[n]);
  
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
