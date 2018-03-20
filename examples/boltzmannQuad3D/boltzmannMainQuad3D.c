#include "boltzmannQuad3D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  //flags will change quickly for now
  /*  if(argc!=3){
    // to run cavity test case with degree N elements
    printf("usage: ./main meshes/cavityH005.msh N\n");
    exit(-1);
    }*/

  // int specify polynomial degree 
  int N = atoi(argv[2]);

  // set up mesh stuff
  dfloat sphereRadius = 1;
  mesh_t *mesh = meshSetupQuad3D(argv[1], N, sphereRadius);

  // set up boltzmann stuff
  solver_t *solver = boltzmannSetupMRQuad3D(mesh,atoi(argv[6]),atoi(argv[7]));

  //load flags into solver
  solver->lserk = atoi(argv[3]); //0 or 1
  solver->mrsaab = atoi(argv[4]); //0 or 1
  solver->filter = atoi(argv[5]); //0 or 1
  //these are also used in setup
  solver->cfl = atoi(argv[6]); //cfl x_
  solver->force_type = atoi(argv[7]); //1->1 2->t 3->original values 
  
  // time step Boltzmann equations
  if (solver->lserk)
    boltzmannRunLSERKQuad3D(solver);
  if (solver->mrsaab)
    boltzmannRunMRSAABQuad3D(solver);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
