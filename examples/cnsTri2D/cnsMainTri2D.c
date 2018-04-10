#include "cnsTri2D.h"

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

  // SET OPTIONS
  // integrator = LSERK, DOPRI5
  // out  = REPORT, REPORT+VTU
  // adv  = CUBATURE, COLLOCATION
  char *options = strdup("integrator = DOPRI5, out=VTU, adv=CUBATURE");
  //    char *options = strdup("integrator = LSERK, out=VTU, adv=CUBATURE"); 

  // set up mesh stuff
  mesh2D *mesh = meshSetupTri2D(argv[1], N);

  char *boundaryHeaderFileName = strdup(DHOLMES "/examples/cnsTri2D/cnsUniform2D.h"); // default

  // set up cns stuff
  cns_t *cns = cnsSetupTri2D(mesh, options, boundaryHeaderFileName);

  // run
  cnsRunTri2D(cns, options);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
