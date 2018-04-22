#include "cnsTri2D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2 && argc!=3){
    //    printf("usage1: ./main meshes/cavityH005.msh N\n");
    printf("usage2: ./main setupfile\n");
    exit(-1);
  }

  // int specify polynomial degree 
  int N = atoi(argv[2]);

  // if argv > 2 then should load from argv
  setupAide newOptions("setuprc");
  
  // set up mesh stuff
  mesh2D *mesh = meshSetupTri2D(argv[1], N);

  char *boundaryHeaderFileName = strdup(DHOLMES "/examples/cnsTri2D/cnsUniform2D.h"); // default

  // set up cns stuff
  cns_t *cns = cnsSetupTri2D(mesh, newOptions, boundaryHeaderFileName);

  // run
  cnsRunTri2D(cns, newOptions);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
