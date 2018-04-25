#include "cnsQuad2D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2){
    printf("usage2: ./cnsMainQuad2D setupfile\n");
    exit(-1);
  }

  // if argv > 2 then should load input data from argv
  setupAide newOptions(argv[1]);

  // set up mesh stuff
  string fileName;
  int N;

  newOptions.getArgs("MESH FILE", fileName);
  newOptions.getArgs("POLYNOMIAL DEGREE", N);
  mesh2D *mesh = meshSetupQuad2D((char*)fileName.c_str(), N);

  char *boundaryHeaderFileName = strdup(DHOLMES "/examples/cnsQuad2D/cnsUniform2D.h"); // default

  // set up cns stuff
  cns_t *cns = cnsSetupQuad2D(mesh, newOptions, boundaryHeaderFileName);

  // run
  cnsRunQuad2D(cns, newOptions);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
