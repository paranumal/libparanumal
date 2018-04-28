#include "cnsHex3D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2){
    printf("usage2: ./cnsMainHex3D setupfile\n");
    exit(-1);
  }

  // if argv > 2 then should load input data from argv
  setupAide newOptions(argv[1]);

  // set up mesh stuff
  string fileName;
  int N;

  newOptions.getArgs("MESH FILE", fileName);
  newOptions.getArgs("POLYNOMIAL DEGREE", N);
  mesh3D *mesh = meshSetupHex3D((char*)fileName.c_str(), N);

  char *boundaryHeaderFileName = strdup(DHOLMES "/examples/cnsHex3D/cnsUniform3D.h"); // default

  // set up cns stuff
  cns_t *cns = cnsSetupHex3D(mesh, newOptions, boundaryHeaderFileName);

  // run
  cnsRunHex3D(cns, newOptions);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
