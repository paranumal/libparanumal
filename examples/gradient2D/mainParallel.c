#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"

void message(char *s){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    printf("%s", s);
  MPI_Barrier(MPI_COMM_WORLD);
}


int main(int argc, char **argv){

  MPI_Init(&argc, &argv);

  if(argc!=3){
    printf("usage: ./main meshname.msh N\n");
    exit(-1);
  }

  mesh2D *mesh;

  // read chunk of elements
  mesh = meshParallelReader2D(argv[1]);

  // int specify polynomial degree 
  int N = atoi(argv[2]);
  
  // partition elements using Morton ordering & parallel sort
  meshGeometricPartition2D(mesh);

  // print out connectivity statistics
  meshPartitionStatistics2D(mesh);

  // connect elements using parallel sort
  meshParallelConnect2D(mesh);

  // load reference (r,s) element nodes
  meshLoadReferenceNodes2D(mesh, N);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodes2D(mesh);

  // compute geometric factors
  meshGeometricFactors2D(mesh);

  // compute samples of q at interpolation nodes
  dfloat *q = (dfloat*) calloc(mesh->Nelements*mesh->Np,
			       sizeof(dfloat));
  dfloat *dqdx = (dfloat*) calloc(mesh->Nelements*mesh->Np,
				  sizeof(dfloat));
  dfloat *dqdy = (dfloat*) calloc(mesh->Nelements*mesh->Np,
				  sizeof(dfloat));
  
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      q[cnt] = mesh->y[cnt];
      ++cnt;
    }
  }

  // compute interpolant and collocation derivative of
  // physical x,y derivatives of q data 
  // interpolated at interpolation nodes
  meshGradient2D(mesh, q, dqdx, dqdy);

  // run a test kernel
  //occaTest2D(mesh, q, dqdx, dqdy);

  // run timings for 9 versions of meshGradient kernel 
  occaOptimizeGradient2D(mesh, q, dqdx, dqdy);

  // estimate pointwise error at interpolation nodes
  dfloat maxError = 0;
  cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      maxError = mymax(maxError, fabs(dqdy[cnt]-1));
      ++cnt;
    }
  }

  printf("maxError = %6.5E\n", maxError);

  // find face-node to face-node connectivity
  meshConnectFaceNodes2D(mesh);

  // create single vtu file name
  char vtuName[BUFSIZ];
  sprintf(vtuName, "%s.vtu", strtok(argv[1], "."));
  meshVTU2D(mesh, vtuName);
  
  MPI_Finalize();

  exit(0);
  return 0;
}
