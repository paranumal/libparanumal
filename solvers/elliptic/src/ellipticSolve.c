#include "elliptic.h"

int ellipticSolve(elliptic_t *elliptic, dfloat lambda, dfloat tol,
                  occa::memory &o_r, occa::memory &o_x){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int Niter = 0;
  int maxIter = 5000; 

  double start = 0.0, end =0.0;

  if(options.compareArgs("VERBOSE","TRUE")){
    mesh->device.finish();
    start = MPI_Wtime(); 
  }

  occaTimerTic(mesh->device,"Linear Solve");
  Niter = pcg (elliptic, lambda, o_r, o_x, tol, maxIter);
  occaTimerToc(mesh->device,"Linear Solve");

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(options.compareArgs("VERBOSE","TRUE")){
    mesh->device.finish();
    end = MPI_Wtime();
    double localElapsed = end-start;

    occa::printTimer();

    if(rank==0) printf("Solver converged in %d iters \n", Niter );

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    hlong   localDofs = (hlong) mesh->Np*mesh->Nelements;
    hlong   localElements = (hlong) mesh->Nelements;
    double globalElapsed;
    hlong   globalDofs;
    hlong   globalElements;

    MPI_Reduce(&localElapsed, &globalElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_HLONG,   MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce(&localElements,&globalElements,1, MPI_HLONG,   MPI_SUM, 0, MPI_COMM_WORLD );

    if (rank==0){
      printf("%02d %02d "hlongFormat" "hlongFormat" %d %17.15lg %3.5g \t [ RANKS N NELEMENTS DOFS ITERATIONS ELAPSEDTIME PRECONMEMORY] \n",
             size, mesh->N, globalElements, globalDofs, Niter, globalElapsed, elliptic->precon->preconBytes/(1E9));
    }
  }
  return Niter;

}
