#include "ellipticTet3D.h"

void diagnostic(int N, occa::memory &o_x, const char *message){
#if 1
  dfloat *x = (dfloat*) calloc(N, sizeof(dfloat));

  o_x.copyTo(x, N*sizeof(dfloat), 0);
  
  int n;
  dfloat normX = 0;
  for(n=0;n<N;++n)
    normX += x[n]*x[n];

  dfloat globalNormX;
  MPI_Allreduce(&normX, &globalNormX, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0)
    printf("rank %d reports norm %17.15lf for %s\n", rank, sqrt(globalNormX), message);

  free(x);
#endif
}

