#include "ellipticHex3D.h"
 

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=3){
    // to run cavity test case with degree N elements
    printf("usage: ./main meshes/cavityH005.msh N\n");
    exit(-1);
  }

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  // int specify polynomial degree 
  int N = atoi(argv[2]);

  // solver can be CG or PCG
  // preconditioner can be JACOBI, OAS, NONE
  // method can be CONTINUOUS or IPDG
  char *options =
    strdup("solver=PCG method=CONTINUOUS preconditioner=NONE");
  //    strdup("solver=PCG method=IPDG preconditioner=NONE");

  // set up mesh stuff
  mesh3D *mesh = meshSetupHex3D(argv[1], N);
  ogs_t *ogs;
  precon_t *precon;
  
  // parameter for elliptic problem (-laplacian + lambda)*q = f
  dfloat lambda = 0;
  
  // set up
  occa::kernelInfo kernelInfo;
  ellipticSetupHex3D(mesh, kernelInfo);

  solver_t *solver = ellipticSolveSetupHex3D(mesh, lambda, kernelInfo, options);

  iint Nall = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  dfloat *r   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  dfloat *x   = (dfloat*) calloc(Nall,   sizeof(dfloat));

  // load rhs into r
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){

      iint ggid = e*mesh->Np*mesh->Nggeo + n;
      dfloat wJ = mesh->ggeo[ggid+mesh->Np*GWJID];

      iint   id = e*mesh->Np+n;
      dfloat xn = mesh->x[id];
      dfloat yn = mesh->y[id];
      dfloat zn = mesh->z[id];

      dfloat f = -(3*M_PI*M_PI+lambda)*cos(M_PI*xn)*cos(M_PI*yn)*cos(M_PI*zn);
      //dfloat f=1.0;

      r[id] = -wJ*f;

      x[id] = 0; // initial guess
    }
  }

  occa::memory o_r   = mesh->device.malloc(Nall*sizeof(dfloat), r);
  occa::memory o_x   = mesh->device.malloc(Nall*sizeof(dfloat), x);

  // sync processes
  mesh->device.finish();
  MPI_Barrier(MPI_COMM_WORLD);
  
  double tic = MPI_Wtime();
  iint maxIterations = 10;
  
  iint iterations = ellipticSolveHex3D(solver, lambda, o_r, o_x, maxIterations, options);

  mesh->device.finish();
  double toc = MPI_Wtime();

  double localElapsed = toc-tic;
  iint   localDofs = mesh->Np*mesh->Nelements;
  double globalElapsed;
  iint   globalDofs;
  int    root = 0;
  
  MPI_Reduce(&localElapsed, &globalElapsed, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD );
  MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_IINT,   MPI_SUM, root, MPI_COMM_WORLD );
  
  if(rank==root){
    printf("%02d %02d %17.15lg %d %17.15E %17.15E \t [ N DOFS ELAPSEDTIME ITERATIONS (DOFS/RANKS) (DOFS/TIME/ITERATIONS/RANKS)]\n",
	   mesh->N, globalDofs, globalElapsed, iterations, globalDofs/(double)size, globalDofs/(globalElapsed/(iterations/size)));
  }
  // copy solution from DEVICE to HOST
  o_x.copyTo(mesh->q);
  
  dfloat maxError = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      iint   id = e*mesh->Np+n;
      dfloat xn = mesh->x[id];
      dfloat yn = mesh->y[id];
      dfloat zn = mesh->z[id];
      dfloat exact = cos(M_PI*xn)*cos(M_PI*yn)*cos(M_PI*zn);
      dfloat error = fabs(exact-mesh->q[id]);
      
      maxError = mymax(maxError, error);

      //mesh->q[id] -= exact;
    }
  }

  dfloat globalMaxError = 0;
  MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  if(rank==0)
    printf("globalMaxError = %g\n", globalMaxError);
  
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
