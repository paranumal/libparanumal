#include "ellipticHex3D.h"


void timeAxOperator(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const char *options){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh_t *mesh = solver->mesh;

  // sync processes
  mesh->device.finish();
  MPI_Barrier(MPI_COMM_WORLD);

  double tic = MPI_Wtime();
  double AxTime;

  iint iterations = 10;

  occa::streamTag start = mesh->device.tagStream();


#if 1

  void ellipticOperator3D(solver_t *solver, dfloat lambda,
			  occa::memory &o_q, occa::memory &o_Aq, const char *options);

    // assume 1 mpi process
  for(int it=0;it<iterations;++it){

    ellipticOperator3D(solver, lambda, o_r, o_x, options);
  }
#else
  // assume 1 mpi process
  for(int it=0;it<iterations;++it)
    solver->partialAxKernel(solver->NlocalGatherElements,
			    solver->o_localGatherElementList,
			    solver->o_gjGeo,
			    solver->o_gjD,
			    solver->o_gjI,
			    lambda, o_r,
			    solver->o_grad,
			    o_x);
#endif

  occa::streamTag end = mesh->device.tagStream();

  mesh->device.finish();
  double toc = MPI_Wtime();

  double localElapsed = toc-tic;
  //  localElapsed = mesh->device.timeBetween(start, end);

  iint   localDofs = mesh->Np*mesh->Nelements;
  iint localElements = mesh->Nelements;
  double globalElapsed;
  iint   globalDofs;
  iint   globalElements;
  int    root = 0;

  MPI_Reduce(&localElapsed, &globalElapsed, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD );
  MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_IINT,   MPI_SUM, root, MPI_COMM_WORLD );
  MPI_Reduce(&localElements,&globalElements,1, MPI_IINT,   MPI_SUM, root, MPI_COMM_WORLD );

  iint gjNq = mesh->gjNq;
  iint Nq = mesh->Nq;

  double flops = gjNq*Nq*Nq*Nq*4 +
    gjNq*gjNq*Nq*Nq*6 +
    gjNq*gjNq*gjNq*Nq*8 +
    gjNq*gjNq*gjNq*17 +
    gjNq*gjNq*gjNq*Nq*8 +
    gjNq*gjNq*Nq*Nq*6 +
    gjNq*Nq*Nq*Nq*4; // excludes inner product

  double gflops = globalElements*flops*iterations/(1024*1024*1024.*globalElapsed);

  if(rank==root){
    printf("%02d %02d %02d %17.15lg %d %17.15E %17.15E %17.15E \t [ RANKS N DOFS ELAPSEDTIME ITERATIONS (DOFS/RANKS) (DOFS/TIME/ITERATIONS/RANKS) (Ax GFLOPS)]\n",
	   size, mesh->N, globalDofs, globalElapsed, iterations, globalDofs/(double)size, (globalDofs*iterations)/(globalElapsed*size), gflops);
  }


}

void timeSolver(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const char *options){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh_t *mesh = solver->mesh;

  // sync processes
  mesh->device.finish();
  MPI_Barrier(MPI_COMM_WORLD);

  double tic = MPI_Wtime();
  iint maxIterations = 30;
  double AxTime;

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
    printf("%02d %02d %02d %17.15lg %d %17.15E %17.15E \t [ RANKS N DOFS ELAPSEDTIME ITERATIONS (DOFS/RANKS) (DOFS/TIME/ITERATIONS/RANKS) \n",
	   size, mesh->N, globalDofs, globalElapsed, iterations, globalDofs/(double)size, (globalDofs*(double)iterations)/(globalElapsed*size));
  }


}




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

  char customCache[BUFSIZ];
  sprintf(customCache, "/home/tcew/._occa_cache_rank_%05d", rank);
  setenv("OCCA_CACHE_DIR", customCache, 1);

  char *check = getenv("OCCA_CACHE_DIR");
  printf("found OCD: %s\n", check);

  // int specify polynomial degree
  int N = atoi(argv[2]);

  // solver can be CG or PCG
  // preconditioner can be JACOBI, OAS, NONE
  // method can be CONTINUOUS or IPDG
  char *options =
    strdup("solver=CG method=CONTINUOUS preconditioner=NONE");
    //strdup("solver=CG method=IPDG preconditioner=NONE");

  // set up mesh stuff
  mesh3D *mesh = meshSetupHex3D(argv[1], N);
  ogs_t *ogs;
  precon_t *precon;

  // parameter for elliptic problem (-laplacian + lambda)*q = f
  dfloat lambda = 1;

  // set up
  occa::kernelInfo kernelInfo;
  ellipticSetupHex3D(mesh, kernelInfo);

  solver_t *solver = ellipticSolveSetupHex3D(mesh, lambda, kernelInfo, options);

  iint Nall = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  printf("Nall = %d, mesh->Np = %d mesh->Nelements = %d mesh->totalHaloPairs = %d \n", Nall, mesh->Np, mesh->Nelements, mesh->totalHaloPairs);
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

  timeAxOperator(solver, lambda, o_r, o_x, options);

  //  timeSolver(solver, lambda, o_r, o_x, options);

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
    printf("globalMaxError = %17.15g\n", globalMaxError);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
