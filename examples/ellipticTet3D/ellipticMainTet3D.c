#include "ellipticTet3D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=3){
    // to run cavity test case with degree N elements
    printf("usage: ./main meshes/cavityTetH02.msh N\n");
    exit(-1);
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // int specify polynomial degree 
  int N = atoi(argv[2]);

  // solver can be CG or PCG
  // preconditioner can be JACOBI, OAS, NONE
  // method can be IPDG
  //char *options = strdup("solver=PCG preconditioner=OAS method=IPDG");
  //char *options = strdup("solver=PCG preconditioner=OAS,PROJECT method=IPDG coarse=NONE");
  //char *options = strdup("solver=PCG,FLEXIBLE preconditioner=OAS,GLOBALALMOND,UBERGRID method=IPDG,PROJECT coarse=COARSEGRID");
  char *options = strdup("solver=PCG,FLEXIBLE preconditioner=FULLALMOND,UBERGRID,MATRIXFREE method=IPDG,PROJECT");
  //char *options = strdup("solver=PCG preconditioner=FULLALMOND method=IPDG,PROJECT");
  
  // set up mesh stuff

  mesh3D *mesh = meshSetupTet3D(argv[1], N);
  ogs_t *ogs;
  precon_t *precon;
  
  // parameter for elliptic problem (-laplacian + lambda)*q = f
  dfloat lambda = 1;
  
  // set up
  //  ellipticSetupTet3D(mesh, &ogs, &precon, lambda);
  occa::kernelInfo kernelInfo;
  ellipticSetupTet3D(mesh, kernelInfo);

  solver_t *solver = ellipticSolveSetupTet3D(mesh, lambda, kernelInfo, options);

  iint Nall = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  dfloat *r   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  dfloat *x   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  
  // convergence tolerance (currently absolute)
  const dfloat tol = 1e-6;

  // load rhs into r
  dfloat *nrhs = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(iint e=0;e<mesh->Nelements;++e){
    dfloat J = mesh->vgeo[e*mesh->Nvgeo+JID];
    for(iint n=0;n<mesh->Np;++n){
      dfloat xn = mesh->x[n+e*mesh->Np];
      dfloat yn = mesh->y[n+e*mesh->Np];
      dfloat zn = mesh->z[n+e*mesh->Np];
      //x[n+e*mesh->Np] = cos(M_PI*xn)*cos(M_PI*yn)*cos(M_PI*zn);
      nrhs[n] = -(3*M_PI*M_PI+lambda)*cos(M_PI*xn)*cos(M_PI*yn)*cos(M_PI*zn);
    }
    for(iint n=0;n<mesh->Np;++n){
      dfloat rhs = 0;
      for(iint m=0;m<mesh->Np;++m){
	      rhs += mesh->MM[n+m*mesh->Np]*nrhs[m];
      }
      iint id = n+e*mesh->Np;
      
      r[id] = -rhs*J;
      x[id] = 0.;
      mesh->q[id] = nrhs[n];
    }
  }
  free(nrhs);

  occa::memory o_r   = mesh->device.malloc(Nall*sizeof(dfloat), r);
  occa::memory o_x   = mesh->device.malloc(Nall*sizeof(dfloat), x);

  ellipticSolveTet3D(solver, lambda, o_r, o_x, options);
  
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
    }
  }
  
  dfloat globalMaxError = 0;
  MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  if(rank==0)
    printf("globalMaxError = %g\n", globalMaxError);
  
  meshPlotVTU3D(mesh, "foo", 0);
  
  // close down MPI
  MPI_Finalize();
  
  exit(0);
  return 0;
}
