#include "ellipticTri2D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(argc!=3 && argc!=4){
    printf("usage 1: ./main meshes/cavityH005.msh N\n");
    printf("usage 2: ./main meshes/cavityH005.msh N BoundaryConditions.h\n");
    exit(-1);
  }

  // int specify polynomial degree
  int N = atoi(argv[2]);

  // solver can be CG or PCG
  // preconditioner can be JACOBI, OAS, NONE, BLOCKJACOBI
  // method can be IPDG
  char *options =
    //strdup("solver=PCG method=IPDG preconditioner=OAS coarse=COARSEGRID,ALMOND");
    //strdup("solver=PCG,FLEXIBLE preconditioner=OAS,PROJECT,GLOBALALMOND,UBERGRID method=IPDG coarse=COARSEGRID");
    strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG,PROJECT preconditioner=FULLALMOND");
    //strdup("solver=PCG,FLEXIBLE method=IPDG preconditioner=BLOCKJACOBI");

  // set up mesh stuff
  mesh2D *mesh = meshSetupTri2D(argv[1], N);
  ogs_t *ogs;
  precon_t *precon;

  // parameter for elliptic problem (-laplacian + lambda)*q = f
  dfloat lambda = 1;
  //dfloat lambda = 0;

  dfloat tau = 2.f*(mesh->N+1)*(mesh->N+1);

  // capture header file
  char *boundaryHeaderFileName;
  if(argc==3)
    boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticTri2D/homogeneous2D.h"); // default
  else
    boundaryHeaderFileName = strdup(argv[3]);

  // set up
  occa::kernelInfo kernelInfo;
  ellipticSetupTri2D(mesh, kernelInfo);

  //add user defined boundary data
  kernelInfo.addIncludeDefine(boundaryHeaderFileName);

  solver_t *solver = ellipticSolveSetupTri2D(mesh, tau, lambda, mesh->EToB, kernelInfo, options);

  iint Nall = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  dfloat *r   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  dfloat *x   = (dfloat*) calloc(Nall,   sizeof(dfloat));

  // load rhs into r
  dfloat *nrhs = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(iint e=0;e<mesh->Nelements;++e){
    dfloat J = mesh->vgeo[e*mesh->Nvgeo+JID];
    for(iint n=0;n<mesh->Np;++n){
      dfloat xn = mesh->x[n+e*mesh->Np];
      dfloat yn = mesh->y[n+e*mesh->Np];
      nrhs[n] = -(2*M_PI*M_PI+lambda)*cos(M_PI*xn)*cos(M_PI*yn);
    }
    for(iint n=0;n<mesh->Np;++n){
      dfloat rhs = 0;
      for(iint m=0;m<mesh->Np;++m){
	      rhs += mesh->MM[n+m*mesh->Np]*nrhs[m];
      }
      iint id = n+e*mesh->Np;

      r[id] = -rhs*J;
      x[id] = 0;
      mesh->q[id] = nrhs[n];
    }
  }
  free(nrhs);

  occa::memory o_r   = mesh->device.malloc(Nall*sizeof(dfloat), r);
  occa::memory o_x   = mesh->device.malloc(Nall*sizeof(dfloat), x);

  ellipticSolveTri2D(solver, lambda, o_r, o_x, options);

  // copy solution from DEVICE to HOST
  o_x.copyTo(mesh->q);

  dfloat maxError = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      iint   id = e*mesh->Np+n;
      dfloat xn = mesh->x[id];
      dfloat yn = mesh->y[id];
      dfloat exact = cos(M_PI*xn)*cos(M_PI*yn);
      dfloat error = fabs(exact-mesh->q[id]);

      maxError = mymax(maxError, error);
    }
  }

  dfloat globalMaxError = 0;
  MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  if(rank==0)
    printf("globalMaxError = %g\n", globalMaxError);

  meshPlotVTU2D(mesh, "foo", 0);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
