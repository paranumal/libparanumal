#include "ellipticTri2D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=3){
    // to run cavity test case with degree N elements
    printf("usage: ./main meshes/cavityTriH02.msh N\n");
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
  //  char *options = strdup("solver=PCG preconditioner=OAS,PROJECT method=IPDG coarse=COARSEGRID");
  //char *options = strdup("solver=PCG,FLEXIBLE preconditioner=OAS,PROJECT,GLOBALALMOND,UBERGRID method=IPDG coarse=COARSEGRID");
  char *options = strdup("solver=PCG,FLEXIBLE preconditioner=FULLALMOND method=CONTINUOUS");
  //char *options = strdup("solver=PCG preconditioner=NONE method=IPDG");
  
  // set up mesh stuff
  mesh2D *mesh = meshSetupTri2D(argv[1], N);
  ogs_t *ogs;
  precon_t *precon;
  
  // parameter for elliptic problem (-laplacian + lambda)*q = f
  dfloat lambda = 1;
  
  // set up
  occa::kernelInfo kernelInfo;
  ellipticSetupTri2D(mesh, kernelInfo);

  solver_t *solver = ellipticSolveSetupTri2D(mesh, lambda, kernelInfo, options);

  iint Nall = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  dfloat *r   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  dfloat *x   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  
  // load rhs into r
  dfloat *cf = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *nrhs = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(iint e=0;e<mesh->Nelements;++e){

#if 0
    for(iint n=0;n<mesh->cubNp;++n){
      dfloat cx = 0, cy = 0;
      for(iint m=0;m<mesh->Np;++m){
	cx += mesh->cubInterp[m+n*mesh->Np]*mesh->x[m+e*mesh->Np];
	cy += mesh->cubInterp[m+n*mesh->Np]*mesh->y[m+e*mesh->Np];
      }
      dfloat J = mesh->vgeo[e*mesh->Nvgeo+JID];
      dfloat w = mesh->cubw[n];
      
      cf[n] = -J*w*(2*M_PI*M_PI+lambda)*cos(M_PI*cx)*cos(M_PI*cy);
    }
    for(iint n=0;n<mesh->Np;++n){
      dfloat rhs = 0;
      for(iint m=0;m<mesh->cubNp;++m){
	rhs += mesh->cubInterp[n+m*mesh->Np]*cf[m];
      }
      iint id = n+e*mesh->Np;
      r[id] = -rhs;
      x[id] = 0; // initial guess
    }
#else
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
#endif
  }
  free(nrhs);
  free(cf);

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
