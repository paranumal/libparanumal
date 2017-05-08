#include "ellipticQuad2D.h"


int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=3){
    // to run cavity test case with degree N elements
    printf("usage: ./main meshes/cavityQuadH02.msh N\n");
    exit(-1);
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // int specify polynomial degree 
  int N = atoi(argv[2]);

  // solver can be CG or PCG
  // preconditioner can be JACOBI, OAS, NONE
  // method can be CONTINUOUS or IPDG
  // opt: coarse=COARSEGRID with XXT or AMG
  char *options =
    //strdup("solver=PCG,FLEXIBLE preconditioner=OAS method=IPDG,PROJECT coarse=COARSEGRID,GLOBALALMOND,UBERGRID");
    strdup("solver=PCG,FLEXIBLE preconditioner=FULLALMOND method=IPDG,PROJECT");
    //strdup("solver=PCG,FLEXIBLE preconditioner=OAS method=IPDG,PROJECT coarse=COARSEGRID,XXT");
  
  
  // set up mesh stuff
  mesh2D *mesh = meshSetupQuad2D(argv[1], N);
  ogs_t *ogs;
  precon_t *precon;
  
  // parameter for elliptic problem (-laplacian + lambda)*q = f
  dfloat lambda = 1;
  
  // set up
  occa::kernelInfo kernelInfo;
  ellipticSetupQuad2D(mesh, kernelInfo);

  solver_t *solver = ellipticSolveSetupQuad2D(mesh, lambda, kernelInfo, options);

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

      dfloat f = -(2*M_PI*M_PI+lambda)*cos(M_PI*xn)*cos(M_PI*yn);
      
      r[id] = -wJ*f;
      x[id] = 0; // initial guess
    }
  }

  occa::memory o_r   = mesh->device.malloc(Nall*sizeof(dfloat), r);
  occa::memory o_x   = mesh->device.malloc(Nall*sizeof(dfloat), x);

  ellipticSolveQuad2D(solver, lambda, o_r, o_x, options);

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

      mesh->q[id] -= exact;
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
