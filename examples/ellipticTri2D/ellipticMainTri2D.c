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
  // can add FLEXIBLE and VERBOSE options
  // preconditioner can be NONE, JACOBI, OAS, BLOCKJACOBI, FULLALMOND, or MULTIGRID
  // OAS and MULTIGRID: smoothers can be EXACTFULLPATCH, APPROXFULLPATCH, OVERLAPPINGPATCH, or DAMPEDJACOBI
  // FULLALMOND: can include MATRIXFREE option
  // method can be IPDG or CONTINUOUS
  char *options =
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=OAS smoother=EXACTFULLPATCH");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=OAS smoother=APPROXFULLPATCH");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=OAS smoother=OVERLAPPINGPATCH");
    strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=MULTIGRID,APPROXFULLPATCH");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=FULLALMOND,MATRIXFREE");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=NONE");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=BLOCKJACOBI");

  //FULLALMOND, OAS, and MULTIGRID will use the parAlmondOptions in setup
  // solver can be EXACT, KCYCLE, or VCYCLE
  // can add GATHER to build a gsop
  // partition can be STRONGNODES, DISTRIBUTED, SATURATE
  char *parAlmondOptions =
    strdup("solver=KCYCLE,VERBOSE partition=STRONGNODES");
    //strdup("solver=EXACT,VERBOSE partition=STRONGNODES");

  // set up mesh stuff
  mesh2D *mesh = meshSetupTri2D(argv[1], N);
  ogs_t *ogs;
  precon_t *precon;

  // parameter for elliptic problem (-laplacian + lambda)*q = f
  dfloat lambda = 1;
  //dfloat lambda = 0;

  // set up
  occa::kernelInfo kernelInfo;
  ellipticSetupTri2D(mesh, kernelInfo);

  // capture header file
  char *boundaryHeaderFileName;
  if(argc==3)
    boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticTri2D/homogeneous2D.h"); // default
  else
    boundaryHeaderFileName = strdup(argv[3]);
  //add user defined boundary data
  kernelInfo.addInclude(boundaryHeaderFileName);

  //add standard boundary functions
  boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticTri2D/ellipticBoundary2D.h");
  kernelInfo.addInclude(boundaryHeaderFileName);

  // Boundary Type translation. Just default from the mesh file.
  int BCType[3] = {0,1,2};

  dfloat tau = 2.0*(mesh->N+1)*(mesh->N+2)/2.0;
  solver_t *solver = ellipticSolveSetupTri2D(mesh, tau, lambda, BCType, kernelInfo, options, parAlmondOptions);

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
      nrhs[n] = -(2*M_PI*M_PI+lambda)*sin(M_PI*xn)*sin(M_PI*yn);
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

  //add boundary condition contribution to rhs
  if (strstr(options,"IPDG")) {
    dfloat zero = 0.f;
    solver->rhsBCIpdgKernel(mesh->Nelements,
                           mesh->o_vmapM,
                           mesh->o_vmapP,
                           solver->tau,
                           zero,
                           mesh->o_x,
                           mesh->o_y,
                           mesh->o_vgeo,
                           mesh->o_sgeo,
                           mesh->o_EToB,
                           mesh->o_DrT,
                           mesh->o_DsT,
                           mesh->o_LIFTT,
                           mesh->o_MM,
                           o_r);
  }

  ellipticSolveTri2D(solver, lambda, o_r, o_x, options);

  // copy solution from DEVICE to HOST
  o_x.copyTo(mesh->q);

  dfloat maxError = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      iint   id = e*mesh->Np+n;
      dfloat xn = mesh->x[id];
      dfloat yn = mesh->y[id];
      dfloat exact = sin(M_PI*xn)*sin(M_PI*yn);
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
