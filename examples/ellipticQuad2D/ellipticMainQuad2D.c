#include "ellipticQuad2D.h"


int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(argc<3){
    printf("usage 1: ./main meshes/cavityH005.msh N\n");
    printf("usage 2: ./main meshes/cavityH005.msh N BoundaryConditions.h\n");
    exit(-1);
  }

  // int specify polynomial degree
  int N = atoi(argv[2]);

  // solver can be CG or PCG
  // can add FLEXIBLE and VERBOSE options
  // method can be IPDG or CONTINUOUS
  // preconditioner can be NONE, JACOBI, OAS, MASSMATRIX, FULLALMOND, or MULTIGRID
  // OAS and MULTIGRID: smoothers can be FULLPATCH, FACEPATCH, LOCALPATCH, OVERLAPPINGPATCH, or DAMPEDJACOBI
  //                      patch smoothers can include EXACT        
  // MULTIGRID: smoothers can include CHEBYSHEV for smoother acceleration
  // MULTIGRID: levels can be ALLDEGREES, HALFDEGREES, HALFDOFS
  // FULLALMOND: can include MATRIXFREE option
  char *options =
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=OAS smoother=FULLPATCH");
    strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=FULLPATCH,EXACT");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=FULLALMOND");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=NONE");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=JACOBI");

  //FULLALMOND, OAS, and MULTIGRID will use the parAlmondOptions in setup
  // solver can be EXACT, KCYCLE, or VCYCLE
  // smoother can be DAMPEDJACOBI or CHEBYSHEV
  // can add GATHER to build a gsop
  // partition can be STRONGNODES, DISTRIBUTED, SATURATE
  char *parAlmondOptions =
    strdup("solver=KCYCLE smoother=CHEBYSHEV partition=STRONGNODES");
    //strdup("solver=EXACT,VERBOSE smoother=DAMPEDJACOBI partition=STRONGNODES");

  //this is strictly for testing, to do repeated runs. Will be removed later
  if (argc==6) {
    options = strdup(argv[4]);
    parAlmondOptions = strdup(argv[5]);
  }

  // set up mesh stuff
  mesh2D *mesh = meshSetupQuad2D(argv[1], N);
  ogs_t *ogs;
  precon_t *precon;

  // parameter for elliptic problem (-laplacian + lambda)*q = f
  dfloat lambda = 1;

  // set up
  occa::kernelInfo kernelInfo;
  ellipticSetupQuad2D(mesh, kernelInfo);

  // capture header file
  char *boundaryHeaderFileName;
  if(argc==3)
    boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticQuad2D/homogeneous2D.h"); // default
  else
    boundaryHeaderFileName = strdup(argv[3]);
  //add user defined boundary data
  kernelInfo.addInclude(boundaryHeaderFileName);

  //add standard boundary functions
  boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticQuad2D/ellipticBoundary2D.h");
  kernelInfo.addInclude(boundaryHeaderFileName);

  // Boundary Type translation. Just default from the mesh file.
  int BCType[3] = {0,1,2};

  dfloat tau = (mesh->N+1)*(mesh->N+1);
  solver_t *solver = ellipticSolveSetupQuad2D(mesh, tau, lambda, BCType, kernelInfo, options);

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

      dfloat f = -(2*M_PI*M_PI+lambda)*sin(M_PI*xn)*sin(M_PI*yn);
      //dfloat f = 1.0;

      r[id] = -wJ*f;
      x[id] = 0; // initial guess
    }
  }

  occa::memory o_r   = mesh->device.malloc(Nall*sizeof(dfloat), r);
  occa::memory o_x   = mesh->device.malloc(Nall*sizeof(dfloat), x);

  //add boundary condition contribution to rhs
  if (strstr(options,"IPDG")) {
    dfloat zero = 0.f;
    solver->rhsBCIpdgKernel(mesh->Nelements,
                           mesh->o_vmapM,
                           solver->tau,
                           zero,
                           mesh->o_x,
                           mesh->o_y,
                           mesh->o_vgeo,
                           mesh->o_sgeo,
                           mesh->o_EToB,
                           mesh->o_D,
                           o_r);
  }

  ellipticSolveQuad2D(solver, lambda, o_r, o_x, options);

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

  meshPlotVTU2D(mesh, "foo.vtu", 0);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
