#include "ellipticTet3D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc<3){
    // to run cavity test case with degree N elements
    printf("usage 1: ./main meshes/cavityTetH02.msh N\n");
    printf("usage 2: ./main meshes/cavityTetH02.msh N BoundaryConditions.h\n");
    exit(-1);
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // int specify polynomial degree
  int N = atoi(argv[2]);

  // solver can be CG or PCG
  // can add FLEXIBLE and VERBOSE options
  // method can be IPDG or CONTINUOUS
  // preconditioner can be NONE, JACOBI, MASSMATRIX, FULLALMOND, or MULTIGRID
  // MULTIGRID: smoothers can be FULLPATCH, FACEPATCH, LOCALPATCH, or DAMPEDJACOBI, patch smoothers can include EXACT
  // MULTIGRID: smoothers can include CHEBYSHEV for smoother acceleration
  // MULTIGRID: levels can be ALLDEGREES, HALFDEGREES, HALFDOFS
  // FULLALMOND: can include MATRIXFREE option
  char *options =
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=DAMPEDJACOBI,CHEBYSHEV");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=FULLALMOND");
    strdup("solver=PCG,FLEXIBLE,VERBOSE method=CONTINUOUS preconditioner=NONE");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=JACOBI");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG preconditioner=MASSMATRIX");

  //FULLALMOND, OAS, and MULTIGRID will use the parAlmondOptions in setup
  // solver can be EXACT, KCYCLE, or VCYCLE
  // smoother can be DAMPEDJACOBI or CHEBYSHEV
  // can add GATHER to build a gsop
  // partition can be STRONGNODES, DISTRIBUTED, SATURATE
  char *parAlmondOptions =
    strdup("solver=KCYCLE,VERBOSE smoother=CHEBYSHEV partition=STRONGNODES");
    //strdup("solver=EXACT,VERBOSE smoother=DAMPEDJACOBI partition=STRONGNODES");


  //this is strictly for testing, to do repeated runs. Will be removed later
  if (argc==6) {
    options = strdup(argv[4]);
    parAlmondOptions = strdup(argv[5]);
  }

  mesh3D *mesh = meshSetupTet3D(argv[1], N);
  ogs_t *ogs;
  precon_t *precon;

  // parameter for elliptic problem (-laplacian + lambda)*q = f
  dfloat lambda = 1;

  // set up
  occa::kernelInfo kernelInfo;
  ellipticSetupTet3D(mesh, kernelInfo);

  // Boundary Type translation. Just default from the mesh file.
  int BCType[3] = {0,1,2};

  dfloat tau = 2.0*(mesh->N+1)*(mesh->N+3);
  solver_t *solver = ellipticSolveSetupTet3D(mesh, tau, lambda, BCType, kernelInfo, options, parAlmondOptions);

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
      dfloat zn = mesh->z[n+e*mesh->Np];

      nrhs[n] = -(3*M_PI*M_PI+lambda)*sin(M_PI*xn)*sin(M_PI*yn)*sin(M_PI*zn);
      x[e*mesh->Np+n] = sin(M_PI*xn)*sin(M_PI*yn)*sin(M_PI*zn);
    }
    for(iint n=0;n<mesh->Np;++n){
      dfloat rhs = 0;
      for(iint m=0;m<mesh->Np;++m){
	      rhs += mesh->MM[n+m*mesh->Np]*nrhs[m];
      }
      iint id = n+e*mesh->Np;

      r[id] = -rhs*J;
      //x[id] = 0.;

      mesh->q[id] = nrhs[n];
    }
  }
  free(nrhs);

  occa::memory o_r   = mesh->device.malloc(Nall*sizeof(dfloat), r);
  occa::memory o_x   = mesh->device.malloc(Nall*sizeof(dfloat), x);

  // capture header file
  char *boundaryHeaderFileName;
  if(argc==3)
    boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticTet3D/homogeneous3D.h"); // default
  else
    boundaryHeaderFileName = strdup(argv[3]);
  //add user defined boundary data
  kernelInfo.addInclude(boundaryHeaderFileName);

  //add boundary condition contribution to rhs
  if (strstr(options,"IPDG")) {
    solver->rhsBCIpdgKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticRhsBCIpdgTet3D.okl",
               "ellipticRhsBCIpdgTet3D",
               kernelInfo);

    dfloat zero = 0.f;
    solver->rhsBCIpdgKernel(mesh->Nelements,
                           mesh->o_vmapM,
                           mesh->o_vmapP,
                           solver->tau,
                           zero,
                           mesh->o_x,
                           mesh->o_y,
                           mesh->o_z,
                           mesh->o_vgeo,
                           mesh->o_sgeo,
                           mesh->o_EToB,
                           mesh->o_DrT,
                           mesh->o_DsT,
                           mesh->o_DtT,
                           mesh->o_LIFTT,
                           mesh->o_MM,
                           o_r);
  }

  if (strstr(options,"CONTINUOUS")) {
    solver->rhsBCKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticRhsBCTet3D.okl",
          "ellipticRhsBCTet3D",
          kernelInfo);

    solver->addBCKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAddBCTet3D.okl",
          "ellipticAddBCTet3D",
          kernelInfo);

    dfloat zero = 0.f;
    solver->rhsBCKernel(mesh->Nelements,
                        mesh->o_ggeo,
                        mesh->o_sgeo,
                        mesh->o_SrrT, mesh->o_SrsT, mesh->o_SrtT,
                        mesh->o_SsrT, mesh->o_SssT, mesh->o_SstT,
                        mesh->o_StrT, mesh->o_StsT, mesh->o_SttT,
                        mesh->o_MM,
                        mesh->o_vmapM,
                        mesh->o_sMT,
                        lambda,
                        zero,
                        mesh->o_x,
                        mesh->o_y,
                        mesh->o_z,
                        mesh->o_mapB,
                        o_r);
  }

  // gather-scatter
  if(strstr(options, "CONTINUOUS")){
    ellipticParallelGatherScatterTet3D(mesh, mesh->ogs, o_r, o_r, dfloatString, "add");  
    if (mesh->Nmasked) mesh->maskKernel(mesh->Nmasked, mesh->o_maskIds, o_r);
  }

  // convergence tolerance
  dfloat tol = 1e-8;
  ellipticSolveTet3D(solver, lambda, tol, o_r, o_x, options);
  //o_x.copyFrom(o_r);
  //o_r.copyTo(r);
  //ellipticOperator3D(solver, lambda, o_x, solver->o_Ax, options);
  //o_x.copyFrom(solver->o_Ax);


  if(strstr(options, "CONTINUOUS")){
    dfloat zero = 0.;
    solver->addBCKernel(mesh->Nelements,
                       zero,
                       mesh->o_x,
                       mesh->o_y,
                       mesh->o_z,
                       mesh->o_mapB,
                       o_x);
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
      dfloat exact = sin(M_PI*xn)*sin(M_PI*yn)*sin(M_PI*zn);
      dfloat error = fabs(exact-mesh->q[id]);

      maxError = mymax(maxError, error);

      //mesh->q[id] -= r[id];
    }
  }

  dfloat globalMaxError = 0;
  MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  if(rank==0)
    printf("globalMaxError = %g\n", globalMaxError);

  char filename[BUFSIZ];
  sprintf(filename, "foo_%d.vtu", rank);
  meshPlotVTU3D(mesh, filename, 0);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
