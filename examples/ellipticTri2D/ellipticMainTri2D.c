#include "ellipticTri2D.h"

void applyElementMatrix(mesh_t *mesh, dfloat *A, dfloat *q, dfloat *Aq) {

  dfloat *Aqn = (dfloat*) calloc(mesh->Np,sizeof(dfloat));
  for (iint e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) {
      Aqn[n] = 0;
      for (int k=0;k<mesh->Np;k++) {
        Aqn[n] += A[k+n*mesh->Np]*q[k+e*mesh->Np];
      }
    }
    for (int n=0;n<mesh->Np;n++) Aq[n+e*mesh->Np] = Aqn[n];
  }
  free(Aqn);
}

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

  // solver can be PCG, PGMRES, or PBiCGStab
  // can add FLEXIBLE and VERBOSE options
  // method can be IPDG or CONTINUOUS
  //  can add NONSYM option
  // basis can be NODAL or BERN
  // preconditioner can be NONE, JACOBI, OAS, MASSMATRIX, FULLALMOND, or MULTIGRID
  // OAS and MULTIGRID: smoothers can be FULLPATCH, FACEPATCH, LOCALPATCH, OVERLAPPINGPATCH, or DAMPEDJACOBI
  //                      patch smoothers can include EXACT
  // MULTIGRID: smoothers can include CHEBYSHEV for smoother acceleration
  // MULTIGRID: levels can be ALLDEGREES, HALFDEGREES, HALFDOFS
  // FULLALMOND: can include MATRIXFREE option
  char *options =
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG basis=NODAL preconditioner=OAS smoother=FULLPATCH");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=BRDG basis=BERN preconditioner=MULTIGRID,HALFDOFS smoother=CHEBYSHEV");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=CONTINUOUS basis=SPARSE preconditioner=MULTIGRID,HALFDOFS smoother=DAMPEDJACOBI,CHEBYSHEV");
    strdup("solver=PCG,FLEXIBLE,VERBOSE method=CONTINUOUS basis=NODAL preconditioner=SEMFEM");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=IPDG basis=NODAL preconditioner=JACOBI");

  //FULLALMOND, OAS, and MULTIGRID will use the parAlmondOptions in setup
  // solver can be KCYCLE, or VCYCLE
  //  can add the EXACT and NONSYM option
  // smoother can be DAMPEDJACOBI or CHEBYSHEV
  // partition can be STRONGNODES, DISTRIBUTED, SATURATE
  char *parAlmondOptions =
    strdup("solver=KCYCLE,VERBOSE smoother=CHEBYSHEV partition=STRONGNODES");
  //strdup("solver=EXACT,VERBOSE smoother=CHEBYSHEV partition=STRONGNODES");

  //this is strictly for testing, to do repeated runs. Will be removed later
  //  if (argc==6) {
  //   options = strdup(argv[4]);
  //   parAlmondOptions = strdup(argv[5]);
  // }

  // set up mesh stuff
  mesh2D *mesh = meshSetupTri2D(argv[1], N);

  // parameter for elliptic problem (-laplacian + lambda)*q = f
  //dfloat lambda = 1;
  dfloat lambda = 0;

  if (strstr(options,"SPARSE")&&(lambda!=0)) { //sanity check
    printf("SPARSE not currently supported for screened Poisson\n");
    exit(-1);
  }

  // set up
  occa::kernelInfo kernelInfo;
  ellipticSetupTri2D(mesh, kernelInfo, options);

  // Boundary Type translation. Just default from the mesh file.
  int BCType[3] = {0,1,2};

  dfloat tau;
  if (strstr(options,"IPDG")) {
    tau = 2.0*(mesh->N+1)*(mesh->N+2)/2.0;
  } else if (strstr(options,"BRDG")) {
    tau = 1.0;
  }

  solver_t *solver = ellipticSolveSetupTri2D(mesh, tau, lambda, BCType, kernelInfo, options, parAlmondOptions);
  
  iint Nall = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  dfloat *r   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  dfloat *x   = (dfloat*) calloc(Nall,   sizeof(dfloat));

  // load forcing into r
  for(iint e=0;e<mesh->Nelements;++e){
    dfloat J = mesh->vgeo[e*mesh->Nvgeo+JID];
    for(iint n=0;n<mesh->Np;++n){
      iint id = n+e*mesh->Np;
      dfloat xn = mesh->x[id];
      dfloat yn = mesh->y[id];
      r[id] = J*(2*M_PI*M_PI+lambda)*sin(M_PI*xn)*sin(M_PI*yn);
      x[id] = 0;
    }
  }

  //Apply some element matrix ops to r depending on our solver
  if (strstr(options,"BERN"))   applyElementMatrix(mesh,mesh->invVB,r,r);
  if (strstr(options,"SPARSE")) applyElementMatrix(mesh,mesh->invSparseV,r,r);

  if (!strstr(options,"NONSYM")) {
    if (strstr(options,"NODAL"))  applyElementMatrix(mesh,mesh->MM,r,r);
    if (strstr(options,"BERN"))   applyElementMatrix(mesh,mesh->BBMM,r,r);
    if (strstr(options,"SPARSE")) applyElementMatrix(mesh,mesh->sparseMM,r,r);
  }

  //copy to occa buffers
  occa::memory o_r   = mesh->device.malloc(Nall*sizeof(dfloat), r);
  occa::memory o_x   = mesh->device.malloc(Nall*sizeof(dfloat), x);

  // capture header file
  char *boundaryHeaderFileName;
  // if(argc==3)
  boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticTri2D/homogeneous2D.h"); // default
  //  else
  //  boundaryHeaderFileName = strdup(argv[3]);
  //add user defined boundary data
  kernelInfo.addInclude(boundaryHeaderFileName);

  //add boundary condition contribution to rhs
  if (strstr(options,"IPDG")) {

    solver->rhsBCIpdgKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticRhsBCIpdgTri2D.okl",
          "ellipticRhsBCIpdgTri2D",
          kernelInfo);

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
                            solver->o_EToB,
                            mesh->o_DrT,
                            mesh->o_DsT,
                            mesh->o_LIFTT,
                            mesh->o_MM,
                            o_r);
  }

  if (strstr(options,"CONTINUOUS")) {

    solver->rhsBCKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticRhsBCTri2D.okl",
          "ellipticRhsBCTri2D",
          kernelInfo);

    solver->addBCKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticAddBCTri2D.okl",
          "ellipticAddBCTri2D",
          kernelInfo);

    dfloat zero = 0.f;
    solver->rhsBCKernel(mesh->Nelements,
                        mesh->o_ggeo,
                        mesh->o_sgeo,
                        mesh->o_SrrT,
                        mesh->o_SrsT,
                        mesh->o_SsrT,
                        mesh->o_SssT,
                        mesh->o_MM,
                        mesh->o_vmapM,
                        mesh->o_sMT,
                        lambda,
                        zero,
                        mesh->o_x,
                        mesh->o_y,
                        mesh->o_mapB,
                        o_r);
  }

  // gather-scatter
  if(strstr(options, "CONTINUOUS")){
    //sign correction for gs
    if (strstr(options,"SPARSE")) solver->dotMultiplyKernel(mesh->Np*mesh->Nelements, o_r, mesh->o_mapSgn, o_r);
    ellipticParallelGatherScatterTri2D(mesh, mesh->ogs, o_r, o_r, dfloatString, "add");  
    if (strstr(options,"SPARSE")) solver->dotMultiplyKernel(mesh->Np*mesh->Nelements, o_r, mesh->o_mapSgn, o_r);       
    //mask
    if (mesh->Nmasked) mesh->maskKernel(mesh->Nmasked, mesh->o_maskIds, o_r);
  }


  // convergence tolerance
  dfloat tol = 1e-8;
  ellipticSolveTri2D(solver, lambda, tol, o_r, o_x, options);



  if(strstr(options, "CONTINUOUS")){
    dfloat zero = 0.;
    solver->addBCKernel(mesh->Nelements,
                       zero,
                       mesh->o_x,
                       mesh->o_y,
                       mesh->o_mapB,
                       o_x);
  }

  // copy solution from DEVICE to HOST
  o_x.copyTo(mesh->q);

  if (strstr(options,"BERN"))   applyElementMatrix(mesh,mesh->VB,mesh->q,mesh->q);
  if (strstr(options,"SPARSE")) applyElementMatrix(mesh,mesh->sparseV,mesh->q,mesh->q);

  dfloat maxError = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      iint   id = e*mesh->Np+n;
      dfloat xn = mesh->x[id];
      dfloat yn = mesh->y[id];
      dfloat exact = sin(M_PI*xn)*sin(M_PI*yn);
      dfloat error = fabs(exact-mesh->q[id]);

      maxError = mymax(maxError, error);
      //mesh->q[id] -= exact;
    }
  }

  dfloat globalMaxError = 0;
  MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  if(rank==0)
    printf("globalMaxError = %g\n", globalMaxError);

  char filename[BUFSIZ];
  sprintf(filename, "foo_%d.vtu", rank);
  meshPlotVTU2D(mesh, filename, 0);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
