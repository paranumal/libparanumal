#include "ellipticTri2D.h"

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
    strdup("solver=PCG,FLEXIBLE,VERBOSE method=CONTINUOUS basis=NODAL preconditioner=FULLALMOND");
    //strdup("solver=PCG,FLEXIBLE,VERBOSE method=CONTINUOUS basis=NODAL preconditioner=NONE");
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

  iint NblockV = 1;
  iint NnodesV = 1;

  if (argc == 5){
    NblockV = atoi(argv[3]);
    NnodesV = atoi(argv[4]);
  }

  // set up mesh stuff
  mesh2D *mesh = meshSetupTri2D(argv[1], N);
  ogs_t *ogs;
  precon_t *precon;

  if (strstr(options,"SPARSE")) {
    //connect the face modes between each element 
    meshConnectFaceModes2D(mesh,mesh->FaceModes,mesh->sparseV);
    //use the mmaps constructed and overwrite vmap and FaceNodes
    for (iint n=0;n<mesh->Nfp*mesh->Nfaces*mesh->Nelements;n++) {
      mesh->vmapM[n] = mesh->mmapM[n];
      mesh->vmapP[n] = mesh->mmapP[n];
    }
    for (int n=0;n<mesh->Nfaces*mesh->Nfp;n++) {
      mesh->faceNodes[n] = mesh->FaceModes[n];
    }
    //free the old gather scatter arrays and re-run the connect nodes function using this updated connectivity
    free(mesh->gatherLocalIds );
    free(mesh->gatherBaseIds  );
    free(mesh->gatherBaseRanks);
    free(mesh->gatherHaloFlags);
    meshParallelConnectNodes(mesh);
  }

  // parameter for elliptic problem (-laplacian + lambda)*q = f
  //dfloat lambda = 1;
  dfloat lambda = 1;

  // set up
  occa::kernelInfo kernelInfo;
  ellipticSetupTri2D(mesh, kernelInfo);

  // Boundary Type translation. Just default from the mesh file.
  int BCType[3] = {0,2,2};

  dfloat tau;
  if (strstr(options,"IPDG")) {
    tau = 2.0*(mesh->N+1)*(mesh->N+2)/2.0;
  } else if (strstr(options,"BRDG")) {
    tau = 1.0;
  }
  if (strstr(options, "SPARSE")){
    //loadElementStiffnessMatricesTri2D(mesh, options, mesh->N);
  } else{
    //loadElementStiffnessMatricesTri2D(mesh, options, mesh->N);
    //buildElementStiffnessMatricesTri2D(mesh, options, mesh->N);
  } 

  solver_t *solver = ellipticSolveSetupTri2D(mesh, tau, lambda, BCType, kernelInfo, options, parAlmondOptions, NblockV, NnodesV);
  iint Nall = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  dfloat *r   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  dfloat *x   = (dfloat*) calloc(Nall,   sizeof(dfloat));

  // load rhs into r
  dfloat *nrhs = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *nrhstmp = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(iint e=0;e<mesh->Nelements;++e){
    dfloat J = mesh->vgeo[e*mesh->Nvgeo+JID];
    for(iint n=0;n<mesh->Np;++n){
      dfloat xn = mesh->x[n+e*mesh->Np];
      dfloat yn = mesh->y[n+e*mesh->Np];
      nrhs[n] = -(2*M_PI*M_PI+lambda)*sin(M_PI*xn)*sin(M_PI*yn);
    }

    if (strstr(options,"BERN")) {
      for(iint n=0;n<mesh->Np;++n){
        nrhstmp[n] = 0.;
        for(iint m=0;m<mesh->Np;++m){
          nrhstmp[n] += mesh->invVB[n*mesh->Np+m]*nrhs[m];
        }
      }
      for(iint n=0;n<mesh->Np;++n){
        dfloat rhs = 0;
        if (strstr(options,"NONSYM")) {
          rhs = nrhs[n];
        } else {
          for(iint m=0;m<mesh->Np;++m){
            rhs += mesh->BBMM[n+m*mesh->Np]*nrhstmp[m];
          }
        }
        iint id = n+e*mesh->Np;

        r[id] = -rhs*J;
        x[id] = 0;
        mesh->q[id] = rhs;
      }
    } else if (strstr(options,"NODAL")) {
      for(iint n=0;n<mesh->Np;++n){
        dfloat rhs = 0;
        if (strstr(options,"NONSYM")) {
          rhs = nrhs[n];
        } else {
          for(iint m=0;m<mesh->Np;++m){
            rhs += mesh->MM[n+m*mesh->Np]*nrhs[m];
          }
        }
        iint id = n+e*mesh->Np;

        r[id] = -rhs*J;
        x[id] = 0;
        mesh->q[id] = rhs;
      }
    }
  }
  free(nrhs);
  free(nrhstmp);

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

  // convergence tolerance
  dfloat tol = 1e-6;
  ellipticSolveTri2D(solver, lambda, tol, o_r, o_x, options, NblockV, NnodesV);

  // copy solution from DEVICE to HOST
  o_x.copyTo(mesh->q);

  if (strstr(options,"BERN")) {
    dfloat *qtmp = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
    for (iint e =0;e<mesh->Nelements;e++){
      iint id = e*mesh->Np;

      for (iint n=0; n<mesh->Np; n++){
        qtmp[n] = mesh->q[id+n];
        mesh->q[id+n] = 0.0;
      }
      for (iint n=0;n<mesh->Np;n++){
        for (iint m=0; m<mesh->Np; m++){
          mesh->q[id+n] += mesh->VB[n*mesh->Np+m]*qtmp[m];
        }
      }
    }
    free(qtmp);
  }

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

  meshPlotVTU2D(mesh, "foo.vtu", 0);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
