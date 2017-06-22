#include "ins2D.h"

// NBN: toggle use of 2nd stream
#define USE_2_STREAMS

ins_t *insSetup2D(mesh2D *mesh, char * options, char *vSolverOptions, char *pSolverOptions, char *boundaryHeaderFileName){

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", (rank)%3);
  //  sprintf(deviceConfig, "mode = CUDA, deviceID = 0");
  //printf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 0");
  //sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");  

  ins_t *ins = (ins_t*) calloc(1, sizeof(ins_t));

  ins->NVfields = 2; //  Total Number of Velocity Fields
  ins->NTfields = 3; // Total Velocity + Pressure
  ins->Nfields  = 1; // Each Velocity Field
  ins->ExplicitOrder = 3; // Order Nonlinear Extrapolation

  mesh->Nfields = ins->Nfields;

  ins->mesh = mesh;
  int Nstages = 4;
  // compute samples of q at interpolation nodes
  ins->U     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->V     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->P     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));

  //rhs storage
  ins->rhsU  = (dfloat*) calloc(2*mesh->Nelements*mesh->Np,sizeof(dfloat));
  ins->rhsV  = (dfloat*) calloc(2*mesh->Nelements*mesh->Np,sizeof(dfloat));
  ins->rhsP  = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));

  //additional field storage (could reduce in the future)
  ins->NU     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->NV     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->Px     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->Py     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->PI     = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));


  if(strstr(options,"SUBCYCLING")){

    ins->Nsubsteps = 1; //was 3

    ins->Ud   = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
    ins->Vd   = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
    ins->Ue   = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
    ins->Ve   = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
    ins->resU = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
    ins->resV = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));

  }


  // SET SOLVER OPTIONS
  // Initial Conditions, Flow Properties
  printf("Starting initial conditions for INS2D\n");
  //
  dfloat ux   = 0.0  ;
  dfloat uy   = 0.0  ;
  dfloat pr   = 0.0  ;
  dfloat nu   = 0.025;   // kinematic viscosity,
  dfloat rho  = 1.0  ;  // Give density for getting actual pressure in nondimensional solve

  dfloat g[2]; g[0] = 0.0; g[1] = 0.0;  // No gravitational acceleration

  // Fill up required fileds
  ins->finalTime = 1.0;
  ins->nu        = nu ;
  ins->rho       = rho;
  ins->tau       = 4.*(mesh->N+1)*(mesh->N+1);

  // Define total DOF per field for INS i.e. (Nelelemts + Nelements_halo)*Np
  ins->NtotalDofs = (mesh->totalHaloPairs+mesh->Nelements)*mesh->Np ;
  ins->NDofs      = mesh->Nelements*mesh->Np;
  // Initialize
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      const iint id = n + mesh->Np*e;
      dfloat t = 0;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];
#if 1
      dfloat lambda = 1./(2.*ins->nu)-sqrt(1./(4.*ins->nu*ins->nu) + 4.*M_PI*M_PI) ;
      //
      ins->U[id] = 1.0 - exp(lambda*x)*cos(2.*M_PI*y);
      ins->V[id] = lambda/(2.*M_PI)*exp(lambda*x)*sin(2.*M_PI*y);
      ins->P[id] = 0.5*(1.0- exp(2.*lambda*x));
#endif

#if 0
      ins->U[id] = y*(4.5f-y)/(2.25f*2.25f);
      ins->V[id] = 0;
      ins->P[id] = (nu*(-2.)/(2.25*2.25))*(x-4.5) ;
#endif

#if 0
      ins->U[id] = -sin(2.0 *M_PI*y)*exp(-ins->nu*4.0*M_PI*M_PI*0.0); ;
      ins->V[id] =  sin(2.0 *M_PI*x)*exp(-ins->nu*4.0*M_PI*M_PI*0.0); 
      ins->P[id] = -cos(2.0 *M_PI*y)*cos(2.f*M_PI*x)*exp(-nu*8.f*M_PI*M_PI*0.0);
#endif


#if 0
      ins->U[id] = 1;
      ins->V[id] = 0;
      ins->P[id] = 0;
#endif
    }
  }

  printf("starting parameters\n");
  // set time step
  dfloat hmin = 1e9, hmax = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint f=0;f<mesh->Nfaces;++f){
      iint sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      dfloat hest = 2./(sJ*invJ);

      hmin = mymin(hmin, hest);
      hmax = mymax(hmax, hest);
    }
  }

  // Find Maximum Velocity
  dfloat umax = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      const iint id = n + mesh->Np*e;
      dfloat t = 0;
      dfloat uxn = ins->U[id];
      dfloat uyn = ins->V[id];

      //Squared maximum velocity
      dfloat numax = uxn*uxn + uyn*uyn;
      umax = mymax(umax, numax);
    }
  }
  // Maximum Velocity
  umax = sqrt(umax);

  dfloat cfl = 0.25; // pretty good estimate (at least for subcycling LSERK4)
  dfloat magVel = mymax(umax,1.0); // Correction for initial zero velocity
  dfloat dt = cfl* hmin/( (mesh->N+1.)*(mesh->N+1.) * magVel) ;

  printf("hmin = %g\n", hmin);
  printf("hmax = %g\n", hmax);
  printf("cfl = %g\n", cfl);
  printf("dt = %g\n", dt);

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(ins->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  if(strstr(options,"SUBCYCLING")){
    ins->NtimeSteps  = ins->finalTime/ins->dt;
    ins->sdt         = ins->finalTime/ins->NtimeSteps;

    ins->dt          = ins->sdt*ins->Nsubsteps;
    ins->NtimeSteps  = ins->NtimeSteps/ins->Nsubsteps;

    printf("dt: %.8f and sdt: %.8f ratio: %.8f \n", ins->dt, ins->sdt, ins->dt/ins->sdt);

  }
  else{
    ins->NtimeSteps = ins->finalTime/ins->dt;
    ins->dt   = ins->finalTime/ins->NtimeSteps;

  }


  // errorStep
  ins->errorStep =1;

  printf("Nsteps = %d NerrStep= %d dt = %.8e\n", ins->NtimeSteps,ins->errorStep, ins->dt);

  occa::kernelInfo kernelInfo;
  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

  //add boundary data to kernel info
  kernelInfo.addInclude(boundaryHeaderFileName);

  occa::kernelInfo kernelInfoV  = kernelInfo;
  occa::kernelInfo kernelInfoP  = kernelInfo;

  printf("==================ELLIPTIC SOLVE SETUP=========================\n");

  // SetUp Boundary Flags types for Elliptic Solve
  int vBCType[4] = {0,1,1,2}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int pBCType[4] = {0,2,2,1}; // bc=3 => outflow => Dirichlet => pBCType[3] = 1, etc.

  // Use third Order Velocity Solve: full rank should converge for low orders
  printf("==================VELOCITY SOLVE SETUP=========================\n");
  ins->lambda = (11.f/6.f) / (ins->dt * ins->nu);
  boundaryHeaderFileName = strdup(DHOLMES "/examples/insTri2D/insVelocityEllipticBC2D.h");
  kernelInfoV.addInclude(boundaryHeaderFileName);
  solver_t *vSolver   = ellipticSolveSetupTri2D(mesh, ins->tau, ins->lambda, vBCType, kernelInfoV, vSolverOptions);
  ins->vSolver        = vSolver;
  ins->vSolverOptions = vSolverOptions;

  printf("==================PRESSURE SOLVE SETUP========================\n");
  // SETUP PRESSURE and VELOCITY SOLVERS
  boundaryHeaderFileName = strdup(DHOLMES "/examples/insTri2D/insPressureEllipticBC2D.h");
  kernelInfoP.addInclude(boundaryHeaderFileName);
  solver_t *pSolver   = ellipticSolveSetupTri2D(mesh, ins->tau, 0.0, pBCType,kernelInfoP, pSolverOptions);
  ins->pSolver        = pSolver;
  ins->pSolverOptions = pSolverOptions;

 

  kernelInfo.addDefine("p_maxNodesVolume", mymax(mesh->cubNp,mesh->Np));
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int maxSurfaceNodes = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  kernelInfo.addDefine("p_maxSurfaceNodes", maxSurfaceNodes);
  printf("maxSurfaceNodes=%d\n", maxSurfaceNodes);


  int NblockV = 256/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 256/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  printf("maxNodes: %d \t NblockV: %d \t NblockS: %d  \n", maxNodes, NblockV, NblockS);

  printf("Np: %d \t Ncub: %d \n", mesh->Np, mesh->cubNp);

  // ADD-DEFINES
  kernelInfo.addDefine("p_Lambda2", 0.5f);
  kernelInfo.addDefine("p_NTfields", ins->NTfields);
  kernelInfo.addDefine("p_NVfields", ins->NVfields);
  kernelInfo.addDefine("p_NfacesNfp",  mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_inu",      (float) 1.f/ins->nu);
  kernelInfo.addDefine("p_nu",      (float) ins->nu);
  kernelInfo.addDefine("p_idt",      (float) 1.f/ins->dt);

  printf("mesh nfields %d\n", mesh->Nfields);
  // MEMORY ALLOCATION
  ins->o_U = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->U);
  ins->o_V = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->V);
  ins->o_P = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->P);

  ins->o_rhsU  = mesh->device.malloc(2*mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsU);
  ins->o_rhsV  = mesh->device.malloc(2*mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsV);
  ins->o_rhsP  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsP);

  ins->o_NU = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->NU);
  ins->o_NV = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->NV);
  ins->o_Px = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->Px);
  ins->o_Py = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->Py);
  ins->o_PI = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->PI);
  ins->o_PIx = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat));
  ins->o_PIy = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat));

  //storage for helmholtz solves. Fix this later
  ins->o_UH = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat));
  ins->o_VH = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat));


  if(strstr(options,"SUBCYCLING")){
    // Note that resU and resV can be replaced with already introduced buffer
    ins->o_Ue   = mesh->device.malloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*sizeof(dfloat), ins->Ue);
    ins->o_Ve   = mesh->device.malloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*sizeof(dfloat), ins->Ve);
    ins->o_resU = mesh->device.malloc((mesh->Nelements)*mesh->Np*sizeof(dfloat), ins->resU);
    ins->o_resV = mesh->device.malloc((mesh->Nelements)*mesh->Np*sizeof(dfloat), ins->resV);


    printf("Compiling SubCycle Advection volume kernel \n");
    ins->subCycleVolumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",
					 "insSubCycleVolume2D",
					 kernelInfo);

    printf("Compiling SubCycle Advection surface kernel\n");
    ins->subCycleSurfaceKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",
					 "insSubCycleSurface2D",
					 kernelInfo);

    printf("Compiling SubCycle Advection cubature  volume kernel \n");
    ins->subCycleCubatureVolumeKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",
					 "insSubCycleCubatureVolume2D",
					 kernelInfo);

    printf("Compiling SubCycle Advection cubature surface kernel\n");
    ins->subCycleCubatureSurfaceKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",
					 "insSubCycleCubatureSurface2D",
					 kernelInfo);


    printf("Compiling SubCycle Advection RK update kernel\n");
    ins->subCycleRKUpdateKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",
					 "insSubCycleRKUpdate2D",
					 kernelInfo);



    ins->subCycleExtKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",
					 "insSubCycleExt2D",
					 kernelInfo);

  }


  if(mesh->totalHaloPairs){
    ins->o_tHaloBuffer = mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*(ins->NTfields)*sizeof(dfloat));
    ins->o_vHaloBuffer = mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*(ins->NVfields)*sizeof(dfloat));
    ins->o_pHaloBuffer = mesh->device.malloc(mesh->totalHaloPairs*mesh->Np *sizeof(dfloat));
  }

  ins->mesh = mesh;

  // ===========================================================================
  printf("Compiling Advection volume kernel with cubature integration\n");
  ins->advectionCubatureVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection2D.okl",
				       "insAdvectionCubatureVolume2D",
				       kernelInfo);

  printf("Compiling Advection surface kernel with cubature integration\n");
  ins->advectionCubatureSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection2D.okl",
				       "insAdvectionCubatureSurface2D",
				       kernelInfo);

  printf("Compiling Advection volume kernel with collocation integration\n");
  ins->advectionVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection2D.okl",
				       "insAdvectionVolume2D",
				       kernelInfo);

  printf("Compiling Advection surface kernel with collocation integration\n");
  ins->advectionSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection2D.okl",
				       "insAdvectionSurface2D",
				       kernelInfo);

  // ===========================================================================
  printf("Compiling Gradient volume kernel with collocation integration\n");
  ins->gradientVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insGradient2D.okl",
				       "insGradientVolume2D",
				       kernelInfo);

  printf("Compiling Gradient volume kernel with collocation integration\n");
  ins->gradientSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insGradient2D.okl",
				       "insGradientSurface2D",
				       kernelInfo);

  // ===========================================================================
  printf("Compiling Divergence volume kernel with collocation integration\n");
  ins->divergenceVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insDivergence2D.okl",
				       "insDivergenceVolume2D",
				       kernelInfo);

  printf("Compiling Divergencesurface kernel with collocation integration\n");
  ins->divergenceSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insDivergence2D.okl",
				       "insDivergenceSurface2D",
				       kernelInfo);

  // ===========================================================================
  printf("Compiling Helmholtz volume update kernel\n");
  ins->helmholtzRhsForcingKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHelmholtzRhs2D.okl",
				       "insHelmholtzRhsForcing2D",
				       kernelInfo);

  printf("Compiling Helmholtz IPDG RHS kernel with collocation integration\n");
  ins->helmholtzRhsIpdgBCKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHelmholtzRhs2D.okl",
				       "insHelmholtzRhsIpdgBC2D",
				       kernelInfo);


  printf("Compiling INS Helmholtz Halo Extract Kernel\n");
  ins->totalHaloExtractKernel=
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
				       "insTotalHaloExtract2D",
				       kernelInfo);

  printf("Compiling INS Helmholtz Halo Extract Kernel\n");
  ins->totalHaloScatterKernel=
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
				       "insTotalHaloScatter2D",
				       kernelInfo);


  // ===========================================================================
  printf("Compiling Helmholtz volume update kernel\n");
  ins->poissonRhsForcingKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonRhs2D.okl",
				       "insPoissonRhsForcing2D",
				       kernelInfo);

  printf("Compiling Poisson IPDG surface kernel with collocation integration\n");
  ins->poissonRhsIpdgBCKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonRhs2D.okl",
				       "insPoissonRhsIpdgBC2D",
				       kernelInfo);

  printf("Compiling Poisson penalty surface kernel\n");
  ins->poissonPenaltyKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonPenalty2D.okl",
				       "insPoissonPenalty2D",
				       kernelInfo);

  printf("Compiling INS Poisson Halo Extract Kernel\n");
  ins->velocityHaloExtractKernel=
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
				       "insVelocityHaloExtract2D",
				       kernelInfo);

  printf("Compiling INS PoissonHalo Extract Kernel\n");
  ins->velocityHaloScatterKernel=
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
				       "insVelocityHaloScatter2D",
				       kernelInfo);

  printf("Compiling Poisson surface kernel with collocation integration\n");
  ins->updateUpdateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insUpdate2D.okl",
				       "insUpdateUpdate2D",
				       kernelInfo);


  printf("Compiling INS Poisson Halo Extract Kernel\n");
  ins->pressureHaloExtractKernel=
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
				       "insPressureHaloExtract2D",
				       kernelInfo);

  printf("Compiling INS PoissonHalo Extract Kernel\n");
  ins->pressureHaloScatterKernel=
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
				       "insPressureHaloScatter2D",
				       kernelInfo);
  // ===========================================================================//

  void insBuildVectorIpdgTri2D(mesh2D *mesh, dfloat tau, dfloat sigma, dfloat lambda,
			       iint *BCType, nonZero_t **A, iint *nnzA,
			       hgs_t **hgs, iint *globalStarts, const char *options);

  nonZero_t *A;
  iint nnz;
  hgs_t *hgs;
  iint *globalStarts = (iint*) calloc(size+1, sizeof(iint));
  dfloat sigma = 100;

  MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT, globalStarts+1, 1, MPI_IINT, MPI_COMM_WORLD);
  for(iint r=0;r<size;++r)
    globalStarts[r+1] = globalStarts[r]+globalStarts[r+1]*mesh->Np*2;
  
  insBuildVectorIpdgTri2D(mesh, ins->tau, sigma, ins->lambda, vBCType, &A, &nnz,&hgs,globalStarts, vSolverOptions);

  //collect global assembled matrix
  iint *globalnnz       = (iint *) calloc(size  ,sizeof(iint));
  iint *globalnnzOffset = (iint *) calloc(size+1,sizeof(iint));
  MPI_Allgather(&nnz, 1, MPI_IINT, 
                globalnnz, 1, MPI_IINT, MPI_COMM_WORLD);
  globalnnzOffset[0] = 0;
  for (iint n=0;n<size;n++)
    globalnnzOffset[n+1] = globalnnzOffset[n]+globalnnz[n];

  iint globalnnzTotal = globalnnzOffset[size];

  iint *globalRecvCounts  = (iint *) calloc(size,sizeof(iint));
  iint *globalRecvOffsets = (iint *) calloc(size,sizeof(iint));
  for (iint n=0;n<size;n++){
    globalRecvCounts[n] = globalnnz[n]*sizeof(nonZero_t);
    globalRecvOffsets[n] = globalnnzOffset[n]*sizeof(nonZero_t);
  }
  nonZero_t *globalNonZero = (nonZero_t*) calloc(globalnnzTotal, sizeof(nonZero_t));

  MPI_Allgatherv(A, nnz*sizeof(nonZero_t), MPI_CHAR, 
                globalNonZero, globalRecvCounts, globalRecvOffsets, MPI_CHAR, MPI_COMM_WORLD);
  
  iint *globalRows = (iint *) calloc(globalnnzTotal, sizeof(iint));
  iint *globalCols = (iint *) calloc(globalnnzTotal, sizeof(iint));
  dfloat *globalVals = (dfloat*) calloc(globalnnzTotal,sizeof(dfloat));

  for (iint n=0;n<globalnnzTotal;n++) {
    globalRows[n] = globalNonZero[n].row;
    globalCols[n] = globalNonZero[n].col;
    globalVals[n] = globalNonZero[n].val;
  }

  ins->precon = (precon_t *) calloc(1,sizeof(precon_t));
  ins->precon->parAlmond = parAlmondSetup(mesh, 
					  globalStarts, 
					  globalnnzTotal,      
					  globalRows,        
					  globalCols,        
					  globalVals,
					  hgs,
					  vSolverOptions);       //rhs will be passed gather-scattered
  
  return ins;
}







