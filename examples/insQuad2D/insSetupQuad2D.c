#include "insQuad2D.h"

ins_t *insSetupQuad2D(mesh2D *mesh, int Ns, char * options, 
                  char *vSolverOptions, char *vParAlmondOptions,
                  char *pSolverOptions, char *pParAlmondOptions,
                  char *boundaryHeaderFileName){

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  long int hostId = gethostid();

  long int* hostIds = (long int*) calloc(size,sizeof(long int));
  MPI_Allgather(&hostId,1,MPI_LONG,hostIds,1,MPI_LONG,MPI_COMM_WORLD);

  int deviceID = 0;
  for (int r=0;r<rank;r++) {
    if (hostIds[r]==hostId) deviceID++;
  }

  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", deviceID);
  //sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 0");
  //sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");  

  ins_t *ins = (ins_t*) calloc(1, sizeof(ins_t));

  ins->NVfields = 2; //  Total Number of Velocity Fields
  ins->NTfields = 3; // Total Velocity + Pressure
  ins->Nfields  = 1; // Each Velocity Field
  //ins->ExplicitOrder = 3; // Order Nonlinear Extrapolation

  mesh->Nfields = ins->Nfields;

  ins->mesh = mesh;

  dlong Ntotal = mesh->Np*mesh->Nelements;
  ins->Nblock = (Ntotal+blockSize-1)/blockSize;

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

  ins->Vx     = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  ins->Vy     = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  ins->Div     = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));

  
  ins->Nsubsteps = Ns;
  if(strstr(options,"SUBCYCLING")){
    ins->Ud   = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
    ins->Vd   = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
    ins->Ue   = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
    ins->Ve   = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
    ins->resU = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
    ins->resV = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  }

  // SET SOLVER OPTIONS
  dfloat ux   = 0.0  ;
  dfloat uy   = 0.0  ;
  dfloat pr   = 0.0  ;
  dfloat nu   = 0.01;   // kinematic viscosity,
  dfloat rho  = 1.0  ;  // Give density for getting actual pressure in nondimensional solve

  dfloat g[2]; g[0] = 0.0; g[1] = 0.0;  // No gravitational acceleration

  // Fill up required fileds
  ins->finalTime = 50.0;
  ins->nu        = nu ;
  ins->rho       = rho;
  ins->tau       = 2.0* (mesh->N+1)*(mesh->N+2)/2.0f; 

  // Define total DOF per field for INS i.e. (Nelelemts + Nelements_halo)*Np
  ins->NtotalDofs = (mesh->totalHaloPairs+mesh->Nelements)*mesh->Np ;
  ins->NDofs      = mesh->Nelements*mesh->Np;
  // Initialize
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      const dlong id = n + mesh->Np*e;
      dfloat t = 0;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];
      #if 0
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
            ins->P[id] = -cos(2.0 *M_PI*y)*cos(2.f*M_PI*x)*exp(-ins->nu*8.f*M_PI*M_PI*0.0);
      #endif


      #if 1 // Zero flow
            ins->U[id] = 1.0;
            ins->V[id] = 0.0;
            ins->P[id] = 0.0;
      #endif
    }
  }

  // set time step
  dfloat hmin = 1e9, hmax = 0;
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      dfloat hest = 2./(sJ*invJ);

      hmin = mymin(hmin, hest);
      hmax = mymax(hmax, hest);
    }
  }

  // Find Maximum Velocity
  dfloat umax = 0;
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      const dlong id = n + mesh->Np*e;
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

 
  dfloat cfl = 1.0; // pretty good estimate (at least for subcycling LSERK4)
  dfloat magVel = mymax(umax,1.0); // Correction for initial zero velocity
  dfloat dt     = cfl* hmin/( (mesh->N+1.)*(mesh->N+1.) * magVel) ;

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(ins->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  if(strstr(options,"SUBCYCLING")){
    ins->dt         = ins->Nsubsteps*ins->dt;
    ins->NtimeSteps = ins->finalTime/ins->dt;
    ins->dt         = ins->finalTime/ins->NtimeSteps;
    ins->sdt        = ins->dt/ins->Nsubsteps;
  } else{
    ins->NtimeSteps = ins->finalTime/ins->dt;
    ins->dt         = ins->finalTime/ins->NtimeSteps;
  }

  if (rank==0) {
    printf("hmin = %g\n", hmin);
    printf("hmax = %g\n", hmax);
    printf("cfl = %g\n",  cfl);
    printf("dt = %g\n",   ins->dt);
  }
  
  if (strstr(options,"SUBCYCLING")&&rank==0) printf("dt: %.8f and sdt: %.8f ratio: %.8f \n", ins->dt, ins->sdt, ins->dt/ins->sdt);
  
  // Hold some inverses for kernels
  ins->inu = 1.0/ins->nu; 
  ins->idt = 1.0/ins->dt;
  
  // errorStep
  if(strstr(options,"SUBCYCLING"))
    // ins->errorStep =100*32/ins->Nsubsteps;
    //ins->errorStep =800/ins->Nsubsteps;
    ins->errorStep = 400;
  else
    ins->errorStep = 400;

  if (rank==0) printf("Nsteps = %d NerrStep= %d dt = %.8e\n", ins->NtimeSteps,ins->errorStep, ins->dt);

  occa::kernelInfo kernelInfo;
  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

  //add boundary data to kernel info
  kernelInfo.addInclude(boundaryHeaderFileName);

  occa::kernelInfo kernelInfoV  = kernelInfo;
  occa::kernelInfo kernelInfoP  = kernelInfo;

  //if (rank!=0) 
    occa::setVerboseCompilation(false);

  if (rank==0) printf("==================ELLIPTIC SOLVE SETUP=========================\n");

  // SetUp Boundary Flags types for Elliptic Solve
  int vBCType[4] = {0,1,1,2}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int pBCType[4] = {0,2,2,1}; // bc=3 => outflow => Dirichlet => pBCType[3] = 1, etc.

  //Solver tolerances 
  ins->presTOL = 1E-8;
  ins->velTOL  = 1E-8;

  // Use third Order Velocity Solve: full rank should converge for low orders
  if (rank==0) printf("==================VELOCITY SOLVE SETUP=========================\n");
  // ins->lambda = (11.f/6.f) / (ins->dt * ins->nu);
  ins->lambda = (1.5f) / (ins->dt * ins->nu);
  solver_t *vSolver   = ellipticSolveSetupQuad2D(mesh, ins->tau, ins->lambda, vBCType, kernelInfoV, vSolverOptions,vParAlmondOptions);
  ins->vSolver        = vSolver;
  ins->vSolverOptions = vSolverOptions;

  if (rank==0) printf("==================PRESSURE SOLVE SETUP========================\n");
  // SETUP PRESSURE and VELOCITY SOLVERS
  dfloat zero =0.0;
  solver_t *pSolver   = ellipticSolveSetupQuad2D(mesh, ins->tau, zero, pBCType,kernelInfoP, pSolverOptions,pParAlmondOptions);
  ins->pSolver        = pSolver;
  ins->pSolverOptions = pSolverOptions;


  ins->VmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  ins->PmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  for (int e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) ins->VmapB[n+e*mesh->Np] = 1E9;
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[f+e*mesh->Nfaces];
      if (bc>0) {
        for (int n=0;n<mesh->Nfp;n++) {
          int fid = mesh->faceNodes[n+f*mesh->Nfp];
          ins->VmapB[fid+e*mesh->Np] = mymin(bc,ins->VmapB[fid+e*mesh->Np]);
          ins->PmapB[fid+e*mesh->Np] = mymax(bc,ins->PmapB[fid+e*mesh->Np]);
        }
      }
    }
  }
  gsParallelGatherScatter(mesh->hostGsh, ins->VmapB, "int", "min"); 
  gsParallelGatherScatter(mesh->hostGsh, ins->PmapB, "int", "max"); 

  for (int n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (ins->VmapB[n] == 1E9) {
      ins->VmapB[n] = 0.;
    }
  }
  ins->o_VmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), ins->VmapB);
  ins->o_PmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), ins->PmapB);

  kernelInfo.addDefine("p_blockSize", blockSize);
  kernelInfo.addParserFlag("automate-add-barriers", "disabled");


  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 256/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 256/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);


  int maxNodesVolumeCub = mymax(mesh->cubNp,mesh->Np);  
  kernelInfo.addDefine("p_maxNodesVolumeCub", maxNodesVolumeCub);
  int cubNblockV = mymax(1,256/maxNodesVolumeCub);
  //
  int maxNodesSurfaceCub = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  kernelInfo.addDefine("p_maxNodesSurfaceCub",maxNodesSurfaceCub);
  int cubNblockS = mymax(256/maxNodesSurfaceCub,1); // works for CUDA
  //
  kernelInfo.addDefine("p_cubNblockV",cubNblockV);
  kernelInfo.addDefine("p_cubNblockS",cubNblockS);

  if (rank==0) {
    printf("maxNodes: %d \t  NblockV: %d \t NblockS: %d  \n", maxNodes, NblockV, NblockS);
    printf("maxNodesVolCub: %d \t maxNodesSurCub: %d \t NblockVCub: %d \t NblockSCub: %d  \n", maxNodesVolumeCub,maxNodesSurfaceCub, cubNblockV, cubNblockS);

    printf("Np: %d \t Ncub: %d \n", mesh->Np, mesh->cubNp);
  }

  // ADD-DEFINES
  kernelInfo.addDefine("p_Lambda2", 0.5f);
  kernelInfo.addDefine("p_NTfields", ins->NTfields);
  kernelInfo.addDefine("p_NVfields", ins->NVfields);
  kernelInfo.addDefine("p_NfacesNfp",  mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_nu",      (float) ins->nu);
  //kernelInfo.addDefine("p_inu",      (float) 1.f/ins->nu);
  //kernelInfo.addDefine("p_idt",      (float) 1.f/ins->dt);


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
    ins->o_Ud   = mesh->device.malloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*sizeof(dfloat), ins->Ud);
    ins->o_Vd   = mesh->device.malloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*sizeof(dfloat), ins->Vd);
    ins->o_resU = mesh->device.malloc((mesh->Nelements)*mesh->Np*sizeof(dfloat), ins->resU);
    ins->o_resV = mesh->device.malloc((mesh->Nelements)*mesh->Np*sizeof(dfloat), ins->resV);

    for (int r=0;r<size;r++) {
      if (r==rank) {
        ins->subCycleVolumeKernel =
          mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",
    					 "insSubCycleVolume2D",
    					 kernelInfo);

        ins->subCycleSurfaceKernel =
          mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",
    					 "insSubCycleSurface2D",
    					 kernelInfo);

        ins->subCycleCubatureVolumeKernel =
          mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",
    					 "insSubCycleCubatureVolume2D",
    					 kernelInfo);

        ins->subCycleCubatureSurfaceKernel =
          mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",
    					 "insSubCycleCubatureSurface2D",
    					 kernelInfo);

        ins->subCycleRKUpdateKernel =
          mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",
    					 "insSubCycleRKUpdate2D",
    					 kernelInfo);

        ins->subCycleExtKernel =
          mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",
    					 "insSubCycleExt2D",
    					 kernelInfo);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  if(mesh->totalHaloPairs){//halo setup
    dlong tHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NTfields)*sizeof(dfloat);
    occa::memory o_tsendBuffer = mesh->device.mappedAlloc(tHaloBytes, NULL);
    occa::memory o_trecvBuffer = mesh->device.mappedAlloc(tHaloBytes, NULL);
    ins->o_tHaloBuffer = mesh->device.malloc(tHaloBytes);
    ins->tSendBuffer = (dfloat*) o_tsendBuffer.getMappedPointer();
    ins->tRecvBuffer = (dfloat*) o_trecvBuffer.getMappedPointer();

    dlong vHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NVfields)*sizeof(dfloat);
    occa::memory o_vsendBuffer = mesh->device.mappedAlloc(vHaloBytes, NULL);
    occa::memory o_vrecvBuffer = mesh->device.mappedAlloc(vHaloBytes, NULL);
    ins->o_vHaloBuffer = mesh->device.malloc(vHaloBytes);
    ins->vSendBuffer = (dfloat*) o_vsendBuffer.getMappedPointer();
    ins->vRecvBuffer = (dfloat*) o_vrecvBuffer.getMappedPointer();

    dlong pHaloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
    ins->o_pHaloBuffer = mesh->device.malloc(pHaloBytes);
    occa::memory o_psendBuffer = mesh->device.mappedAlloc(pHaloBytes, NULL);
    occa::memory o_precvBuffer = mesh->device.mappedAlloc(pHaloBytes, NULL);
    ins->pSendBuffer = (dfloat*) o_psendBuffer.getMappedPointer();
    ins->pRecvBuffer = (dfloat*) o_precvBuffer.getMappedPointer();
  }


  for (int r=0;r<size;r++) {
    if (r==rank) {
      // ===========================================================================
      ins->scaledAddKernel =
          mesh->device.buildKernelFromSource(DHOLMES "/okl/scaledAdd.okl",
               "scaledAddwOffset",
               kernelInfo);
     
      ins->totalHaloExtractKernel=
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
                 "insTotalHaloExtract2D",
                 kernelInfo);

      ins->totalHaloScatterKernel=
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
                 "insTotalHaloScatter2D",
                 kernelInfo);
      ins->velocityHaloExtractKernel=
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
                "insVelocityHaloExtract2D",
                 kernelInfo);

      ins->velocityHaloScatterKernel=
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
                 "insVelocityHaloScatter2D",
                 kernelInfo);

      ins->pressureHaloExtractKernel=
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
                 "insPressureHaloExtract",
                 kernelInfo);

      ins->pressureHaloScatterKernel=
       mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
                 "insPressureHaloScatter",
                 kernelInfo);  
     
      // ===========================================================================
      // ins->advectionCubatureVolumeKernel =
      //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection2D.okl",
    		// 		       "insAdvectionCubatureVolume2D",
    		// 		       kernelInfo);

      // ins->advectionCubatureSurfaceKernel =
      //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection2D.okl",
    		// 		       "insAdvectionCubatureSurface2D",
    		// 		       kernelInfo);

      ins->advectionVolumeKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvectionQuad2D.okl",
    				       "insAdvectionVolumeQuad2D",
    				       kernelInfo);

      ins->advectionSurfaceKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvectionQuad2D.okl",
    				       "insAdvectionSurfaceQuad2D",
    				       kernelInfo);

        // ins->subCycleRKUpdateKernel =
        //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",
        //        "insSubCycleRKUpdate2D",
        //        kernelInfo);
      // ===========================================================================
      ins->gradientVolumeKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insGradientQuad2D.okl",
    				       "insGradientVolumeQuad2D",
    				       kernelInfo);

      ins->gradientSurfaceKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insGradientQuad2D.okl",
    				       "insGradientSurfaceQuad2D",
    				       kernelInfo);

      // ===========================================================================
      ins->divergenceVolumeKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insDivergenceQuad2D.okl",
    				       "insDivergenceVolumeQuad2D",
    				       kernelInfo);

      ins->divergenceSurfaceKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insDivergenceQuad2D.okl",
    				       "insDivergenceSurfaceQuad2D",
    				       kernelInfo);

      // ===========================================================================
      ins->helmholtzRhsForcingKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insHelmholtzRhsQuad2D.okl",
    				       "insHelmholtzRhsForcingQuad2D",
    				       kernelInfo);

      ins->helmholtzRhsIpdgBCKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insHelmholtzRhsQuad2D.okl",
    				       "insHelmholtzRhsIpdgBCQuad2D",
    				       kernelInfo);

      // ===========================================================================
      ins->poissonRhsForcingKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonRhsQuad2D.okl",
    				       "insPoissonRhsForcingQuad2D",
    				       kernelInfo);

      ins->poissonRhsIpdgBCKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonRhsQuad2D.okl",
    				       "insPoissonRhsIpdgBCQuad2D",
    				       kernelInfo);

      ins->poissonRhsBCKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonRhsQuad2D.okl",
                   "insPoissonRhsBCQuad2D",
                   kernelInfo);

      ins->poissonAddBCKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonRhsQuad2D.okl",
                   "insPoissonAddBCQuad2D",
                   kernelInfo);

      // ins->poissonPenaltyKernel =
      //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonPenalty2D.okl",
    		// 		       "insPoissonPenalty2D",
    		// 		       kernelInfo);

      ins->updateUpdateKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insUpdate2D.okl",
    				       "insUpdateUpdate2D",
    				       kernelInfo);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  return ins;
}







