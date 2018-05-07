#include "ins.h"
#include "omp.h"
#include <unistd.h>

ins_t *insSetup(mesh_t *mesh, setupAide options){

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  long int hostId = gethostid();

  long int* hostIds = (long int*) calloc(size,sizeof(long int));
  MPI_Allgather(&hostId,1,MPI_LONG,hostIds,1,MPI_LONG,MPI_COMM_WORLD);

  int deviceID = 0;
  int totalDevices = 0;
  for (int r=0;r<rank;r++) {
    if (hostIds[r]==hostId) deviceID++;
  }
  for (int r=0;r<size;r++) {
    if (hostIds[r]==hostId) totalDevices++;
  }

  // read thread model/device/platform from options
  if(options.compareArgs("THREAD MODEL", "CUDA")){
    int dev;
    options.getArgs("DEVICE NUMBER" ,dev);
    sprintf(deviceConfig, "mode = CUDA, deviceID = %d",dev);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenCL")){
    int dev, plat;
    options.getArgs("DEVICE NUMBER", dev);
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "mode = OpenCL, deviceID = %d, platformID = %d", dev, plat);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenMP")){
    sprintf(deviceConfig, "mode = OpenMP");
  }
  else{
    sprintf(deviceConfig, "mode = Serial");
  }
        
  //set number of omp threads to use
  int Ncores = sysconf(_SC_NPROCESSORS_ONLN);
  int Nthreads = Ncores/totalDevices;
  omp_set_num_threads(Nthreads);
  if (rank==0 && options.compareArgs("VERBOSE","TRUE")) 
    printf("Rank %d: Ncores = %d, Nthreads = %d\n", rank, Ncores, Nthreads);

  ins_t *ins = (ins_t*) calloc(1, sizeof(ins_t));
  ins->mesh = mesh;
  ins->options = options;

  options.getArgs("MESH DIMENSION", ins->dim);
  options.getArgs("ELEMENT TYPE", ins->elementType);

  ins->NVfields = (ins->dim==3) ? 3:2; //  Total Number of Velocity Fields
  ins->NTfields = (ins->dim==3) ? 4:3; // Total Velocity + Pressure
  ins->Nfields  = 1; // Each Velocity Field

  mesh->Nfields = ins->Nfields;

  dlong Ntotal = mesh->Np*mesh->Nelements;
  ins->Nblock = (Ntotal+blockSize-1)/blockSize;

  int Nstages = 3;

  // compute samples of q at interpolation nodes
  ins->U     = (dfloat*) calloc(ins->NVfields*Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->P     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));

  //rhs storage
  ins->rhsU  = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  ins->rhsV  = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  ins->rhsW  = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  ins->rhsP  = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));

  //additional field storage (could reduce in the future)
  ins->NU     = (dfloat*) calloc(ins->NVfields*Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->gradP  = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->PI     = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));

  ins->Vort     = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  ins->Div     = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));

  if(ins->dim==3) {
    ins->Vx     = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
    ins->Vy     = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
    ins->Vz     = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  }

  options.getArgs("SUCYCLING STEPS",ins->Nsubsteps)

  if(ins->Nsubsteps){
    ins->Ud   = (dfloat*) calloc(ins->NVfields*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
    ins->Ue   = (dfloat*) calloc(ins->NVfields*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
    ins->resU = (dfloat*) calloc(ins->NVfields*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  }

  dfloat rho  = 1.0  ;  // Give density for getting actual pressure in nondimensional solve
  dfloat g[2]; g[0] = 0.0; g[1] = 0.0;  // No gravitational acceleration

  options.getArgs("UBAR", cns->ubar);
  options.getArgs("VBAR", cns->vbar);
  if (ins->dim==3)
    options.getArgs("WBAR", cns->wbar);
  options.getArgs("PBAR", cns->pbar);
  options.getArgs("VISCOSITY", cns->nu);

  //Reynolds number
  ins->Re = ins->ubar/ins->nu;

  // Define total DOF per field for INS i.e. (Nelelemts + Nelements_halo)*Np
  ins->NtotalDofs = (mesh->totalHaloPairs+mesh->Nelements)*mesh->Np ;
  ins->NDofs      = mesh->Nelements*mesh->Np;
  
  // Initialize
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      const int id = n + mesh->Np*e;
      dfloat t = 0;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];
      dfloat z = mesh->z[id];
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

      #if 1 // mean flow
        ins->U[id+0*ins->NtotalDofs] = ins->ubar;
        ins->U[id+1*ins->NtotalDofs] = ins->vbar;
        if (ins->dim==3)
          ins->U[id+2*ins->NtotalDofs] = ins->wbar;

        ins->P[id] = pbar;
      #endif
    }
  }

  // set time step
  dfloat hmin = 1e9, hmax = 0;
  dfloat umax = 0;
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      dfloat hest = 2./(sJ*invJ);

      hmin = mymin(hmin, hest);
      hmax = mymax(hmax, hest);
    }
    for(int n=0;n<mesh->Np;++n){
      const int id = n + mesh->Np*e;
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
  dfloat magVel = mymax(umax,1.0); // Correction for initial zero velocity
 
  dfloat cfl; 
  options.getArgs("CFL", cfl);
  dfloat dt     = cfl* hmin/( (mesh->N+1.)*(mesh->N+1.) * magVel) ;

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(ins->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  options.getArgs("FINAL TIME", mesh->finalTime);
  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  if (options.compareArgs("TIME INTEGRATOR","LSERK4")){
    mesh->dt = mesh->finalTime/mesh->NtimeSteps;
  }

  if(ins->Nsubsteps){
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
    printf("dt = %g\n",   dt);
  }

  if (strstr(options,"SUBCYCLING")&&rank==0) printf("dt: %.8f and sdt: %.8f ratio: %.8f \n", ins->dt, ins->sdt, ins->dt/ins->sdt);
  
  // Hold some inverses for kernels
  ins->inu = 1.0/ins->nu; 
  ins->idt = 1.0/ins->dt;
  
  options.getArgs("TSTEPS FOR SOLUTION OUTPUT", ins->outputStep);
  if (rank==0) printf("Nsteps = %d NerrStep= %d dt = %.8e\n", ins->NtimeSteps,ins->outputStep, ins->dt);

  occa::kernelInfo kernelInfo;
  if(cns->dim==3)
    meshOccaSetup3D(mesh, deviceConfig, kernelInfo);
  else
    meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

  //make option objects for elliptc solvers
  ins->vOptions = options;
  ins->vOptions.setArgs("KRYLOV SOLVER",        options.getArgs("VELOCITY KRYLOV SOLVER"));
  ins->vOptions.setArgs("DISCRETIZATION",       options.getArgs("VELOCITY DISCRETIZATION"));
  ins->vOptions.setArgs("BASIS",                options.getArgs("VELOCITY BASIS"));
  ins->vOptions.setArgs("PRECONDITIONER",       options.getArgs("VELOCITY PRECONDITIONER"));
  ins->vOptions.setArgs("MULTIGRID COARSENING", options.getArgs("VELOCITY MULTIGRID COARSENING"));
  ins->vOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("VELOCITY MULTIGRID SMOOTHER"));
  ins->vOptions.setArgs("PARALMOND CYCLE",      options.getArgs("VELOCITY PARALMOND CYCLE"));
  ins->vOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("VELOCITY PARALMOND SMOOTHER"));
  ins->vOptions.setArgs("PARALMOND PARTITION",  options.getArgs("VELOCITY PARALMOND PARTITION"));

  ins->pOptions = options; 
  ins->pOptions.setArgs("KRYLOV SOLVER",        options.getArgs("PRESSURE KRYLOV SOLVER"));
  ins->pOptions.setArgs("DISCRETIZATION",       options.getArgs("PRESSURE DISCRETIZATION"));
  ins->pOptions.setArgs("BASIS",                options.getArgs("PRESSURE BASIS"));
  ins->pOptions.setArgs("PRECONDITIONER",       options.getArgs("PRESSURE PRECONDITIONER"));
  ins->pOptions.setArgs("MULTIGRID COARSENING", options.getArgs("PRESSURE MULTIGRID COARSENING"));
  ins->pOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("PRESSURE MULTIGRID SMOOTHER"));
  ins->pOptions.setArgs("PARALMOND CYCLE",      options.getArgs("PRESSURE PARALMOND CYCLE"));
  ins->pOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("PRESSURE PARALMOND SMOOTHER"));
  ins->pOptions.setArgs("PARALMOND PARTITION",  options.getArgs("PRESSURE PARALMOND PARTITION"));

  occa::kernelInfo kernelInfoV  = kernelInfo;
  occa::kernelInfo kernelInfoP  = kernelInfo;

  if (rank==0) printf("==================ELLIPTIC SOLVE SETUP=========================\n");

  // SetUp Boundary Flags types for Elliptic Solve
  // bc = 1 -> wall
  // bc = 2 -> inflow
  // bc = 3 -> outflow
  // bc = 4 -> x-aligned slip
  // bc = 5 -> y-aligned slip
  // bc = 6 -> z-aligned slip
  int uBCType[7] = {0,1,1,2,1,2,2}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int vBCType[7] = {0,1,1,2,2,1,2}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int wBCType[7] = {0,1,1,2,2,2,1}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int pBCType[7] = {0,2,2,1,2,2,2}; // bc=3 => outflow => Dirichlet => pBCType[3] = 1, etc.

  //Solver tolerances 
  ins->presTOL = 1E-8;
  ins->velTOL  = 1E-8;

  // Use third Order Velocity Solve: full rank should converge for low orders
  if (rank==0) printf("==================VELOCITY SOLVE SETUP=========================\n");
  // ins->lambda = (11.f/6.f) / (ins->dt * ins->nu);
  ins->lambda = (1.5f) / (ins->dt * ins->nu);

  ins->uSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  ins->uSolver->mesh = mesh;
  ins->uSolver->options = ins->vOptions;
  ins->uSolver->dim = ins->dim;
  ins->uSolver->elementType = ins->elementType;
  ins->uSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(ins->uSolver->BCType,uBCType,7*sizeof(int));
  ellipticSolveSetup(ins->uSolver, ins->lambda, kernelInfoV);

  ins->vSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  ins->vSolver->mesh = mesh;
  ins->vSolver->options = ins->vOptions;
  ins->vSolver->dim = ins->dim;
  ins->vSolver->elementType = ins->elementType;
  ins->vSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(ins->vSolver->BCType,vBCType,7*sizeof(int));
  ellipticSolveSetup(ins->vSolver, ins->lambda, kernelInfoV);

  if (ins->dim==3) {
    ins->wSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
    ins->wSolver->mesh = mesh;
    ins->wSolver->options = ins->vOptions;
    ins->wSolver->dim = ins->dim;
    ins->wSolver->elementType = ins->elementType;
    ins->wSolver->BCType = (int*) calloc(7,sizeof(int));
    memcpy(ins->wSolver->BCType,wBCType,7*sizeof(int));
    ellipticSolveSetup(ins->wSolver, ins->lambda, kernelInfoV);  
  }
  
  if (rank==0) printf("==================PRESSURE SOLVE SETUP========================\n");
  ins->pSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  ins->pSolver->mesh = mesh;
  ins->pSolver->options = ins->vOptions;
  ins->pSolver->dim = ins->dim;
  ins->pSolver->elementType = ins->elementType;
  ins->pSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(ins->pSolver->BCType,pBCType,7*sizeof(int));
  ellipticSolveSetup(ins->pSolver, 0.0, kernelInfoV);


  //make node-wise boundary flags
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

  if (rank==0 && options.compareArgs("VERBOSE","TRUE")) 
    occa::setVerboseCompilation(true);
  else 
    occa::setVerboseCompilation(false);

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

    // capture header file
  char *boundaryHeaderFileName;
  if(argc==3)
    boundaryHeaderFileName = strdup(DHOLMES "/examples/insTri2D/insUniform2D.h"); // default
  else
    boundaryHeaderFileName = strdup(argv[3]);

  //add boundary data to kernel info
  kernelInfo.addInclude(boundaryHeaderFileName);
  
  // MEMORY ALLOCATION

  ins->o_U = mesh->device.malloc(ins->NVfields*Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->U);
  ins->o_P = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->P);

  ins->o_rhsU  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsU);
  ins->o_rhsV  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsV);
  ins->o_rhsW  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsW);
  ins->o_rhsP  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsP);

  ins->o_NU    = mesh->device.malloc(ins->NVfields*Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->NU);
  ins->o_gradP = mesh->device.malloc(ins->dim*Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->Px);
  ins->o_PI = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->PI);
  ins->o_gradPx = mesh->device.malloc(ins->dim*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat));

  //storage for helmholtz solves. Fix this later
  ins->o_UH = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat));
  ins->o_VH = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat));
  ins->o_WH = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat));

  ins->o_Vort = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), ins->Vort);
  ins->o_Div = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), ins->Div);

  if (ins->dim==3) {
    ins->o_Vx = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), ins->Vx);
    ins->o_Vy = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), ins->Vy);
    ins->o_Vz = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), ins->Vz);
  }

  if(strstr(options,"SUBCYCLING")){
    // Note that resU and resV can be replaced with already introduced buffer
    ins->o_Ue   = mesh->device.malloc(ins->NVfields*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*sizeof(dfloat), ins->Ue);
    ins->o_Ud   = mesh->device.malloc(ins->NVfields*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*sizeof(dfloat), ins->Ud);
    ins->o_resU = mesh->device.malloc(ins->NVfields*(mesh->Nelements)*mesh->Np*sizeof(dfloat), ins->resU);

    for (int r=0;r<size;r++) {
      if (r==rank) {
        ins->subCycleVolumeKernel =
          mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycleTri2D.okl",
    					 "insSubCycleVolumeTri2D",
    					 kernelInfo);

        ins->subCycleSurfaceKernel =
          mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycleTri2D.okl",
    					 "insSubCycleSurfaceTri2D",
    					 kernelInfo);

        ins->subCycleCubatureVolumeKernel =
          mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycleTri2D.okl",
    					 "insSubCycleCubatureVolumeTri2D",
    					 kernelInfo);

        ins->subCycleCubatureSurfaceKernel =
          mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycleTri2D.okl",
    					 "insSubCycleCubatureSurfaceTri2D",
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
      ins->advectionCubatureVolumeKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvectionTri2D.okl",
    				       "insAdvectionCubatureVolumeTri2D",
    				       kernelInfo);

      ins->advectionCubatureSurfaceKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvectionTri2D.okl",
    				       "insAdvectionCubatureSurfaceTri2D",
    				       kernelInfo);

      ins->advectionVolumeKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvectionTri2D.okl",
    				       "insAdvectionVolumeTri2D",
    				       kernelInfo);

      ins->advectionSurfaceKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvectionTri2D.okl",
    				       "insAdvectionSurfaceTri2D",
    				       kernelInfo);

      // ===========================================================================
      ins->gradientVolumeKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insGradientTri2D.okl",
    				       "insGradientVolumeTri2D",
    				       kernelInfo);

      ins->gradientSurfaceKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insGradientTri2D.okl",
    				       "insGradientSurfaceTri2D",
    				       kernelInfo);

      // ===========================================================================
      ins->divergenceVolumeKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insDivergenceTri2D.okl",
    				       "insDivergenceVolumeTri2D",
    				       kernelInfo);

      ins->divergenceSurfaceKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insDivergenceTri2D.okl",
    				       "insDivergenceSurfaceTri2D",
    				       kernelInfo);

      // ===========================================================================
      ins->helmholtzRhsForcingKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insHelmholtzRhsTri2D.okl",
    				       "insHelmholtzRhsForcingTri2D",
    				       kernelInfo);

      ins->helmholtzRhsIpdgBCKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insHelmholtzRhsTri2D.okl",
    				       "insHelmholtzRhsIpdgBCTri2D",
    				       kernelInfo);

      ins->helmholtzRhsBCKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insHelmholtzRhsTri2D.okl",
                   "insHelmholtzRhsBCTri2D",
                   kernelInfo);

      ins->helmholtzAddBCKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insHelmholtzRhsTri2D.okl",
                   "insHelmholtzAddBCTri2D",
                   kernelInfo);

      // ===========================================================================
      ins->poissonRhsForcingKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonRhsTri2D.okl",
    				       "insPoissonRhsForcingTri2D",
    				       kernelInfo);

      ins->poissonRhsIpdgBCKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonRhsTri2D.okl",
    				       "insPoissonRhsIpdgBCTri2D",
    				       kernelInfo);

      ins->poissonRhsBCKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonRhsTri2D.okl",
                   "insPoissonRhsBCTri2D",
                   kernelInfo);

      ins->poissonAddBCKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonRhsTri2D.okl",
                   "insPoissonAddBCTri2D",
                   kernelInfo);

      // ins->poissonPenaltyKernel =
      //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonPenaltyTri2D.okl",
    		// 		       "insPoissonPenaltyTri2D",
    		// 		       kernelInfo);

      ins->updateUpdateKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insUpdate2D.okl",
    				       "insUpdateUpdate2D",
    				       kernelInfo);

      ins->vorticityKernel=
        mesh->device.buildKernelFromSource(DHOLMES "/okl/insVorticityTri2D.okl",
                   "insVorticityTri2D",
                   kernelInfo);     
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  return ins;
}







