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

  if (size==1) options.getArgs("DEVICE NUMBER" ,deviceID);

  // read thread model/device/platform from options
  if(options.compareArgs("THREAD MODEL", "CUDA")){
    
    sprintf(deviceConfig, "mode = CUDA, deviceID = %d",deviceID);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenCL")){
    int plat;
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "mode = OpenCL, deviceID = %d, platformID = %d", deviceID, plat);
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

  mesh->Nfields = 1; 

  ins->g0 =  1.0;

  if (options.compareArgs("TIME INTEGRATOR", "ARK1")) {
    ins->Nstages = 1;
    int Nrk = 2;
    dfloat rkC[2] = {0.0, 1.0};
    dfloat erkA[2*2] ={0.0, 0.0,\
                       1.0, 0.0};
    dfloat irkA[2*2] ={0.0, 0.0,\
                       0.0, 1.0};
    dfloat prkA[2*2] ={0.0, 0.0,\
                       1.0, 0.0};

    ins->Nrk = Nrk;
    ins->rkC = (dfloat*) calloc(ins->Nrk, sizeof(dfloat));
    ins->erkA = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));
    ins->irkA = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));
    ins->prkA = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));

    memcpy(ins->rkC, rkC, ins->Nrk*sizeof(dfloat));
    memcpy(ins->erkA, erkA, ins->Nrk*ins->Nrk*sizeof(dfloat));
    memcpy(ins->irkA, irkA, ins->Nrk*ins->Nrk*sizeof(dfloat));
    memcpy(ins->prkA, prkA, ins->Nrk*ins->Nrk*sizeof(dfloat));

    ins->embeddedRKFlag = 0; //no embedded method
  } 


  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    ins->extbdfA = (dfloat*) calloc(3, sizeof(dfloat));
    ins->extbdfB = (dfloat*) calloc(3, sizeof(dfloat));
    ins->extbdfC = (dfloat*) calloc(3, sizeof(dfloat));

    ins->extC = (dfloat*) calloc(3, sizeof(dfloat));
  }

  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF1")) {
    ins->Nstages = 1;
    ins->temporalOrder = 1;
    ins->g0 = 1.0;
  } else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF2")) {
    ins->Nstages = 2;
    ins->temporalOrder = 2;
    ins->g0 = 1.5;
  } else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF3")) {
    ins->Nstages = 3;
    ins->temporalOrder = 3;
    ins->g0 = 11.f/6.f;
  }

  dlong Nlocal = mesh->Np*mesh->Nelements;
  dlong Ntotal = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  
  ins->Ntotal = Ntotal;
  ins->fieldOffset = Ntotal;
  ins->Nblock = (Nlocal+blockSize-1)/blockSize;

  // compute samples of q at interpolation nodes
  ins->U     = (dfloat*) calloc(ins->NVfields*ins->Nstages*Ntotal,sizeof(dfloat));
  ins->P     = (dfloat*) calloc(              ins->Nstages*Ntotal,sizeof(dfloat));

  //rhs storage
  ins->rhsU  = (dfloat*) calloc(Nlocal,sizeof(dfloat));
  ins->rhsV  = (dfloat*) calloc(Nlocal,sizeof(dfloat));
  ins->rhsW  = (dfloat*) calloc(Nlocal,sizeof(dfloat));
  ins->rhsP  = (dfloat*) calloc(Nlocal,sizeof(dfloat));

  //additional field storage
  ins->NU   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));
  ins->LU   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));
  ins->GP   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));

  ins->GU   = (dfloat*) calloc(ins->NVfields*Ntotal*4,sizeof(dfloat));
  
  ins->rkU  = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkP  = (dfloat*) calloc(              Ntotal,sizeof(dfloat));
  ins->PI   = (dfloat*) calloc(              Ntotal,sizeof(dfloat));
  
  ins->rkNU = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkLU = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkGP = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));

  //plotting fields
  ins->Vort = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->Div  = (dfloat*) calloc(              Nlocal,sizeof(dfloat));

  //extra storage for interpolated fields
  if(ins->elementType==HEXAHEDRA)
    ins->cU = (dfloat *) calloc(ins->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
  else 
    ins->cU = ins->U;

  ins->Nsubsteps = 0;
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF"))
    options.getArgs("SUBCYCLING STEPS",ins->Nsubsteps);

  if(ins->Nsubsteps){
    ins->Ud    = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
    ins->Ue    = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
    ins->resU  = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
    ins->rhsUd = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));

    if(ins->elementType==HEXAHEDRA)
      ins->cUd = (dfloat *) calloc(ins->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
    else 
      ins->cUd = ins->U;
  }

  dfloat rho  = 1.0 ;  // Give density for getting actual pressure in nondimensional solve
  dfloat g[3]; g[0] = 0.0; g[1] = 0.0; g[2] = 0.0;  // No gravitational acceleration

  options.getArgs("UBAR", ins->ubar);
  options.getArgs("VBAR", ins->vbar);
  if (ins->dim==3)
    options.getArgs("WBAR", ins->wbar);
  options.getArgs("PBAR", ins->pbar);
  options.getArgs("VISCOSITY", ins->nu);

  //Reynolds number
  ins->Re = ins->ubar/ins->nu;

  occa::kernelInfo kernelInfo;
  if(ins->dim==3)
    meshOccaSetup3D(mesh, deviceConfig, kernelInfo);
  else
    meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

  occa::kernelInfo kernelInfoV  = kernelInfo;
  occa::kernelInfo kernelInfoP  = kernelInfo;

  // ADD-DEFINES
  kernelInfo.addDefine("p_pbar", ins->pbar);
  kernelInfo.addDefine("p_ubar", ins->ubar);
  kernelInfo.addDefine("p_vbar", ins->vbar);
  kernelInfo.addDefine("p_wbar", ins->wbar);
  kernelInfo.addDefine("p_nu", ins->nu);

  kernelInfo.addDefine("p_NTfields", ins->NTfields);
  kernelInfo.addDefine("p_NVfields", ins->NVfields);
  kernelInfo.addDefine("p_NfacesNfp",  mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_Nstages",  ins->Nstages);
  kernelInfo.addDefine("p_SUBCYCLING",  ins->Nsubsteps);

  if (options.compareArgs("TIME INTEGRATOR", "ARK")) 
    ins->ARKswitch = 1;   
  else 
    ins->ARKswitch = 0;

  //add boundary data to kernel info
  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo.addInclude((char*)boundaryHeaderFileName.c_str());

  ins->o_U = mesh->device.malloc(ins->NVfields*ins->Nstages*Ntotal*sizeof(dfloat), ins->U);
  ins->o_P = mesh->device.malloc(              ins->Nstages*Ntotal*sizeof(dfloat), ins->P);

  for (int r=0;r<size;r++) {
    if (r==rank) {
      if (ins->dim==2) 
        ins->setFlowFieldKernel =  mesh->device.buildKernelFromSource(DINS "/okl/insSetFlowField2D.okl", "insSetFlowField2D", kernelInfo);  
      else
        ins->setFlowFieldKernel =  mesh->device.buildKernelFromSource(DINS "/okl/insSetFlowField3D.okl", "insSetFlowField3D", kernelInfo);  
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  ins->startTime =0.0;
  options.getArgs("START TIME", ins->startTime);
  ins->setFlowFieldKernel(mesh->Nelements,
                          ins->startTime,
                          mesh->o_x,
                          mesh->o_y,
                          mesh->o_z,
                          ins->fieldOffset,
                          ins->o_U,
                          ins->o_P);
  ins->o_U.copyTo(ins->U);

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
      const dlong id = n + mesh->Np*e;
      dfloat t = 0;
      dfloat uxn = ins->U[id+0*ins->fieldOffset];
      dfloat uyn = ins->U[id+1*ins->fieldOffset];
      dfloat uzn = 0.0;
      if (ins->dim==3) uzn = ins->U[id+2*ins->fieldOffset];


      //Squared maximum velocity
      dfloat numax;
      if (ins->dim==2)
        numax = uxn*uxn + uyn*uyn;
      else 
        numax = uxn*uxn + uyn*uyn + uzn*uzn;
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

  options.getArgs("FINAL TIME", ins->finalTime);
  
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    ins->NtimeSteps = ins->finalTime/ins->dt;

    if(ins->Nsubsteps){
      ins->dt         = ins->Nsubsteps*ins->dt;
      ins->NtimeSteps = ins->finalTime/ins->dt;
      ins->dt         = ins->finalTime/ins->NtimeSteps;
      ins->sdt        = ins->dt/ins->Nsubsteps;
    } else{
      ins->NtimeSteps = ins->finalTime/ins->dt;
      ins->dt         = ins->finalTime/ins->NtimeSteps;
    }
  }

  ins->dtMIN = 1E-2*ins->dt; //minumum allowed timestep

  if (rank==0) {
    printf("hmin = %g\n", hmin);
    printf("hmax = %g\n", hmax);
    printf("cfl = %g\n",  cfl);
    printf("dt = %g\n",   dt);
  }

  if (ins->Nsubsteps && rank==0) printf("dt: %.8f and sdt: %.8f ratio: %.8f \n", ins->dt, ins->sdt, ins->dt/ins->sdt);
  
  // Hold some inverses for kernels
  ins->inu = 1.0/ins->nu; 
  ins->idt = 1.0/ins->dt;
  
  ins->lambda = ins->g0 / (ins->dt * ins->nu);

  options.getArgs("TSTEPS FOR SOLUTION OUTPUT", ins->outputStep);
  if (rank==0) printf("Nsteps = %d NerrStep= %d dt = %.8e\n", ins->NtimeSteps,ins->outputStep, ins->dt);


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
  
  if (rank==0) printf("==================PRESSURE SOLVE SETUP=========================\n");
  ins->pSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  ins->pSolver->mesh = mesh;
  ins->pSolver->options = ins->pOptions;
  ins->pSolver->dim = ins->dim;
  ins->pSolver->elementType = ins->elementType;
  ins->pSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(ins->pSolver->BCType,pBCType,7*sizeof(int));
  ellipticSolveSetup(ins->pSolver, 0.0, kernelInfoP);


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

  // if (rank==0) {
  //   printf("maxNodes: %d \t  NblockV: %d \t NblockS: %d  \n", maxNodes, NblockV, NblockS);
  //   printf("maxNodesVolCub: %d \t maxNodesSurCub: %d \t NblockVCub: %d \t NblockSCub: %d  \n", maxNodesVolumeCub,maxNodesSurfaceCub, cubNblockV, cubNblockS);

  //   printf("Np: %d \t Ncub: %d \n", mesh->Np, mesh->cubNp);
  // }
  
  if (options.compareArgs("TIME INTEGRATOR", "ARK")) {
    ins->o_rkC  = mesh->device.malloc(         ins->Nrk*sizeof(dfloat),ins->rkC );
    ins->o_erkA = mesh->device.malloc(ins->Nrk*ins->Nrk*sizeof(dfloat),ins->erkA);
    ins->o_irkA = mesh->device.malloc(ins->Nrk*ins->Nrk*sizeof(dfloat),ins->irkA);
    ins->o_prkA = mesh->device.malloc(ins->Nrk*ins->Nrk*sizeof(dfloat),ins->prkA);
  }

  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    dfloat rkC[4] = {1.0, 0.0, -1.0, -2.0};

    ins->o_rkC  = mesh->device.malloc(4*sizeof(dfloat),rkC);
    ins->o_extbdfA = mesh->device.malloc(3*sizeof(dfloat));
    ins->o_extbdfB = mesh->device.malloc(3*sizeof(dfloat));
    ins->o_extbdfC = mesh->device.malloc(3*sizeof(dfloat)); 

    ins->o_extC = mesh->device.malloc(3*sizeof(dfloat)); 

    ins->o_prkA = ins->o_extbdfC;
  }

  // MEMORY ALLOCATION
  ins->o_rhsU  = mesh->device.malloc(Nlocal*sizeof(dfloat), ins->rhsU);
  ins->o_rhsV  = mesh->device.malloc(Nlocal*sizeof(dfloat), ins->rhsV);
  ins->o_rhsW  = mesh->device.malloc(Nlocal*sizeof(dfloat), ins->rhsW);
  ins->o_rhsP  = mesh->device.malloc(Nlocal*sizeof(dfloat), ins->rhsP);

  ins->o_NU    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->NU);
  ins->o_LU    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->LU);
  ins->o_GP    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->GP);
  
  ins->o_GU    = mesh->device.malloc(ins->NVfields*Ntotal*4*sizeof(dfloat), ins->GU);
  
  ins->o_rkU   = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkU);
  ins->o_rkP   = mesh->device.malloc(              Ntotal*sizeof(dfloat), ins->rkP);
  ins->o_PI    = mesh->device.malloc(              Ntotal*sizeof(dfloat), ins->PI);
  
  ins->o_rkNU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkNU);
  ins->o_rkLU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkLU);
  ins->o_rkGP  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkGP);

  //storage for helmholtz solves
  ins->o_UH = mesh->device.malloc(Ntotal*sizeof(dfloat));
  ins->o_VH = mesh->device.malloc(Ntotal*sizeof(dfloat));
  ins->o_WH = mesh->device.malloc(Ntotal*sizeof(dfloat));

  //plotting fields
  ins->o_Vort = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Vort);
  ins->o_Div  = mesh->device.malloc(              Nlocal*sizeof(dfloat), ins->Div);

  if(ins->elementType==HEXAHEDRA)
    ins->o_cU = mesh->device.malloc(ins->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->cU);
  else 
    ins->o_cU = ins->o_U;

  if(mesh->totalHaloPairs){//halo setup
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

    occa::memory o_gatherTmpPinned = mesh->device.mappedAlloc(ins->NVfields*mesh->ogs->NhaloGather*sizeof(dfloat), NULL);
    ins->velocityHaloGatherTmp = (dfloat*) o_gatherTmpPinned.getMappedPointer();
    ins->o_velocityHaloGatherTmp = mesh->device.malloc(ins->NVfields*mesh->ogs->NhaloGather*sizeof(dfloat),  ins->velocityHaloGatherTmp);
  }

  // set kernel name suffix
  char *suffix;
  
  if(ins->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(ins->elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(ins->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(ins->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  for (int r=0;r<size;r++) {
    if (r==rank) {
      sprintf(fileName, DINS "/okl/insHaloExchange.okl");
      sprintf(kernelName, "insVelocityHaloExtract");
      ins->velocityHaloExtractKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insVelocityHaloScatter");
      ins->velocityHaloScatterKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insPressureHaloExtract");
      ins->pressureHaloExtractKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insPressureHaloScatter");
      ins->pressureHaloScatterKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // ===========================================================================

      sprintf(fileName, DINS "/okl/insAdvection%s.okl", suffix);
      sprintf(kernelName, "insAdvectionCubatureVolume%s", suffix);
      ins->advectionCubatureVolumeKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insAdvectionCubatureSurface%s", suffix);
      ins->advectionCubatureSurfaceKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insAdvectionVolume%s", suffix);
      ins->advectionVolumeKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insAdvectionSurface%s", suffix);
      ins->advectionSurfaceKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // ===========================================================================
      
      sprintf(fileName, DINS "/okl/insDiffusion%s.okl", suffix);
      sprintf(kernelName, "insDiffusion%s", suffix);
      ins->diffusionKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(fileName, DINS "/okl/insDiffusionIpdg%s.okl", suffix);
      sprintf(kernelName, "insDiffusionIpdg%s", suffix);
      ins->diffusionIpdgKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(fileName, DINS "/okl/insVelocityGradient%s.okl", suffix);
      sprintf(kernelName, "insVelocityGradient%s", suffix);
      ins->velocityGradientKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // ===========================================================================

      sprintf(fileName, DINS "/okl/insGradient%s.okl", suffix);
      sprintf(kernelName, "insGradientVolume%s", suffix);
      ins->gradientVolumeKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insGradientSurface%s", suffix);
      ins->gradientSurfaceKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // ===========================================================================
      
      sprintf(fileName, DINS "/okl/insDivergence%s.okl", suffix);
      sprintf(kernelName, "insDivergenceVolume%s", suffix);
      ins->divergenceVolumeKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insDivergenceSurface%s", suffix);
      ins->divergenceSurfaceKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // ===========================================================================
      
      sprintf(fileName, DINS "/okl/insVelocityRhs%s.okl", suffix);
      if (options.compareArgs("TIME INTEGRATOR", "ARK")) 
        sprintf(kernelName, "insVelocityRhsARK%s", suffix);
      else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) 
        sprintf(kernelName, "insVelocityRhsEXTBDF%s", suffix);
      ins->velocityRhsKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(fileName, DINS "/okl/insVelocityBC%s.okl", suffix);
      sprintf(kernelName, "insVelocityIpdgBC%s", suffix);
      ins->velocityRhsIpdgBCKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insVelocityBC%s", suffix);
      ins->velocityRhsBCKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insVelocityAddBC%s", suffix);
      ins->velocityAddBCKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // ===========================================================================
      
      sprintf(fileName, DINS "/okl/insPressureRhs%s.okl", suffix);
      sprintf(kernelName, "insPressureRhs%s", suffix);
      ins->pressureRhsKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(fileName, DINS "/okl/insPressureBC%s.okl", suffix);
      sprintf(kernelName, "insPressureIpdgBC%s", suffix);
      ins->pressureRhsIpdgBCKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insPressureBC%s", suffix);
      ins->pressureRhsBCKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insPressureAddBC%s", suffix);
      ins->pressureAddBCKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // ===========================================================================

      sprintf(fileName, DINS "/okl/insPressureUpdate.okl");
      sprintf(kernelName, "insPressureUpdate");
      ins->pressureUpdateKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      sprintf(fileName, DINS "/okl/insVelocityUpdate.okl");
      sprintf(kernelName, "insVelocityUpdate");
      ins->velocityUpdateKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);      

      // ===========================================================================

      sprintf(fileName, DINS "/okl/insVorticity%s.okl", suffix);
      sprintf(kernelName, "insVorticity%s", suffix);
      ins->vorticityKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);
    

      if(ins->Nsubsteps){
        // Note that resU and resV can be replaced with already introduced buffer
        ins->o_Ue    = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Ue);
        ins->o_Ud    = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Ud);
        ins->o_resU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->resU);
        ins->o_rhsUd = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rhsUd);

        if(ins->elementType==HEXAHEDRA)
          ins->o_cUd = mesh->device.malloc(ins->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->cUd);
        else 
          ins->o_cUd = ins->o_Ud;

        sprintf(fileName, DHOLMES "/okl/scaledAdd.okl");
        sprintf(kernelName, "scaledAddwOffset");
        ins->scaledAddKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

        sprintf(fileName, DINS "/okl/insSubCycle%s.okl", suffix);
        sprintf(kernelName, "insSubCycleVolume%s", suffix);
        ins->subCycleVolumeKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "insSubCycleSurface%s", suffix);
        ins->subCycleSurfaceKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "insSubCycleCubatureVolume%s", suffix);
        ins->subCycleCubatureVolumeKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "insSubCycleCubatureSurface%s", suffix);
        ins->subCycleCubatureSurfaceKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

        sprintf(fileName, DINS "/okl/insSubCycle.okl");
        sprintf(kernelName, "insSubCycleRKUpdate");
        ins->subCycleRKUpdateKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "insSubCycleExt");
        ins->subCycleExtKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  return ins;
}







