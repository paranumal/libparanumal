#include "ins2D.h"

// NBN: toggle use of 2nd stream
#define USE_2_STREAMS

ins_t *insSetup2D(mesh2D *mesh, char * options, char *vSolverOptions, char *pSolverOptions){

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  // use rank to choose DEVICE
  // sprintf(deviceConfig, "mode = CUDA, deviceID = %d", (rank+1)%3);
  //sprintf(deviceConfig, "mode = CUDA, deviceID = 0");
  // sprintf(deviceConfig, "mode = OpenCL, deviceID = 1, platformID = 0");
  // sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 0");
  
  sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  // sprintf(deviceConfig, "mode = Serial");  

  // compute samples of q at interpolation nodes
 
  ins_t *ins = (ins_t*) calloc(1, sizeof(ins_t));
  
  ins->NVfields = 2; //  Total Number of Velocity Fields 
  ins->NTfields = 3; // Total Velocity + Pressure
  ins->Nfields  = 1; // Each Velocity Field
  ins->ExplicitOrder = 3; // Order Nonlinear Extrapolation
  
  mesh->Nfields = ins->NTfields;
  
  ins->mesh = mesh;  
  // compute samples of q at interpolation nodes
  ins->U     = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->V     = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->P     = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->PI    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*ins->Nfields,sizeof(dfloat));
  
  //
  ins->rhsU  = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->rhsV  = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->rhsP  = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));
  //
  ins->NU    = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->NV    = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));  //
  
  ins->UO     = (dfloat*) calloc(mesh->Nelements*mesh->Np*(ins->ExplicitOrder-1)*ins->NVfields,sizeof(dfloat));
  ins->NO     = (dfloat*) calloc(mesh->Nelements*mesh->Np*(ins->ExplicitOrder-1)*ins->NVfields,sizeof(dfloat));
  
  // Hold UO in [uxn uyn ux(n-1) uy(n-1) ] format
    
  // SET SOLVER OPTIONS 
  // Initial Conditions, Flow Properties 
  printf("Starting initial conditions for INS2D\n");
  // 
  dfloat ux   = 0.0  ; 
  dfloat uy   = 0.0  ; 
  dfloat pr   = 0.0  ;   
  dfloat nu   = 1.0/40.0;  // kinematic viscosity,
  dfloat rho  = 1.0  ;  // Give density for getting actual pressure in nondimensional solve
  
  dfloat g[2]; g[0] = 0.0; g[1] = 0.0;  // No gravitational acceleration
  
#if 0
  // Create Periodic Boundaries
  printf("Creating periodic connections if exist \n");
  dfloat xper = 1.0;   dfloat yper = 0.0;
  insMakePeriodic2D(mesh,xper,yper);
#endif
  // Fill up required fileds
  ins->finalTime = 0.01;
  ins->nu        = nu ;
  ins->rho       = rho;
  ins->tau       = 2.f*(mesh->N+1)*(mesh->N+1); 
  //
  //memcpy(ins->g, g,2*sizeof(dfloat));

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

      dfloat lamda = 1./(2. * ins->nu) - sqrt(1./(4.*ins->nu * ins->nu) + 4.*M_PI*M_PI) ;  
      //
      ins->U[id] = 1.0 - exp(lamda*x)*cos(2.*M_PI*y); 
      ins->V[id] = lamda/(2.*M_PI)*exp(lamda*x)*sin(2.*M_PI*y);
      ins->P[id] = 0.5*(1.0- exp(2.*lamda*x));
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

  dfloat cfl = 0.25; 
  dfloat magVel = mymax(umax,1.0); // Correction for initial zero velocity
  dfloat dt = cfl* hmin/( (mesh->N+1.)*(mesh->N+1.) * magVel) ;
  
  printf("hmin = %g\n", hmin);
  printf("hmax = %g\n", hmax);
  printf("cfl = %g\n", cfl);
  printf("dt = %g\n", dt);
   
  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(ins->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
 
  ins->NtimeSteps = ins->finalTime/ins->dt;
  ins->dt         = ins->finalTime/ins->NtimeSteps;

  // errorStep
  ins->errorStep =100;

  printf("Nsteps = %d NerrStep= %d dt = %.8e\n", ins->NtimeSteps,ins->errorStep, ins->dt);
  
   // specialization for Boltzmann

  #if 1

  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields, sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields, sizeof(dfloat));
  mesh->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields, sizeof(dfloat));
  

  occa::kernelInfo kernelInfo;
  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

  occa::kernelInfo kernelInfoV  = kernelInfo;
  occa::kernelInfo kernelInfoP  = kernelInfo;


   printf("==================ELLIPTIC SOLVE SETUP=========================\n");
  // SetUp Boundary Conditions Flags for Elliptic Solve
 
  iint *vEToB = (iint*) calloc(mesh->Nelements*mesh->Nfaces, sizeof(iint));
  iint *pEToB = (iint*) calloc(mesh->Nelements*mesh->Nfaces, sizeof(iint));
  //
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint f=0;f<mesh->Nfaces;++f){
      const iint id = f + mesh->Nfaces*e; 
      const iint bc = mesh->EToB[id];

      if(bc==1 || bc == 2) {// Wall or Inflow
       vEToB[id] = 1;  // Drichlet for Velocity
       pEToB[id] = 2;  // Neumann for Pressure
      }

      if(bc==3) {// OutFlow
       vEToB[id] = 2; // Neumann for Velocity
       pEToB[id] = 1; // Dirichlet for Pressure
      }
    }
  }
  // Use third Order Velocity Solve: full rank should converge for low orders
   printf("==================VELOCITY SOLVE SETUP=========================\n");
  ins->lambda = (1.5f) / (ins->dt * ins->nu);
  solver_t *vSolver   = ellipticSolveSetupTri2D(mesh, ins->tau, ins->lambda, vEToB, kernelInfoV, vSolverOptions); 
  ins->vSolver        = vSolver;  
  ins->vSolverOptions = vSolverOptions;
  
  printf("==================PRESSURE SOLVE SETUP========================\n");
  // SETUP PRESSURE and VELOCITY SOLVERS
  solver_t *pSolver   = ellipticSolveSetupTri2D(mesh, ins->tau, 0.0, pEToB,kernelInfoP, pSolverOptions); 
  ins->pSolver        = pSolver; 
  ins->pSolverOptions = pSolverOptions;


  free(vEToB);
  free(pEToB);
  #endif

  kernelInfo.addDefine("p_maxNodesVolume", mymax(mesh->cubNp,mesh->Np));
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);


  int maxSurfaceNodes = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  kernelInfo.addDefine("p_maxSurfaceNodes", maxSurfaceNodes);
  printf("maxSurfaceNodes=%d\n", maxSurfaceNodes);


  #if 1
  int NblockV = 256/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 128/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  printf("maxNodes: %d \t NblockV: %d \t NblockS: %d  \n", maxNodes, NblockV, NblockS);
  
  

  ins->PID  = 0; 
  ins->PIID = 1; 
  #endif
  
  printf("Np: %d \t Ncub: %d \n", mesh->Np, mesh->cubNp);



  // ADD-DEFINES 
  kernelInfo.addDefine("p_Lambda2", 0.5f);
  kernelInfo.addDefine("p_NTfields", ins->NTfields);
  kernelInfo.addDefine("p_NVfields", ins->NVfields);
  kernelInfo.addDefine("p_NfacesNfp",  mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_inu",      (float) 1.f/ins->nu);
  kernelInfo.addDefine("p_nu",      (float) ins->nu);
  kernelInfo.addDefine("p_idt",      (float) 1.f/ins->dt);

  // MEMORY ALLOCATION
  ins->o_U  =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*ins->Nfields*sizeof(dfloat), ins->U);
  ins->o_V  =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*ins->Nfields*sizeof(dfloat), ins->V);  
  ins->o_P =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat),  ins->P);  
  ins->o_PI =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->PI);      


  ins->o_rhsU  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsU);
  ins->o_rhsV  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsV);
  ins->o_rhsP  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsP); 
  ins->o_NU    = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->NU); 
  ins->o_NV    = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->NV);  
  //  
  ins->o_UO     = mesh->device.malloc(mesh->Nelements*mesh->Np*(ins->ExplicitOrder-1)*ins->NVfields*sizeof(dfloat),ins->UO);
  ins->o_NO     = mesh->device.malloc(mesh->Nelements*mesh->Np*(ins->ExplicitOrder-1)*ins->NVfields*sizeof(dfloat),ins->NO);

  
  // Need to be changed !!!!!
  if(mesh->totalHaloPairs){
    ins->o_totHaloBuffer = mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*(ins->NTfields)*sizeof(dfloat));
    ins->o_velHaloBuffer = mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*(ins->NVfields)*sizeof(dfloat));
    ins->o_prHaloBuffer  = mesh->device.malloc(mesh->totalHaloPairs*mesh->Np *sizeof(dfloat));
  }

 // =========================================================================== 
  if(1) { // strstr(options, "CUBATURE")){ 
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
  }
  if(1) { // else{
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
  }
  
// =========================================================================== 
  printf("Compiling Gradient volume kernel with collocation integration\n");
  ins->gradientVolumeKernel = 
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insGradient2D.okl",
      "insGradientVolume2D",
        kernelInfo);
   
     // KERNEL DEFINITIONS
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
  
  
  // printf("Compiling INS Helmholtz Halo Extract Kernel\n");
  // ins->helmholtzHaloExtractKernel= 
  //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
  //     "insHelmholtzHaloExtract2D",
  //       kernelInfo);

  //  printf("Compiling INS Helmholtz Halo Extract Kernel\n");
  // ins->helmholtzHaloScatterKernel= 
  //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
  //     "insHelmholtzHaloScatter2D",
  //       kernelInfo);  


 // ===========================================================================
  

  printf("Compiling Helmholtz volume update kernel\n");
  ins->poissonRhsForcingKernel = 
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonRhs2D.okl",
      "insPoissonRhsForcing2D",
        kernelInfo);  

  // printf("Compiling Poisson IPDG surface kernel with collocation integration\n");
  // ins->poissonRhsIpdgBCKernel = 
  //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonRhs2D.okl",
  //     "insPoissonRhsIpdgBC2D",
  //       kernelInfo);




  //   printf("Compiling INS Poisson Halo Extract Kernel\n");
  // ins->poissonHaloExtractKernel= 
  //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
  //     "insPoissonHaloExtract2D",
  //       kernelInfo);

  //  printf("Compiling INS PoissonHalo Extract Kernel\n");
  // ins->poissonHaloScatterKernel= 
  //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
  //     "insPoissonHaloScatter2D",
  //       kernelInfo);


// ===========================================================================//
  // printf("Compiling Poisson volume kernel with collocation integration\n");
  // ins->updateVolumeKernel = 
  //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insUpdate2D.okl",
  //     "insUpdateVolume2D",
  //       kernelInfo);

  // printf("Compiling Poisson surface kernel with collocation integration\n");
  // ins->updateSurfaceKernel = 
  //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insUpdate2D.okl",
  //     "insUpdateSurface2D",
  //       kernelInfo);


  printf("Compiling Poisson surface kernel with collocation integration\n");
  ins->updateUpdateKernel = 
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insUpdate2D.okl",
      "insUpdateUpdate2D",
        kernelInfo);


 // printf("Compiling INS Poisson Halo Extract Kernel\n");
 //  ins->updateHaloExtractKernel= 
 //    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
 //      "insUpdateHaloExtract2D",
 //        kernelInfo);

 //   printf("Compiling INS PoissonHalo Extract Kernel\n");
 //  ins->updateHaloScatterKernel= 
 //    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
 //      "insUpdateHaloScatter2D",
 //        kernelInfo); 
// ===========================================================================//


 // INS->mesh = mesh; // No modificatin of mesh after this point 
  ins->mesh       = mesh;

  return ins;

}


  




