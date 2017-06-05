#include "ins2D.h"

// NBN: toggle use of 2nd stream
#define USE_2_STREAMS

ins_t *insSetup2D(mesh2D *mesh, char * options, char *velSolverOptions, char *prSolverOptions){

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
  ins->Ux     = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->Uy     = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->Pr     = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->PrI    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*ins->Nfields,sizeof(dfloat));
  
  //
  ins->rhsUx  = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->rhsUy  = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->rhsPr  = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));
  //
  ins->NUx    = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->NUy    = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));  //
  
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
  ins->finalTime = 1. ;
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
      ins->Ux[id] = 1.0 - exp(lamda*x)*cos(2.*M_PI*y); 
      ins->Uy[id] = lamda/(2.*M_PI)*exp(lamda*x)*sin(2.*M_PI*y);
      ins->Pr[id] = 0.5*(1.0- exp(2.*lamda*x));
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
      dfloat uxn = ins->Ux[id];
      dfloat uyn = ins->Uy[id];

      //Squared maximum velocity
      dfloat numax = uxn*uxn + uyn*uyn; 
      umax = mymax(umax, numax);
     
    }
  }
  // Maximum Velocity
  umax = sqrt(umax);

  dfloat cfl = 0.5; 
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
  ins->errorStep = 100;

  printf("Nsteps = %d NerrStep= %d dt = %.8e\n", ins->NtimeSteps,ins->errorStep, ins->dt);
  
   // specialization for Boltzmann

  #if 1

  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields, sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields, sizeof(dfloat));
  mesh->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields, sizeof(dfloat));
  

  occa::kernelInfo kernelInfo;
  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

  occa::kernelInfo kernelInfoVel = kernelInfo;
  occa::kernelInfo kernelInfoPr  = kernelInfo;


   printf("==================PRESSURE SOLVE SETUP========================\n");
  // SETUP PRESSURE and VELOCITY SOLVERS
  solver_t *prsolver   = ellipticSolveSetupTri2D(mesh,0.0, kernelInfoPr, prSolverOptions); 
  ins->prsolver        = prsolver; 
  ins->prsolverOptions = prSolverOptions;

  // Use third Order Velocity Solve: full rank should converge for low orders
  //ins->lamda = (11./ 6.) / (ins->dt * ins->nu);
   printf("==================VELOCITY SOLVE SETUP=========================\n");
  ins->lamda = (1.0) / (ins->dt * ins->nu);
  solver_t *velsolver   = ellipticSolveSetupTri2D(mesh, ins->lamda, kernelInfoVel, velSolverOptions); 
  ins->velsolver        = velsolver;  
  ins->velsolverOptions = velSolverOptions;

 // mesh->Nfields = ins->NTfields; 
  #endif

  kernelInfo.addDefine("p_maxNodesVolume", mymax(mesh->cubNp,mesh->Np));
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);
  #if 1
  int NblockV = 256/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 128/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  printf("maxNodes: %d \t NblockV: %d \t NblockS: %d  \n", maxNodes, NblockV, NblockS);
  
  

  ins->PrSolverID  = 0; 
  ins->PrISolverID = 1; 
  #endif
  
  printf("Np: %d \t Ncub: %d \n", mesh->Np, mesh->cubNp);



  // ADD-DEFINES 
  kernelInfo.addDefine("p_Lambda2", 0.5f);
  kernelInfo.addDefine("p_NTfields", ins->NTfields);
  kernelInfo.addDefine("p_NVfields", ins->NVfields);
  kernelInfo.addDefine("p_NfacesNfp",  mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_inu",      (float) 1.f/ins->nu);
  kernelInfo.addDefine("p_idt",      (float) 1.f/ins->dt);

  // MEMORY ALLOCATION
  ins->o_Ux  =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*ins->Nfields*sizeof(dfloat), ins->Ux);
  ins->o_Uy  =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*ins->Nfields*sizeof(dfloat), ins->Uy);  
  ins->o_Pr =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat),  ins->Pr);  
  ins->o_PrI =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->PrI);      


  ins->o_rhsUx  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsUx);
  ins->o_rhsUy  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsUy);
  ins->o_rhsPr  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsPr); 
  ins->o_NUx    = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->NUx); 
  ins->o_NUy    = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->NUy);  
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
   // KERNEL DEFINITIONS
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


  




