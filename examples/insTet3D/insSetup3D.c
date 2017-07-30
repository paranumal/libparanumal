#include "ins3D.h"

// NBN: toggle use of 2nd stream
#define USE_2_STREAMS

ins_t *insSetup3D(mesh3D *mesh,char * options, char *vSolverOptions, char *pSolverOptions, char *boundaryHeaderFileName){

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

  ins->NVfields = 3; // Total Number of Velocity Fields
  ins->NTfields = 4; // Total Velocity + Pressure
  ins->Nfields  = 1; // Each Velocity Field
  ins->ExplicitOrder = 3; // Order Nonlinear Extrapolation

  mesh->Nfields = ins->Nfields; 

  ins->mesh = mesh;
  int Nstages = 4;
  // compute samples at interpolation nodes
  ins->U     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->V     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->W     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->P     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));

  //rhs storage *3 is for velocity projection, need to be checked
  ins->rhsU  = (dfloat*) calloc(3*mesh->Nelements*mesh->Np,sizeof(dfloat));
  ins->rhsV  = (dfloat*) calloc(3*mesh->Nelements*mesh->Np,sizeof(dfloat));
  ins->rhsW  = (dfloat*) calloc(3*mesh->Nelements*mesh->Np,sizeof(dfloat));
  ins->rhsP  = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));

  //additional field storage (could reduce in the future)
  ins->NU     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->NV     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->NW     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->Px     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->Py     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->Pz     = (dfloat*) calloc(Nstages*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  ins->PI     = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));


//   if(strstr(options,"SUBCYCLING")){

//     ins->Nsubsteps = 8; //was 3

//     ins->Ud   = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
//     ins->Vd   = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
//     ins->Ue   = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
//     ins->Ve   = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
//     ins->resU = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
//     ins->resV = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));

//   }
 
  // SET SOLVER OPTIONS
  printf("Starting initial conditions for INS3D\n");
  //
 
  dfloat nu   = 0.01;   // kinematic viscosity,
  dfloat rho  = 1.0  ;  // Give density for getting actual pressure in nondimensional solve

  dfloat g[3]; g[0] = 0.0; g[1] = 0.0; g[2] = 0.0;   // No gravitational acceleration

  // Fill up required fileds
  ins->finalTime = 1.0;
  ins->nu        = nu ;
  ins->rho       = rho;
  ins->tau       = 2.f*(mesh->N+1)*(mesh->N+1); // (mesh->N)*(mesh->N+1);

  // Define total DOF per field for INS i.e. (Nelm + Nelm_halo)*Np
  ins->NtotalDofs = (mesh->totalHaloPairs+mesh->Nelements)*mesh->Np ;
  ins->NDofs      = mesh->Nelements*mesh->Np;
  // Initialize
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      const iint id = n + mesh->Np*e;
      dfloat u= 0.f, v= 0.f, w=0.f, p =0.f, t = 0.f;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];
      dfloat z = mesh->z[id];

     #if 1 //Beltrami Flow
      dfloat a = M_PI/4.0f, d = M_PI/2.0f; 
      u = -a*( exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y) )* exp(-d*d*t);
      v = -a*( exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z) )* exp(-d*d*t);
      w = -a*( exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x) )* exp(-d*d*t);

      p = -a*a*exp(-2.f*d*d*t)*( exp(2.f*a*x) +exp(2.f*a*y)+exp(2.f*a*z))*( 
          sin(a*x+d*y)*cos(a*z+d*x)*exp(a*(y+z))+
          sin(a*y+d*z)*cos(a*x+d*y)*exp(a*(x+z))+
          sin(a*z+d*x)*cos(a*y+d*z)*exp(a*(x+y))   );   
     #endif

      ins->U[id] = u;
      ins->V[id] = v;
      ins->W[id] = w;
      ins->P[id] = p;
    }
  }

  dfloat hmin = 1e9, hmax = 0;
  for(iint e=0;e<mesh->Nelements;++e){ 
     for(iint f=0;f<mesh->Nfaces;++f){
       iint sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
       dfloat sJ   = mesh->sgeo[sid + SJID];
       dfloat invJ = mesh->sgeo[sid + IJID];

       dfloat hest = 2.0/(sJ*invJ); 
       hmin = mymin(hmin, hest);
       hmax = mymax(hmax, hest);
     }
  }

  // Find Maximum Velocity
  dfloat umax = 0.f;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      const iint id = n + mesh->Np*e;
      dfloat u = ins->U[id];
      dfloat v = ins->V[id];
      dfloat w = ins->W[id];

      //Squared maximum velocity
      dfloat unmax = u*u + v*v + w*w;
      umax = mymax(umax, unmax);
    }
  }
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

  ins->NtimeSteps = ins->finalTime/ins->dt;
  ins->dt   = ins->finalTime/ins->NtimeSteps;

  // errorStep
  ins->errorStep =50;

  printf("Nsteps = %d NerrStep= %d dt = %.8e\n", ins->NtimeSteps,ins->errorStep, ins->dt);

  occa::kernelInfo kernelInfo;
  meshOccaSetup3D(mesh, deviceConfig, kernelInfo);

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
  //ins->lambda = (11.f/6.f) / (ins->dt * ins->nu);
  // ins->lambda = (1.5f) / (ins->dt * ins->nu);
  // boundaryHeaderFileName = strdup(DHOLMES "/examples/insTet3D/insVelocityEllipticBC3D.h");
  // kernelInfoV.addInclude(boundaryHeaderFileName);
  // solver_t *vSolver   = ellipticSolveSetupTet3D(mesh, ins->tau, ins->lambda, vBCType, kernelInfoV, vSolverOptions);
  // ins->vSolver        = vSolver;
  // ins->vSolverOptions = vSolverOptions;

  printf("==================PRESSURE SOLVE SETUP========================\n");
  //SETUP PRESSURE and VELOCITY SOLVERS
  // boundaryHeaderFileName = strdup(DHOLMES "/examples/insTet3D/insPressureEllipticBC3D.h");
  // kernelInfoP.addInclude(boundaryHeaderFileName);
  // solver_t *pSolver   = ellipticSolveSetupTet3D(mesh, ins->tau, 0.0, pBCType,kernelInfoP, pSolverOptions);
  // ins->pSolver        = pSolver;
  // ins->pSolverOptions = pSolverOptions;



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
  ins->o_W = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->W);
  ins->o_P = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->P);

  ins->o_rhsU  = mesh->device.malloc(3*mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsU);
  ins->o_rhsV  = mesh->device.malloc(3*mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsV);
  ins->o_rhsW  = mesh->device.malloc(3*mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsW);
  ins->o_rhsP  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsP);

  ins->o_NU = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->NU);
  ins->o_NV = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->NV);
  ins->o_NW = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->NW);
  
  ins->o_Px = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->Px);
  ins->o_Py = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->Py);
  ins->o_Pz = mesh->device.malloc(Nstages*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->Pz);
  
  ins->o_PI = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->PI);
  
  ins->o_PIx = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat));
  ins->o_PIy = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat));
  ins->o_PIz = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat));

  //storage for helmholtz solves. Fix this later !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ins->o_UH = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat));
  ins->o_VH = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat));
  ins->o_WH = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat));

  if(mesh->totalHaloPairs){
    ins->o_tHaloBuffer = mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*(ins->NTfields)*sizeof(dfloat));
    ins->o_vHaloBuffer = mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*(ins->NVfields)*sizeof(dfloat));
    ins->o_pHaloBuffer = mesh->device.malloc(mesh->totalHaloPairs*mesh->Np *sizeof(dfloat));
  }

 ins->mesh = mesh;

  // ===========================================================================
  printf("Compiling Advection volume kernel with cubature integration\n");
  ins->advectionCubatureVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection3D.okl",
				       "insAdvectionCubatureVolume3D",
				       kernelInfo);

  printf("Compiling Advection surface kernel with cubature integration\n");
  ins->advectionCubatureSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection3D.okl",
				       "insAdvectionCubatureSurface3D",
				       kernelInfo);

  printf("Compiling Advection volume kernel with collocation integration\n");
  ins->advectionVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection3D.okl",
				       "insAdvectionVolume3D",
				       kernelInfo);

  printf("Compiling Advection surface kernel with collocation integration\n");
  ins->advectionSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection3D.okl",
				       "insAdvectionSurface3D",
				       kernelInfo);

//    printf("Compiling SubCycle Advection RK update kernel\n");
//     ins->subCycleRKUpdateKernel =
//       mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",
//            "insSubCycleRKUpdate2D",
//            kernelInfo);
  // ===========================================================================
  printf("Compiling Gradient volume kernel with collocation integration\n");
  ins->gradientVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insGradient3D.okl",
				       "insGradientVolume3D",
				       kernelInfo);

  printf("Compiling Gradient volume kernel with collocation integration\n");
  ins->gradientSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insGradient3D.okl",
				       "insGradientSurface3D",
				       kernelInfo);

  // ===========================================================================
  printf("Compiling Divergence volume kernel with collocation integration\n");
  ins->divergenceVolumeKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insDivergence3D.okl",
				       "insDivergenceVolume3D",
				       kernelInfo);

  printf("Compiling Divergencesurface kernel with collocation integration\n");
  ins->divergenceSurfaceKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insDivergence3D.okl",
				       "insDivergenceSurface3D",
				       kernelInfo);

  // ===========================================================================
  printf("Compiling Helmholtz volume update kernel\n");
  ins->helmholtzRhsForcingKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHelmholtzRhs3D.okl",
				       "insHelmholtzRhsForcing3D",
				       kernelInfo);

  printf("Compiling Helmholtz IPDG RHS kernel with collocation integration\n");
  ins->helmholtzRhsIpdgBCKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHelmholtzRhs3D.okl",
				       "insHelmholtzRhsIpdgBC3D",
				       kernelInfo);


  printf("Compiling INS Helmholtz Halo Extract Kernel\n");
  ins->totalHaloExtractKernel=
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
				       "insTotalHaloExtract3D",
				       kernelInfo);

  printf("Compiling INS Helmholtz Halo Extract Kernel\n");
  ins->totalHaloScatterKernel=
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
				       "insTotalHaloScatter3D",
				       kernelInfo);


  // ===========================================================================
  printf("Compiling Helmholtz volume update kernel\n");
  ins->poissonRhsForcingKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonRhs3D.okl",
				       "insPoissonRhsForcing3D",
				       kernelInfo);

  printf("Compiling Poisson IPDG surface kernel with collocation integration\n");
  ins->poissonRhsIpdgBCKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonRhs3D.okl",
				       "insPoissonRhsIpdgBC3D",
				       kernelInfo);

  printf("Compiling Poisson penalty surface kernel\n");
  ins->poissonPenaltyKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insPoissonPenalty3D.okl",
				       "insPoissonPenalty3D",
				       kernelInfo);

  printf("Compiling INS Poisson Halo Extract Kernel\n");
  ins->velocityHaloExtractKernel=
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
				       "insVelocityHaloExtract3D",
				       kernelInfo);

  printf("Compiling INS PoissonHalo Extract Kernel\n");
  ins->velocityHaloScatterKernel=
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
				       "insVelocityHaloScatter3D",
				       kernelInfo);

  printf("Compiling Poisson surface kernel with collocation integration\n");
  ins->updateUpdateKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insUpdate3D.okl",
				       "insUpdateUpdate3D",
				       kernelInfo);


  printf("Compiling INS Poisson Halo Extract Kernel\n");
  ins->pressureHaloExtractKernel=
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
				       "insPressureHaloExtract",
				       kernelInfo);

  printf("Compiling INS PoissonHalo Extract Kernel\n");
  ins->pressureHaloScatterKernel=
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
				       "insPressureHaloScatter",
				       kernelInfo);
  // ===========================================================================//

  return ins;
}







