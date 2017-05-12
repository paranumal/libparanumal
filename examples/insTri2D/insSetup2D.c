#include "ins2D.h"

// NBN: toggle use of 2nd stream
#define USE_2_STREAMS

  solver_t *insSetup2D(mesh2D *mesh, char * options){
  
  solver_t *ins = (solver_t*) calloc(1, sizeof(solver_t));
  ins->mesh = mesh;


  ins->Nfields = 2; // Hold Each   
  // compute samples of q at interpolation nodes
  ins->U     = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->Pr    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,sizeof(dfloat));
  
  //
  ins->rhsU  = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->rhsPr = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));
  //
  ins->UO    = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));
  ins->UI    = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));
  //
  ins->NU    = (dfloat*) calloc(mesh->Nelements*mesh->Np*ins->Nfields,sizeof(dfloat));

  
  // SET SOLVER OPTIONS 
  // Initial Conditions, Flow Properties 
  printf("Starting initial conditions for INS2D\n");
  // 
  dfloat ux   = 0.0  ; 
  dfloat uy   = 0.0  ; 
  dfloat pr   = 0.0  ;   
  dfloat nu   = 1.0/40.0;  // kinematic viscosity,
  dfloat rho  = 1.0  ;  // Give density for getting actual pressure in nondimensional solve
  
  dfloat g[2]; g[0] = 0.0; g[1] = 1.0; 
  
  
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
  //
   //memcpy(ins->g, g,2*sizeof(dfloat));



  // Define total DOF per field for INS i.e. Nelelemts*Np + Nelements_halo
  ins->NtotalDofs = mesh->Nelements*mesh->Np + mesh->totalHaloPairs ; 
  // Initialize
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      const iint id = n + mesh->Np*e;
      dfloat t = 0;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];

      dfloat lamda = 1./(2. * ins->nu) - sqrt(1./(4.*ins->nu * ins->nu) + 4.*M_PI*M_PI) ;  
      ins->U[id + UXID*ins->NtotalDofs] = 1.0 - exp(lamda*x)*cos(2.*M_PI*y); 
      ins->U[id + UYID*ins->NtotalDofs] = lamda/(2.*M_PI)*exp(lamda*x)*sin(2.*M_PI*y);
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
      dfloat uxn = ins->U[id + UXID*ins->NtotalDofs];
      dfloat uyn = ins->U[id + UYID*ins->NtotalDofs];

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
  ins->errorStep = 1000;

  printf("Nsteps = %d NerrStep= %d dt = %.8e\n", ins->NtimeSteps,ins->errorStep, ins->dt);


    // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // use rank to choose DEVICE
 // sprintf(deviceConfig, "mode = CUDA, deviceID = %d", (rank+1)%3);
  sprintf(deviceConfig, "mode = OpenCL, deviceID = 1, platformID = 0");
  // sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 0");

  // sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");  


  occa::kernelInfo kernelInfo;

  meshOccaSetup2D(mesh, deviceConfig,  kernelInfo);

    // specialization for Boltzmann

  kernelInfo.addDefine("p_maxNodesVolume", mymax(mesh->cubNp,mesh->Np));
    
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 128/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 128/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  printf("Np: %d \t Ncub: %d \n", mesh->Np, mesh->cubNp);

  printf("maxNodes: %d \t NblockV: %d \t NblockS: %d  \n", maxNodes, NblockV, NblockS);

  // ADD-DEFINES 
  kernelInfo.addDefine("p_Lambda2", 0.5f);
  kernelInfo.addDefine("p_NtotalDofs", ins->NtotalDofs);

  // kernelInfo.addDefine("p_sqrtRT", mesh->sqrtRT);
  // kernelInfo.addDefine("p_sqrt2", (float)sqrt(2.));
  // kernelInfo.addDefine("p_isq12", (float)sqrt(1./12.));
  // kernelInfo.addDefine("p_isq6", (float)sqrt(1./6.));
  // kernelInfo.addDefine("p_invsqrt2", (float)sqrt(1./2.));
  // kernelInfo.addDefine("p_tauInv", mesh->tauInv);


  // MEMORY ALLOCATION
  ins->o_U  =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*ins->Nfields*sizeof(dfloat), ins->U);

  ins->o_rhsU  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->rhsU);


  ins->o_Pr =    
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), ins->Pr);

  ins->o_NU  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->NU);  
   
  ins->o_UO  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->UO);

  ins->o_UI  = mesh->device.malloc(mesh->Np*mesh->Nelements*ins->Nfields*sizeof(dfloat), ins->UI);



   // KERNEL DEFINITIONS
  // if(strstr(options, "CUBATURE")){ 
  // printf("Compiling Advection volume kernel with cubature integration\n");
  // mesh->volumeKernel =
  //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection2D.okl",
  //     "insAdvectionVolumeCub2D",
  //       kernelInfo);
  // }
  // else{
  printf("Compiling Advection volume kernel with collocation integration\n");
  ins->advectionVolumeKernel = 
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection2D.okl",
      "insAdvectionVolume2D",
        kernelInfo);
  // }

   // KERNEL DEFINITIONS
  // if(strstr(options, "CUBATURE")){ 
  // printf("Compiling Advection volume kernel with cubature integration\n");
  // mesh->volumeKernel =
  //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection2D.okl",
  //     "insAdvectionVolumeCub2D",
  //       kernelInfo);
  // }
  // else{
  printf("Compiling Advection volume kernel with collocation integration\n");
  ins->advectionSurfaceKernel = 
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection2D.okl",
      "insAdvectionSurface2D",
        kernelInfo);
  // }


  printf("Compiling Advection volume kernel with collocation integration\n");
  ins->advectionUpdateKernel = 
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insAdvection2D.okl",
      "insAdvectionUpdate2D",
        kernelInfo);

  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
      "meshHaloExtract2D",
      kernelInfo);

  printf("Compiling INS Halo Extract Kernel\n");
  ins->haloExtractKernel= 
    mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExtract.okl",
      "insHaloExtract2D",
        kernelInfo);

return ins;

}


  




