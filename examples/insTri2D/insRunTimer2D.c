#include "ins2D.h"

// 1 Advection Volume 2 Advection Surface 3 Ax 4 Gradient
#define KERNEL_TEST 1


void insRunTimer2D(ins_t *ins, char *options){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D *mesh = ins->mesh;
  solver_t *solver = ins->pSolver;

  ins->index = 0; 
  dfloat lambda = 0.0; 
  iint Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;

  dfloat *Z  = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  dfloat *G  = (dfloat*) calloc(Ntotal*4, sizeof(dfloat));
  
  occa::memory o_U, o_V, o_X, o_Y, o_Ud, o_Vd, o_G;
  o_U   = mesh->device.malloc(Ntotal*sizeof(dfloat),Z);
  o_V   = mesh->device.malloc(Ntotal*sizeof(dfloat),Z);
  o_Ud  = mesh->device.malloc(Ntotal*sizeof(dfloat),Z);
  o_Vd  = mesh->device.malloc(Ntotal*sizeof(dfloat),Z);
  o_X   = mesh->device.malloc(Ntotal*sizeof(dfloat),Z);
  o_Y   = mesh->device.malloc(Ntotal*sizeof(dfloat),Z);
  o_G   = mesh->device.malloc(Ntotal*4*sizeof(dfloat),G); 
  
  free(Z); free(G);

  dfloat time = 0.0; 
  iint iterations = 1000;
  // sync processes
  mesh->device.finish();
  MPI_Barrier(MPI_COMM_WORLD);

 //occaTimerTic(mesh->device,"KernelTime");
  double tic = MPI_Wtime();  
  occa::streamTag start = mesh->device.tagStream();

  // assume 1 mpi process
  for(int it=0;it<iterations;++it){
  
  #if KERNEL_TEST==1
    ins->subCycleCubatureVolumeKernel(mesh->Nelements,
          mesh->o_vgeo,
          mesh->o_cubDrWT,
          mesh->o_cubDsWT,
          mesh->o_cubInterpT,
          o_U,
          o_V,
          o_Ud,
          o_Vd,
          o_X,
          o_Y);
  #endif

  #if KERNEL_TEST==2
  ins->subCycleCubatureSurfaceKernel(mesh->Nelements,
              mesh->o_sgeo,
              mesh->o_intInterpT,
              mesh->o_intLIFTT,
              mesh->o_vmapM,
              mesh->o_vmapP,
              mesh->o_EToB,
              time,
              mesh->o_intx,
              mesh->o_inty,
              o_U,
              o_V,
              o_Ud,
              o_Vd,
              o_X,
              o_Y);
  #endif

  #if KERNEL_TEST==3
  solver->partialIpdgKernel(mesh->NinternalElements,
                                  mesh->o_internalElementIds,
                                  mesh->o_vmapM,
                                  mesh->o_vmapP,
                                  lambda,
                                  solver->tau,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  solver->o_EToB,
                                  mesh->o_DrT,
                                  mesh->o_DsT,
                                  mesh->o_LIFTT,
                                  mesh->o_MM,
                                  o_G,
                                  o_X);
  #endif

  #if KERNEL_TEST==4
  solver->partialGradientKernel(mesh->Nelements,
                                  ins->index, //zero
                                  mesh->o_vgeo,
                                  mesh->o_DrT,
                                  mesh->o_DsT,
                                  o_U,
                                  o_G);
  #endif


}
  occa::streamTag end = mesh->device.tagStream();
  mesh->device.finish();  
  double toc = MPI_Wtime();
  double localElapsed    = toc-tic;
  
  // Compute memory copy of the kernel
  iint Nbytes;
  #if KERNEL_TEST==1
  Nbytes =(sizeof(dfloat)*(6*mesh->Np +4*mesh->Np)/2);
  #endif

  #if KERNEL_TEST==2
  Nbytes = (sizeof(dfloat)*(8*mesh->Nfaces*mesh->Nfp 
                           +4*mesh->Np                   //Field variables
                           +4*mesh->Nfaces // Geometric factors - BCType
                           ))/2;
  #endif

  #if KERNEL_TEST==3
  Nbytes = (sizeof(dfloat)*(4*mesh->Np    // float4
                            +8*mesh->Nfaces*mesh->Nfp // float4
                            +5*mesh->Nfaces // cached -> *mesh->Nfp
                            +6*mesh->Np))/2;
  #endif

  #if KERNEL_TEST==4
  Nbytes = (sizeof(dfloat)*(5*mesh->Np    // float4
                            +4*mesh->Np))/2;
  #endif

  occa::memory o_foo = mesh->device.malloc(Nbytes*mesh->Nelements);
  occa::memory o_bah = mesh->device.malloc(Nbytes*mesh->Nelements);
  
  mesh->device.finish(); 
  tic = MPI_Wtime();

  occa::streamTag startCopy = mesh->device.tagStream();
  for(int it=0;it<iterations;++it){
    o_bah.copyTo(o_foo);
  }
  occa::streamTag endCopy = mesh->device.tagStream();
  
  mesh->device.finish();
  toc = MPI_Wtime();
  double copyElapsed = (toc-tic);
  copyElapsed = mesh->device.timeBetween(startCopy, endCopy);
  
  
  double copyBandwidth = mesh->Nelements*((Nbytes*iterations*2)/(1024.*1024.*1024.*copyElapsed));

  iint   localDofs = mesh->Np*mesh->Nelements;
  iint   localElements = mesh->Nelements;
  double globalElapsed;
  iint   globalDofs;
  iint   globalElements;
  int    root = 0;
  
  MPI_Reduce(&localElapsed, &globalElapsed, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD );
  MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_IINT,   MPI_SUM, root, MPI_COMM_WORLD );
  MPI_Reduce(&localElements,&globalElements,1, MPI_IINT,   MPI_SUM, root, MPI_COMM_WORLD );
  
  iint Np      = mesh->Np;
  iint Nc      = mesh->cubNp;
  iint Nfp     = mesh->Nfaces*mesh->Nfp; 
  iint Ntfp    = mesh->Nfaces*mesh->Nfp; 
  iint intNtfp = mesh->Nfaces*mesh->intNfp;
  
  double flops;
  double bw;
  

  #if KERNEL_TEST==1
    flops = Nc*Np*8 + 4*Nc + Np*Nc*16 + Np*14 + Np*2;  // All float ops only
    bw = 2.0*Nbytes;
  #endif

  #if KERNEL_TEST==2
   flops = intNtfp*Nfp*16 + intNtfp*6 + intNtfp*1 + intNtfp*28 + Np*intNtfp*4;  // All float ops only      
   bw = 2.0*Nbytes;
  #endif
  
   #if KERNEL_TEST==3
   flops = 1*Np + 21*Ntfp + Np*(4*Ntfp + 8) +6*Ntfp+ Np*(4*Np + 2)+ Np*(1+2*Ntfp) + Np*(1+2*Np);  // All float ops only      
   bw = 2.0*Nbytes;
  #endif

  #if KERNEL_TEST==4
   flops = Np*(Np*2 + Np*2 + 6);  // All float ops only      
   bw = 2.0*Nbytes;
  #endif
  
  double gflops   = globalElements*flops*iterations/(1024*1024*1024.*globalElapsed);
  double d2dbound = copyBandwidth*flops/bw;  
  bw *= (globalElements*iterations)/(1024.*1024*1024*globalElapsed);

  if(rank==root){
    printf("[ RANKS N DOFS ELAPSEDTIME ITERATIONS (DOFS/RANKS) (DOFS/TIME/ITERATIONS/RANKS) (Kernel GFLOPS) (d2dbound GB/s) (copy GB/s) (achieved GP/s)]\n");
    printf("%02d %02d %02d %17.15lg %d %17.15E %17.15E %17.15E %17.15E %17.15E %17.15E\t\n",
           size, mesh->N, globalDofs, globalElapsed, iterations, globalDofs/(double)size,
           (globalDofs*iterations)/(globalElapsed*size), gflops, d2dbound, copyBandwidth, bw);

     printf("%02d\t%02d\t%17.15lg\t%17.15E\t%17.15E\n",
              mesh->N, globalDofs, globalElapsed, gflops, d2dbound);


  }
  

  
  



  }



