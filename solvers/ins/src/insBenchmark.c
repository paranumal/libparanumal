#include "ins.h"
// 1 Advection Volume 2 Advection Surface 3 Ax  4 Gradient
// for optimization
// 5 advaction volume 6 advection surface 7 Ax kernel 8 Gradient
// #define KERNEL_TEST 7

void insBenchmark(mesh_t *mesh, setupAide options){

  ins_t *ins   = (ins_t*) calloc(1, sizeof(ins_t));
  ins->mesh    = mesh; 
  ins->options = options;
  //
  options.getArgs("MESH DIMENSION", ins->dim);
  options.getArgs("ELEMENT TYPE", ins->elementType);
  //
  ins->NVfields = (ins->dim==3) ? 3:2; //  Total Number of Velocity Fields
  ins->NTfields = (ins->dim==3) ? 4:3; // Total Velocity + Pressure

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh->Nfields = 1;

  dlong Nlocal = mesh->Np*mesh->Nelements;
  dlong Ntotal = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

  ins->Ntotal = Ntotal;
  ins->fieldOffset = Ntotal;
  ins->Nblock = (Nlocal+blockSize-1)/blockSize;


  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();


  if(ins->dim==3)
    meshOccaSetup3D(mesh, options, kernelInfo);
  else
    meshOccaSetup2D(mesh, options, kernelInfo);

  occa::properties kernelInfoV  = kernelInfo;
  occa::properties kernelInfoP  = kernelInfo;

  ins->Nstages = 1; 

  // CPU allocation
  ins->U     = (dfloat*) calloc(ins->NVfields*ins->Nstages*Ntotal,sizeof(dfloat));
  ins->P     = (dfloat*) calloc(              ins->Nstages*Ntotal,sizeof(dfloat));
  ins->GP   =  (dfloat*) calloc(             4*ins->Nstages*Ntotal,sizeof(dfloat));


  //rhs storage
  ins->rhsU  = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  ins->rhsV  = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  ins->rhsW  = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  ins->rhsP  = (dfloat*) calloc(Ntotal,sizeof(dfloat));

  //additional field storage
  ins->NU   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));
  ins->LU   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));
  // ins->GP   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));

  //
  ins->GU   = (dfloat*) calloc(ins->NVfields*Ntotal*4,sizeof(dfloat));
  
  ins->rkU  = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkP  = (dfloat*) calloc(              Ntotal,sizeof(dfloat));
  ins->PI   = (dfloat*) calloc(              Ntotal,sizeof(dfloat));
  
  ins->rkNU = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkLU = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkGP = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));


   // MEMORY ALLOCATION
  ins->o_U = mesh->device.malloc(ins->NVfields*ins->Nstages*Ntotal*sizeof(dfloat), ins->U);
  ins->o_P = mesh->device.malloc(              ins->Nstages*Ntotal*sizeof(dfloat), ins->P);

  ins->o_rhsU  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsU);
  ins->o_rhsV  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsV);
  ins->o_rhsW  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsW);
  ins->o_rhsP  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsP);

  ins->o_NU    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->NU);
  ins->o_LU    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->LU);
  // ins->o_GP    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->GP);
  ins->o_GP    = mesh->device.malloc(             4*ins->Nstages*Ntotal*sizeof(dfloat), ins->GP);
  
  ins->o_GU    = mesh->device.malloc(ins->NVfields*Ntotal*4*sizeof(dfloat), ins->GU);
  
  ins->o_rkU   = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkU);
  ins->o_rkP   = mesh->device.malloc(              Ntotal*sizeof(dfloat), ins->rkP);
  ins->o_PI    = mesh->device.malloc(              Ntotal*sizeof(dfloat), ins->PI);
  
  ins->o_rkNU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkNU);
  ins->o_rkLU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkLU);
  ins->o_rkGP  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkGP);


  ins->Ud    = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->Ue    = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->resU  = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rhsUd = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));

  if(ins->elementType==HEXAHEDRA){
   ins->cUd = (dfloat *) calloc(ins->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
   ins->cU = (dfloat *) calloc(ins->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
  }
  else{ 
   ins->cUd = ins->U;
   ins->cU  = ins->U;
 }


 
  kernelInfo["defines/" "p_blockSize"]= blockSize;
  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = mymax(1,256/mesh->Np); // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = mymax(1,256/maxNodes); // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;


  int maxNodesVolumeCub = mymax(mesh->cubNp,mesh->Np);  
  kernelInfo["defines/" "p_maxNodesVolumeCub"]= maxNodesVolumeCub;
  int cubNblockV = mymax(1,256/maxNodesVolumeCub);
  //
  int maxNodesSurfaceCub = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  kernelInfo["defines/" "p_maxNodesSurfaceCub"]=maxNodesSurfaceCub;
  int cubNblockS = mymax(256/maxNodesSurfaceCub,1); // works for CUDA
  //
  kernelInfo["defines/" "p_cubNblockV"]=cubNblockV;
  kernelInfo["defines/" "p_cubNblockS"]=cubNblockS;

  //
  kernelInfo["defines/" "p_NTfields"]= ins->NTfields;
  kernelInfo["defines/" "p_NVfields"]= ins->NVfields;
  kernelInfo["defines/" "p_NfacesNfp"]=  mesh->Nfaces*mesh->Nfp;
  kernelInfo["defines/" "p_Nstages"]=  ins->Nstages;
  kernelInfo["defines/" "p_SUBCYCLING"]=  ins->Nsubsteps;

  // 
  ins->o_U = mesh->device.malloc(ins->NVfields*ins->Nstages*Ntotal*sizeof(dfloat), ins->U);
  ins->o_P = mesh->device.malloc(              ins->Nstages*Ntotal*sizeof(dfloat), ins->P);

  // Note that resU and resV can be replaced with already introduced buffer
  ins->o_Ue    = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Ue);
  ins->o_Ud    = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Ud);
  ins->o_resU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->resU);
  ins->o_rhsUd = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rhsUd);

  if(ins->elementType==HEXAHEDRA){
    ins->o_cU = mesh->device.malloc(ins->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->cU);
    ins->o_cUd = mesh->device.malloc(ins->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->cUd);
  }
  else{ 
    ins->o_cUd = ins->o_Ud;
    ins->o_cU  = ins->o_U;
  }


  //add boundary data to kernel info
  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo["includes"] += (char*)boundaryHeaderFileName.c_str();

  int iterations =1000;
  dfloat tic =0.f, toc = 0.f, kernelElapsed=0.f; 
  occa::kernel TestKernel; 
  

  int NKernels;
  options.getArgs("BENCHMARK KERNEL NUMBER", NKernels);
  occa::kernel *testKernels = new occa::kernel[NKernels];
  char kernelName[NKernels][BUFSIZ];

  // Define Variables for timing and roofline models
  dfloat *kernelShared  = (dfloat *)calloc(NKernels, sizeof(dfloat));
  dfloat *kernelBytes   = (dfloat *)calloc(NKernels, sizeof(dfloat));
  // dfloat *kernelSharedB = (dfloat *)calloc(NKernels, sizeof(dfloat));
  dfloat *kernelFlops   = (dfloat *)calloc(NKernels, sizeof(dfloat));
  dfloat *kernelTime    = (dfloat *)calloc(NKernels, sizeof(dfloat));
  dfloat *kernelCopyTime= (dfloat *)calloc(NKernels, sizeof(dfloat));


  // set kernel name suffix
  char *suffix;
  char fileName[BUFSIZ]; 
  
  if(ins->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(ins->elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(ins->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(ins->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  // Find kernels and build them
  if(options.compareArgs("BENCHMARK KERNEL", "CUBATURE_VOLUME")){
    printf("Benchmarking cubature volume kernel....\n");
    for(int i=0; i<NKernels; i++){
      sprintf(fileName, DINS "/okl/insOptAdvectionVolume%s.okl", suffix);
      sprintf(kernelName[i], "insSubCycleCubatureVolume%s_v%d", suffix, i);
      
      printf("Building %s kernel...", kernelName[i]);
      testKernels[i] = mesh->device.buildKernel(fileName, kernelName[i], kernelInfo);
      printf("done\n");
    }
  }


  // Find kernels and build them
  if(options.compareArgs("BENCHMARK KERNEL", "CUBATURE_SURFACE")){
    printf("Benchmarking cubature surface kernel....\n");
    for(int i=0; i<NKernels; i++){
      sprintf(fileName, DINS "/okl/insOptAdvectionSurface%s.okl", suffix);
      sprintf(kernelName[i], "insSubCycleCubatureSurface%s_v%d", suffix, i);
      
      printf("Building %s kernel...", kernelName[i]);
      testKernels[i] = mesh->device.buildKernel(fileName, kernelName[i], kernelInfo);
      printf("done\n");
    }
  }

  // Find kernels and build them
  if(options.compareArgs("BENCHMARK KERNEL", "GRADIENT")){
    printf("Benchmarking gradient kernel....\n");
    for(int i=0; i<NKernels; i++){
      sprintf(fileName, DINS "/okl/insOptGradient%s.okl", suffix);
      sprintf(kernelName[i], "ellipticPartialGradient%s_v%d", suffix, i);
      
      printf("Building %s kernel...", kernelName[i]);
      testKernels[i] = mesh->device.buildKernel(fileName, kernelName[i], kernelInfo);
      printf("done\n");
    }
  }


  // Find kernels and build them
  if(options.compareArgs("BENCHMARK KERNEL", "AXIPDG")){
    printf("Benchmarking AxIPDG kernel....\n");
    for(int i=0; i<NKernels; i++){
      sprintf(fileName, DINS "/okl/insOptAxIpdg%s.okl", suffix);
      sprintf(kernelName[i], "ellipticPartialAxIpdg%s_v%d", suffix, i);
      
      printf("Building %s kernel...", kernelName[i]);
      testKernels[i] = mesh->device.buildKernel(fileName, kernelName[i], kernelInfo);
      printf("done\n");
    }
  }

  // Run kernels 
  for(int i=0; i<NKernels; i++){
    MPI_Barrier(MPI_COMM_WORLD);
    mesh->device.finish();

    // #if KERNEL_TEST==1
    if(options.compareArgs("BENCHMARK KERNEL", "CUBATURE_VOLUME")){
      //occaTimerTic(mesh->device,"KernelTime");
      tic = MPI_Wtime();  
      occa::streamTag start = mesh->device.tagStream();

      // assume 1 mpi process
      for(int it=0;it<iterations;++it){
        //printf("Cubature Points: %d", mesh->cubNp);
        testKernels[i](mesh->Nelements,
                      mesh->o_vgeo,
                      mesh->o_cubvgeo,
                      mesh->o_cubDWmatrices,
                      mesh->o_cubInterpT,
                      mesh->o_cubProjectT,
                      ins->fieldOffset,
                      ins->o_Ue,
                      ins->o_Ud,
                      ins->o_cU,     
                      ins->o_cUd,     
                      ins->o_rhsUd);
      }

      occa::streamTag end = mesh->device.tagStream();

      mesh->device.finish();  
      toc = MPI_Wtime();
      kernelTime[i]    = mesh->device.timeBetween(start,end)/iterations;

      // Needs to adjusted for other kernels!!!!!
      if(i==0){
       kernelBytes[i]  = (sizeof(dfloat)*(4*mesh->Np*mesh->cubNp*0 +4*mesh->Np + 2*mesh->Np)/2);
       kernelShared[i] = (sizeof(dfloat)*(4*mesh->cubNp + 4*mesh->Np*mesh->cubNp)); 
       kernelFlops[i]  = mesh->cubNp*mesh->Np*8 + 4*mesh->cubNp + mesh->Np*mesh->cubNp*16 + mesh->Np*14 + mesh->Np*2;
      }else{
       kernelBytes[i]  = (sizeof(dfloat)*(4*mesh->Np +4*mesh->Np + 2*mesh->Np)/2);
       kernelShared[i] = (sizeof(dfloat)*(4*mesh->Np + 4*mesh->Np*mesh->cubNp + 4*mesh->cubNp + 4*mesh->Np*mesh->cubNp)); 
       kernelFlops[i]  = mesh->Np*6 + mesh->Np*mesh->cubNp*8 + 4*mesh->cubNp + 8*mesh->Np*mesh->cubNp + 2*mesh->Np ;
      }
    }

     // #if KERNEL_TEST==1
    if(options.compareArgs("BENCHMARK KERNEL", "CUBATURE_SURFACE")){
       dfloat bScale = 0.0; 
       dfloat time   = 0.0; 
      //occaTimerTic(mesh->device,"KernelTime");
      tic = MPI_Wtime();  
      occa::streamTag start = mesh->device.tagStream();
      // assume 1 mpi process
      for(int it=0;it<iterations;++it){
        // //printf("Cubature Points: %d", mesh->cubNp);
        testKernels[i](mesh->Nelements,
                      mesh->o_vgeo,
                      mesh->o_sgeo,
                      mesh->o_cubsgeo,
                      mesh->o_intInterpT,
                      mesh->o_intLIFTT,
                      mesh->o_cubInterpT,
                      mesh->o_cubProjectT,
                      mesh->o_vmapM,
                      mesh->o_vmapP,
                      mesh->o_EToB,
                      bScale,
                      time,
                      mesh->o_intx,
                      mesh->o_inty,
                      mesh->o_intz,
                      ins->fieldOffset,
                      ins->o_Ue,
                      ins->o_Ud,
                      ins->o_rhsUd);
      }

      occa::streamTag end = mesh->device.tagStream();

      mesh->device.finish();  
      toc = MPI_Wtime();
      // kernelElapsed = toc-tic;
      kernelTime[i]    = mesh->device.timeBetween(start,end)/iterations;
      if(i==0){
       kernelBytes[i]  = (sizeof(dfloat)*(4*mesh->Np*mesh->cubNp*0 +4*mesh->Np + 2*mesh->Np)/2);
       kernelShared[i] = (sizeof(dfloat)*(4*mesh->cubNp + 4*mesh->Np*mesh->cubNp)); 
       kernelFlops[i]  = mesh->cubNp*mesh->Np*8 + 4*mesh->cubNp + mesh->Np*mesh->cubNp*16 + mesh->Np*14 + mesh->Np*2;  // All float ops only
      }else{
       kernelBytes[i]  = (sizeof(dfloat)*(4*mesh->Np +4*mesh->Np + 2*mesh->Np)/2);
       kernelShared[i] = (sizeof(dfloat)*(4*mesh->Np + 4*mesh->Np*mesh->cubNp + 4*mesh->cubNp + 4*mesh->Np*mesh->cubNp)); 
       kernelFlops[i]  = mesh->Np*6 + mesh->Np*mesh->cubNp*8 + 4*mesh->cubNp + 8*mesh->Np*mesh->cubNp + 2*mesh->Np ;  // All float ops only
      }
    }



     // #if KERNEL_TEST==1
    if(options.compareArgs("BENCHMARK KERNEL", "GRADIENT")){
      dlong offset = 0;  
      //occaTimerTic(mesh->device,"KernelTime");
      tic = MPI_Wtime();  
      occa::streamTag start = mesh->device.tagStream();
      // assume 1 mpi process
      for(int it=0;it<iterations;++it){
        // //printf("Cubature Points: %d", mesh->cubNp);
        testKernels[i](mesh->Nelements,
                      offset,
                      mesh->o_vgeo,
                      mesh->o_Dmatrices,
                      ins->o_P,
                      ins->o_GP);
      }

      occa::streamTag end = mesh->device.tagStream();

      mesh->device.finish();  
      toc = MPI_Wtime();
      // kernelElapsed = toc-tic;
      kernelTime[i]    = mesh->device.timeBetween(start,end)/iterations;
       const dfloat alpha = 1.0;
      if(i==5){
       kernelBytes[i]  = (sizeof(dfloat)*(1*mesh->Np +4*alpha*mesh->Np + 4*(1.0-alpha)*1 + 4*mesh->Np)/2);
       kernelShared[i] = (sizeof(dfloat)*(1*mesh->Np + mesh->Np*(2*mesh->Np) + 4 + 4*mesh->Np + mesh->Np));  
       kernelFlops[i]  = mesh->Np*(mesh->Np*4 +6);
      }else{
       kernelBytes[i]  = (sizeof(dfloat)*(1*mesh->Np +4*alpha*mesh->Np + 4*(1.0-alpha)*1 + 4*mesh->Np)/2);
       kernelShared[i] = (sizeof(dfloat)*(1*mesh->Np + mesh->Np*(2*mesh->Np) + 1*mesh->Np)); 
       kernelFlops[i]  = mesh->Np*(mesh->Np*4 +6);  
      }
    }


     // #if KERNEL_TEST==1
    if(options.compareArgs("BENCHMARK KERNEL", "AXIPDG")){
      dlong offset = 0;  
      dfloat lambda = 0.0;
      dfloat tau   = 10.0; 

      //occaTimerTic(mesh->device,"KernelTime");
      tic = MPI_Wtime();  
      occa::streamTag start = mesh->device.tagStream();
      // assume 1 mpi process
      for(int it=0;it<iterations;++it){
        // //printf("Cubature Points: %d", mesh->cubNp);
        testKernels[i]( mesh->NinternalElements,
                        mesh->o_internalElementIds,
                        mesh->o_vmapM,
                        mesh->o_vmapP,
                        lambda,
                        tau,
                        mesh->o_vgeo,
                        mesh->o_sgeo,
                        mesh->o_EToB,
                        mesh->o_Dmatrices,
                        mesh->o_LIFTT,
                        mesh->o_MM,
                        ins->o_GP,
                        ins->o_rhsP);
      }

      occa::streamTag end = mesh->device.tagStream();

      mesh->device.finish();  
      toc = MPI_Wtime();
      // kernelElapsed = toc-tic;
      kernelTime[i]    = mesh->device.timeBetween(start,end)/iterations;
       
       dfloat Ntfp = mesh->Nfaces*mesh->Nfp; 
       dfloat alpha = 1.0; // 0.0; 

       kernelShared[i] = (sizeof(dfloat)*(3*mesh->Np + 3*Ntfp + mesh->Np*(2*Ntfp +2+4) + Ntfp*3 + mesh->Np*(2*mesh->Np + 1) + mesh->Np*(1*Ntfp + 1) + mesh->Np*(1*mesh->Np)));   
       kernelFlops[i]  = 1*mesh->Np + 21*Ntfp + mesh->Np*(4*Ntfp + 8) + 6*Ntfp+ mesh->Np*(4*mesh->Np + 2)+ mesh->Np*(1+2*Ntfp) + mesh->Np*(1+2*mesh->Np); 
       
      if(i==3){
       kernelBytes[i]  = (sizeof(dfloat)*(4*mesh->Np    // float4
                                        +4*Ntfp  // float4
                                        +alpha*5*mesh->Nfaces + (1-alpha)*5*Ntfp   // alpha = 1 means cached face geometric factors
                                        +alpha*5 +(1-alpha)*5*mesh->Np   // alpha 1 means drdx, J etc.. are all cached for an element
                                        +1*mesh->Np))/2           // gloabl write of output
                       + (sizeof(int)*(alpha*1 + (1-alpha)*mesh->Np*1  // element list//  was 3
                                        +2*Ntfp  // idM, idP
                                        +alpha*mesh->Nfaces + (1-alpha)*Ntfp // BC's
                                        ))/2;  
      }else{
       kernelBytes[i]  = (sizeof(dfloat)*(4*mesh->Np    // float4
                                         +8*Ntfp  // float4
                                         +alpha*5*mesh->Nfaces + (1-alpha)*5*Ntfp   // alpha = 1 means cached face geometric factors
                                         +alpha*5 +(1-alpha)*5*mesh->Np   // alpha 1 means drdx, J etc.. are all cached for an element
                                         +1*mesh->Np))/2           // gloabl write of output
                      + (sizeof(int)*(alpha*1 + (1-alpha)*mesh->Np*1  // element list//  was 3
                                          +2*Ntfp  // idM, idP
                                          +alpha*mesh->Nfaces + (1-alpha)*Ntfp // BC's
                                          ))/2;  
      }
    }
    // printf("%02d %02d %02d %02d %12.10E\n",i, mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), kernelTime[i]/iterations);
  }



  // Kernel Copy Times
  for(int i=0; i<NKernels; i++){
    occa::memory o_foo = mesh->device.malloc(kernelBytes[i]*mesh->Nelements);
    occa::memory o_bah = mesh->device.malloc(kernelBytes[i]*mesh->Nelements);

    mesh->device.finish(); 
    occa::streamTag startCopy = mesh->device.tagStream();
    for(int it=0;it<iterations;++it){
      o_bah.copyTo(o_foo);
    }
    occa::streamTag endCopy = mesh->device.tagStream();
    mesh->device.finish();

    kernelCopyTime[i] = mesh->device.timeBetween(startCopy, endCopy)/iterations;
  }


  double max_smb_p100  = 7882.0;
  double max_gfl_p100  = 549.0;
  double max_thg_p100  = 4670; 

    for(int i=0; i<NKernels; i++){

    double copyBandwidth = mesh->Nelements*((kernelBytes[i]*2)/(1e9*kernelCopyTime[i]));
    double  bw           = mesh->Nelements*((kernelBytes[i]*2)/(1e9*kernelTime[i]));

    double gflops        = mesh->Nelements*kernelFlops[i]/(1e9*kernelTime[i]);
    double d2dbound      = copyBandwidth*gflops/bw;

    double smbound       = max_smb_p100*kernelFlops[i]/( (double) kernelShared[i]);
    double intensity     = gflops/bw; 
    double roofline      = mymin(d2dbound, smbound);
    double ach_thg       = mymin(max_gfl_p100*intensity, max_thg_p100);

#if 0
    printf("[ N\tK\tDOFS\tKernelTime\tCopyTime\tIntensity\tGFLOPS/s\t(d2dbound GFLOPS/s)\tBW(GB/s)\t(SMBOUND GFLOPS/s)\t(ROOFLINE GFLOPS/s)\tTH_peak]\n");
    printf("%02d \t%02d\t%02d\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\n",
        mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), kernelTime[i], kernelCopyTime[i], intensity, gflops, d2dbound, bw, smbound, roofline, ach_thg);
#endif

    char fname[BUFSIZ];
    FILE *fp; 
    if(options.compareArgs("BENCHMARK KERNEL", "AXIPDG"))
      sprintf(fname, "AxKernelData.dat");
    if(options.compareArgs("BENCHMARK KERNEL", "GRADIENT"))
      sprintf(fname, "GradientKernelData.dat");
    if(options.compareArgs("BENCHMARK KERNEL", "CUBATURE_SURFACE"))
      sprintf(fname, "CubSurfaceKernelData.dat");
    if(options.compareArgs("BENCHMARK KERNEL", "CUBATURE_VOLUME"))
      sprintf(fname, "CubVolumeKernelData.dat");

    fp = fopen(fname,"a");
    fprintf(fp, "%02d %12.10E %12.10E %12.10E %12.10E\n", mesh->N, kernelTime[i], kernelCopyTime[i], gflops, roofline);
    fclose(fp);
}
        

free(kernelTime);
free(kernelBytes);
free(kernelShared);
free(kernelFlops);
free(kernelCopyTime);

}



// /////////////////////////////////////////////////////
      // int Nbytes = 0; 
      // int NbytesShared = 0; 
      // int NbytesShared2= 0; 
      // int flops = 0; 

      // if(i==0){
      //   Nbytes       = (sizeof(dfloat)*(4*mesh->Np*mesh->cubNp*0 +4*mesh->Np + 2*mesh->Np)/2);
      //   NbytesShared = (sizeof(dfloat)*(4*mesh->cubNp + 4*mesh->Np*mesh->cubNp)); 
      //   flops        = mesh->cubNp*mesh->Np*8 + 4*mesh->cubNp + mesh->Np*mesh->cubNp*16 + mesh->Np*14 + mesh->Np*2;  // All float ops only
      // }
      //  else
      // {
          
      //   Nbytes        = (sizeof(dfloat)*(4*mesh->Np +4*mesh->Np + 2*mesh->Np)/2);
      //   NbytesShared  = (sizeof(dfloat)*(4*mesh->Np + 4*mesh->Np*mesh->cubNp + 4*mesh->cubNp + 4*mesh->Np*mesh->cubNp)); 
      //   NbytesShared2 = (sizeof(dfloat)*(4*mesh->Np + 4*mesh->Np*mesh->cubNp + 4*mesh->cubNp + 4*mesh->Np*mesh->cubNp + 1.*mesh->Np*mesh->cubNp + 2.*mesh->Np*mesh->cubNp )); // Add operators
      //   flops         = mesh->Np*6 + mesh->Np*mesh->cubNp*8 + 4*mesh->cubNp + 8*mesh->Np*mesh->cubNp + 2*mesh->Np ;  // All float ops only

      // }


        // occa::memory o_foo = mesh->device.malloc(Nbytes*mesh->Nelements);
        // occa::memory o_bah = mesh->device.malloc(Nbytes*mesh->Nelements);

        // mesh->device.finish(); 
        // tic = MPI_Wtime();

        // occa::streamTag startCopy = mesh->device.tagStream();
        // for(int it=0;it<iterations;++it){
        //    o_bah.copyTo(o_foo);
        // }
        // occa::streamTag endCopy = mesh->device.tagStream();

        // mesh->device.finish();
        // toc = MPI_Wtime();
        // //      double copyElapsed = (toc-tic);
        // double copyElapsed = mesh->device.timeBetween(startCopy, endCopy);

        //    // Compute Data
        // //double peakSharedStreaming = 
        // double copyBandwidth = mesh->Nelements*((Nbytes*iterations*2)/(1e9*copyElapsed));
        // double  bw           = mesh->Nelements*((Nbytes*iterations*2)/(1e9*kernelElapsed));

        // double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);
        // double d2dbound      = copyBandwidth*gflops/bw;

        // // TW: avoid using hard coded numbers like this in formula, define a new variable peakSharedStreaming
        // //     7882 is defined wrt to 10^9
        // double smbound       = 7882*flops/( (double) NbytesShared);
        // double intensity     = gflops/bw; 
        // // double l1bound      = 7882*flops/( (double) NbytesShared2);

        // // printf("l1bound :%.2e \n", l1bound);



        // double roofline1     = mymin(d2dbound, smbound);
        // // double roofline2     = mymin(d2dbound, l1bound);
        // double max_thg_p100  = 4670; 
        // double ach_thg       = mymin(549*intensity, max_thg_p100);



        // printf("[ N\tK\tDOFS\tKernelTime\tCopyTime\tIntensity\tGFLOPS/s\t(d2dbound GFLOPS/s)\tBW(GB/s)\t(SMBOUND GFLOPS/s)\t(ROOFLINE GFLOPS/s)\t(ROOFLINE2 GFLOPS/s)\tTH_peak]\n");
        // printf("%02d \t%02d\t%02d\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\n",
        //         mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
        //  kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline1, ach_thg);

        // char fname[BUFSIZ];
        // sprintf(fname, "KernelData.dat");
        // FILE *fp;
        // fp = fopen(fname, "a");

        // fprintf(fp, "%02d %02d %02d %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E\n",
        //         mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
        //         kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline1,  ach_thg);
        // fclose(fp);


        // sprintf(fname, "KernelTime.dat");
        // fp = fopen(fname, "a");

        // fprintf(fp, "%02d %02d %02d %12.10E %12.10E\n",
        //         mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), kernelElapsed/iterations, copyElapsed/iterations);
        // fclose(fp);


  //   // #if KERNEL_TEST==1
  // if(options.compareArgs("BENCHMARK_KERNEL", "CUBATURE_SURFACE")){

  //   int NKernels = 7;
  //   occa::kernel *testKernels = new occa::kernel[NKernels];
  //   char kernelName[NKernels][BUFSIZ];
  //   printf("Benchmarking cubature volume kernel....\n");
  //   for(int i=0; i<NKernels; i++){

  //     sprintf(fileName, DINS "/okl/insOptAdvection%s.okl", suffix);
  //     sprintf(kernelName[i], "insSubCycleCubatureVolume%s_v%d", suffix, i);

  //     printf("Building %s kernel...", kernelName[i]);
  //     testKernels[i] = mesh->device.buildKernel(fileName, kernelName[i], kernelInfo);
  //     // testKernels1 = mesh->device.buildKernel(fileName, kernelName[i], kernelInfo);
  //     printf("done\n");
      
  //     // sync processes
  //     MPI_Barrier(MPI_COMM_WORLD);
  //     mesh->device.finish();
      

  //     //occaTimerTic(mesh->device,"KernelTime");
  //     tic = MPI_Wtime();  

  //     occa::streamTag start = mesh->device.tagStream();

  //       // assume 1 mpi process
  //       for(int it=0;it<iterations;++it){
  //         //printf("Cubature Points: %d", mesh->cubNp);
  //         testKernels[i](mesh->Nelements,
  //                        mesh->o_vgeo,
  //                        mesh->o_cubvgeo,
  //                        mesh->o_cubDWmatrices,
  //                        mesh->o_cubInterpT,
  //                        mesh->o_cubProjectT,
  //                        ins->fieldOffset,
  //                        ins->o_Ue,
  //                        ins->o_Ud,
  //                        ins->o_cU,     
  //                        ins->o_cUd,     
  //                        ins->o_rhsUd);
  //       }

  //       occa::streamTag end = mesh->device.tagStream();
        
  //       mesh->device.finish();  
  //       toc = MPI_Wtime();
  //       //      kernelElapsed    = toc-tic;
  //       kernelElapsed = mesh->device.timeBetween(start,end)/iterations;

  //       int Nbytes = 0; 
  //       int NbytesShared = 0; 
  //       int NbytesShared2= 0; 
  //       int flops = 0; 

  //       int Np = mesh->Np; 
  //       int Nc = mesh->cubNp; 


  //       dfloat nmt =1; 

  //       if(i==0){
  //         Nbytes       = (sizeof(dfloat)*(4*Np*Nc*0 +4*Np + 2*Np)/2);
  //         NbytesShared = (sizeof(dfloat)*(4*Nc + 4*Np*Nc)); 
  //         flops        = Nc*Np*8 + 4*Nc + Np*Nc*16 + Np*14 + Np*2;  // All float ops only
  //       }
  //        else
  //       {
            
  //         Nbytes        = (sizeof(dfloat)*(4*Np +4*Np + 2*Np)/2);
  //         NbytesShared  = (sizeof(dfloat)*(4*Np + 4*Np*Nc + 4*Nc + 4*Np*Nc)); 
  //         NbytesShared2 = (sizeof(dfloat)*(4*Np + 4*Np*Nc + 4*Nc + 4*Np*Nc + 1.*Np*Nc/nmt + 2.*Np*Nc/nmt )); // Add operators
  //         flops         = Np*6 + Np*Nc*8 + 4*Nc + 8*Np*Nc + 2*Np ;  // All float ops only

  //       }


  //       occa::memory o_foo = mesh->device.malloc(Nbytes*mesh->Nelements);
  //       occa::memory o_bah = mesh->device.malloc(Nbytes*mesh->Nelements);

  //       mesh->device.finish(); 
  //       tic = MPI_Wtime();

  //       occa::streamTag startCopy = mesh->device.tagStream();
  //       for(int it=0;it<iterations;++it){
  //          o_bah.copyTo(o_foo);
  //       }
  //       occa::streamTag endCopy = mesh->device.tagStream();

  //       mesh->device.finish();
  //       toc = MPI_Wtime();
  //       //      double copyElapsed = (toc-tic);
  //       double copyElapsed = mesh->device.timeBetween(startCopy, endCopy);

  //          // Compute Data
  //       //double peakSharedStreaming = 
  //       double copyBandwidth = mesh->Nelements*((Nbytes*iterations*2)/(1e9*copyElapsed));
  //       double  bw           = mesh->Nelements*((Nbytes*iterations*2)/(1e9*kernelElapsed));

  //       double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);
  //       double d2dbound      = copyBandwidth*gflops/bw;

  //       // TW: avoid using hard coded numbers like this in formula, define a new variable peakSharedStreaming
  //       //     7882 is defined wrt to 10^9
  //       double smbound       = 7882*flops/( (double) NbytesShared);
  //       double intensity     = gflops/bw; 
  //       double l1bound      = 7882*flops/( (double) NbytesShared2);

  //       printf("l1bound :%.2e \n", l1bound);



  //       double roofline1     = mymin(d2dbound, smbound);
  //       double roofline2     = mymin(d2dbound, l1bound);
  //       double max_thg_p100  = 4670; 
  //       double ach_thg       = mymin(549*intensity, max_thg_p100);



  //       printf("[ N\tK\tDOFS\tKernelTime\tCopyTime\tIntensity\tGFLOPS/s\t(d2dbound GFLOPS/s)\tBW(GB/s)\t(SMBOUND GFLOPS/s)\t(ROOFLINE GFLOPS/s)\t(ROOFLINE2 GFLOPS/s)\tTH_peak]\n");
  //       printf("%02d \t%02d\t%02d\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\n",
  //               mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
  //        kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline1, roofline2, ach_thg);

  //       char fname[BUFSIZ];
  //       sprintf(fname, "KernelData.dat");
  //       FILE *fp;
  //       fp = fopen(fname, "a");

  //       fprintf(fp, "%02d %02d %02d %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E\n",
  //               mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
  //               kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline1, roofline2, ach_thg);
  //       fclose(fp);


  //       sprintf(fname, "KernelTime.dat");
  //       fp = fopen(fname, "a");

  //       fprintf(fp, "%02d %02d %02d %12.10E %12.10E\n",
  //               mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), kernelElapsed/iterations, copyElapsed/iterations);
  //       fclose(fp);


  //   }


  



  // // SURFACE KERNEL
  // #if KERNEL_TEST==2
  
  // int NKernels = 8;

  // dfloat  NMT[10] = {2,1,2,2,2,2,2,2,3,2};

  // occa::kernel *testKernels = new occa::kernel[NKernels];

  // char kernelNames[NKernels][BUFSIZ];

  // for(int i=0; i<NKernels; i++)
  // {
        
  //   sprintf(kernelNames[i], "insSubCycleCubatureSurface2D_v%d", i);

  //   testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/optSubcycleSurface2D.okl",kernelNames[i], kernelInfo);

  //   printf("insSubCycleCubatureSurface Kernel #%02d\n", i);
  //   printf("Nblock: %d cubNblock: %d N: %d Np: %d cubNp: %d\n", NblockV, cubNblockV, mesh->N, mesh->Np, mesh->cubNp);

    
  //    // sync processes
  //   mesh->device.finish();
  //   //    MPI_Barrier(MPI_COMM_WORLD);

  //   // if(i==7){

  //   // //occaTimerTic(mesh->device,"KernelTime");
  //   // tic = MPI_Wtime();  

  //   // occa::streamTag start = mesh->device.tagStream();

  //   //   // assume 1 mpi process
  //   //   for(int it=0;it<iterations;++it){
  //   //     //printf("Cubature Points: %d", mesh->cubNp);
  //   //     testKernels[i](mesh->Nelements,
  //   //             mesh->o_sgeo,
  //   //             mesh->o_intInterpT,
  //   //             mesh->o_intLIFTT,
  //   //             mesh->o_vmapM,
  //   //             mesh->o_vmapP,
  //   //             mesh->o_EToB,
  //   //             time,
  //   //             mesh->o_intx,
  //   //             mesh->o_inty,
  //   //             o_UM,
  //   //             o_UP,
  //   //             o_UdM,
  //   //             o_UdP,
  //   //             o_X,
  //   //             o_Y);
  //   //   }

  //   //   occa::streamTag end = mesh->device.tagStream();
  //   //   mesh->device.finish();  
  //   //   toc = MPI_Wtime();
  //   //   //      kernelElapsed    = toc-tic;
  //   //   kernelElapsed = mesh->device.timeBetween(start,end);

  //   // }

  //   // else{
  //   //occaTimerTic(mesh->device,"KernelTime");
  //   tic = MPI_Wtime();  

  //   occa::streamTag start = mesh->device.tagStream();

  //     // assume 1 mpi process
  //     for(int it=0;it<iterations;++it){
  //       //printf("Cubature Points: %d", mesh->cubNp);
  //       testKernels[i](mesh->Nelements,
  //               mesh->o_sgeo,
  //               mesh->o_intInterpT,
  //               mesh->o_intLIFTT,
  //               mesh->o_vmapM,
  //               mesh->o_vmapP,
  //               mesh->o_EToB,
  //               time,
  //               mesh->o_intx,
  //               mesh->o_inty,
  //               o_U,
  //               o_V,
  //               o_Ud,
  //               o_Vd,
  //               o_X,
  //               o_Y);
  //     }

  //     occa::streamTag end = mesh->device.tagStream();
  //     mesh->device.finish();  
  //     toc = MPI_Wtime();
  //     //      kernelElapsed    = toc-tic;
  //     kernelElapsed = mesh->device.timeBetween(start,end);
  //   // }



      
  //     dfloat nmt = 1; 

  //     // // TW - why ?
  //     //  if(i==4 || i==5){
  //     //   nmt = NMT[mesh->N-1];
  //     //   // nmt = cNmtS;
  //     //  }


  //     dfloat alpha =0.0; 


  //      // if(i==8){
  //      //  Nbytes        = (sizeof(dfloat)*(8*intNtfp + 4*intNtfp*alpha + 4*mesh->Nfaces*(1.0-alpha) + 4*mesh->Np) 
  //      //                  +sizeof(int)*(0*Ntfp   +   intNtfp*alpha +   mesh->Nfaces*(1.0-alpha)))/2;

  //      //  NbytesShared  = (sizeof(dfloat)*(0*Ntfp + intNtfp*Nfp*alpha + 2*intNtfp + 2*Np*intNtfp)); 
  //      //  NbytesShared2 = (sizeof(dfloat)*(0*Ntfp + intNtfp*Nfp*alpha + 2*intNtfp + 2*Np*intNtfp 
  //      //                                      +intNtfp*Nfp/nmt + Np*intNtfp/nmt ));

  //      //  flops = intNtfp*( 0*Nfp*16 + 6 + 1 + 28) + Np*intNtfp*4; 

  //      // }
  //      // else{
          
  //      Nbytes           = (sizeof(dfloat)*(8*Ntfp + 4*intNtfp*alpha + 4*mesh->Nfaces*(1.0-alpha) + 4*mesh->Np) 
  //                         +sizeof(int)*(2*Ntfp   +   intNtfp*alpha +   mesh->Nfaces*(1.0-alpha)))/2;

  //      NbytesShared     = (sizeof(dfloat)*(8*Ntfp + 8*intNtfp*Nfp + 2*intNtfp + 2*Np*intNtfp)); 

  //      NbytesShared2     = (sizeof(dfloat)*(8*Ntfp + 8*intNtfp*Nfp + 2*intNtfp + 2*Np*intNtfp 
  //                                                  +intNtfp*Nfp/nmt + Np*intNtfp/nmt ));

  //      // NbytesShared2     = (sizeof(dfloat)*(8*Ntfp + 4*mesh->Nfaces + 
		// 			  //                                 8*intNtfp*Nfp + 
		// 		   //                                   intNtfp*Nfp + 
		// 			  //                                  2*intNtfp +
		// 			  //                                  2*Np*intNtfp +
		// 			  //                                    Np*intNtfp
		// 			  //                                    ) );

  //      flops            = intNtfp*( Nfp*16 + 6 + 1 + 28) + Np*intNtfp*4;
  //     // }

      
            
  //     occa::memory o_foo = mesh->device.malloc(Nbytes*mesh->Nelements);
  //     occa::memory o_bah = mesh->device.malloc(Nbytes*mesh->Nelements);

  //     mesh->device.finish(); 
  //     tic = MPI_Wtime();

  //     occa::streamTag startCopy = mesh->device.tagStream();
  //     for(int it=0;it<iterations;++it){
  //        o_bah.copyTo(o_foo);
  //     }
  //     occa::streamTag endCopy = mesh->device.tagStream();

  //     mesh->device.finish();
  //     toc = MPI_Wtime();
  //     //      double copyElapsed = (toc-tic);
  //     double copyElapsed = mesh->device.timeBetween(startCopy, endCopy);

  //     // Compute Data
  //     double copyBandwidth = mesh->Nelements*((Nbytes*iterations*2)/(1e9*copyElapsed));
  //     double  bw           = mesh->Nelements*((Nbytes*iterations*2)/(1e9*kernelElapsed));

  //     double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);
  //     double d2dbound      = copyBandwidth*gflops/bw;

  //     double smbound       = 7882*flops/( (double) NbytesShared);
  
  //     double intensity     = gflops/bw; 

  //     double l1bound      = 7882*flops/( (double) NbytesShared2);

  //     double roofline1     = mymin(d2dbound, smbound);

  //     double roofline2     = mymin(mymin(d2dbound, smbound),l1bound);

  //     double max_thg_p100 = 4670; 
  //     double ach_thg      = mymin(549*intensity, max_thg_p100);



  //     printf("[ N\tK\tDOFS\tKernelTime\tCopyTime\tIntensity\tGFLOPS/s\t(d2dbound GFLOPS/s)\tBW(GB/s)\t(SMBOUND GFLOPS/s)\t(ROOFLINE GFLOPS/s)\t(ROOFLINE2 GFLOPS/s)\tTH_peak]\n");
  //     printf("%02d \t%02d\t%02d\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\n",
  //             mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
  //      kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline1, roofline2, ach_thg);

  //     char fname[BUFSIZ];
  //     sprintf(fname, "KernelData.dat");
  //     FILE *fp;
  //     fp = fopen(fname, "a");

  //     fprintf(fp, "%02d %02d %02d %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E\n",
  //             mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
  //             kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline1, roofline2, ach_thg);
  //     fclose(fp);


  //   }

  // #endif


  // #if KERNEL_TEST==3
  // int NKernels = 5;

  // occa::kernel *testKernels = new occa::kernel[NKernels];
  // char kernelNames[NKernels][BUFSIZ];

  // for(int i=0; i<NKernels; i++)
  // {
    
  //   sprintf(kernelNames[i], "ellipticPartialAxIpdgTri2D_v%d", i);


  //   testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/optAxIpdgTri2D.okl",kernelNames[i], kernelInfo);
  //   printf("Ax Kernel #%02d\n", i);
  //   printf("N: %d Np: %d Nfp: %d\n", mesh->N, mesh->Np, mesh->Nfaces*mesh->Nfp);


  //   // sync processes
  //   mesh->device.finish();
  //   //    MPI_Barrier(MPI_COMM_WORLD);

  //   //occaTimerTic(mesh->device,"KernelTime");
  //   tic = MPI_Wtime();  

  //   occa::streamTag start = mesh->device.tagStream();

  //     // assume 1 mpi process
  //     for(int it=0;it<iterations;++it){
  //       testKernels[i](mesh->NinternalElements,
  //                                 mesh->o_internalElementIds,
  //                                 mesh->o_vmapM,
  //                                 mesh->o_vmapP,
  //                                 lambda, // lamda = 0.0
  //                                 tau,
  //                                 mesh->o_vgeo,
  //                                 mesh->o_sgeo,
  //                                 mesh->o_EToB,
  //                                 mesh->o_DrT,
  //                                 mesh->o_DsT,
  //                                 mesh->o_LIFTT,
  //                                 mesh->o_MM,
  //                                 o_G,
  //                                 o_X);
      
  //     }

  //     occa::streamTag end = mesh->device.tagStream();
  //     mesh->device.finish();  
  //     toc = MPI_Wtime();
  //     //      kernelElapsed    = toc-tic;
  //     kernelElapsed = mesh->device.timeBetween(start,end);
      
  //     // const dfloat alpha = 1.0;  // assumption on caching element/face values
  //      const dfloat alpha = 0.0;
  //     // const dfloat alpha = -1.0/9.0 + 1.0/9.0*mesh->N;

  //      if(i==3){

  //        Nbytes = (sizeof(dfloat)*(4*Np    // float4
  //                                 +4*Ntfp  // float4
  //                                 +alpha*5*Nfaces + (1-alpha)*5*Ntfp   // alpha = 1 means cached face geometric factors
  //                                 +alpha*5 +(1-alpha)*5*Np   // alpha 1 means drdx, J etc.. are all cached for an element
  //                                 +1*mesh->Np))/2           // gloabl write of output
  //             + (sizeof(int)*(alpha*1 + (1-alpha)*Np*1  // element list//  was 3
  //                             +2*Ntfp  // idM, idP
  //                             +alpha*Nfaces + (1-alpha)*Ntfp // BC's
  //                             ))/2;  



  //      }
  //      else{
  //      Nbytes = (sizeof(dfloat)*(4*Np    // float4
  //                              +8*Ntfp  // float4
  //                              +alpha*5*Nfaces + (1-alpha)*5*Ntfp   // alpha = 1 means cached face geometric factors
  //                              +alpha*5 +(1-alpha)*5*Np   // alpha 1 means drdx, J etc.. are all cached for an element
  //                              +1*mesh->Np))/2           // gloabl write of output
  //             + (sizeof(int)*(alpha*1 + (1-alpha)*Np*1  // element list//  was 3
  //                             +2*Ntfp  // idM, idP
  //                             +alpha*Nfaces + (1-alpha)*Ntfp // BC's
  //                             ))/2;  

  //      }

     

  //      NbytesShared = (sizeof(dfloat)*(3*Np + 3*Ntfp + Np*(2*Ntfp +2+4) + Ntfp*3 + Np*(2*Np + 1) + Np*(1*Ntfp + 1) + Np*(1*Np))); 

  //      flops = 1*Np + 21*Ntfp + Np*(4*Ntfp + 8) + 6*Ntfp+ Np*(4*Np + 2)+ Np*(1+2*Ntfp) + Np*(1+2*Np); 



  //     occa::memory o_foo = mesh->device.malloc(Nbytes*mesh->Nelements);
  //     occa::memory o_bah = mesh->device.malloc(Nbytes*mesh->Nelements);

  //     mesh->device.finish(); 
  //     tic = MPI_Wtime();

  //     occa::streamTag startCopy = mesh->device.tagStream();
  //     for(int it=0;it<iterations;++it){
  //        o_bah.copyTo(o_foo);
  //     }
  //     occa::streamTag endCopy = mesh->device.tagStream();

  //     mesh->device.finish();
  //     toc = MPI_Wtime();
  //     //      double copyElapsed = (toc-tic);
  //     double copyElapsed = mesh->device.timeBetween(startCopy, endCopy);


  //     // Compute Data
  //     double copyBandwidth = mesh->Nelements*((Nbytes*iterations*2)/(1e9*copyElapsed));
  //     double  bw           = mesh->Nelements*((Nbytes*iterations*2)/(1e9*kernelElapsed));

  //     double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);
  //     double d2dbound      = copyBandwidth*gflops/bw;

  //     double smbound       = 7882*flops/( (double) NbytesShared);
  
  //     double intensity    = gflops/bw; 

  //     double roofline     = mymin(d2dbound, smbound);

  //     double max_thg_p100 = 4670; 
  //     double ach_thg      = mymin(549*intensity, max_thg_p100);



  //     printf("[ N\tK\tDOFS\tKernelTime\tCopyTime\tIntensity\tGFLOPS/s\t(d2dbound GFLOPS/s)\tBW(GB/s)\t(SMBOUND GFLOPS/s)\t(ROOFLINE GFLOPS/s)\tTH_peak]\n");
  //     printf("%02d \t%02d\t%02d\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\n",
  //             mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
  //      kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline, ach_thg);

  //     char fname[BUFSIZ];
  //     sprintf(fname, "KernelData.dat");
  //     FILE *fp;
  //     fp = fopen(fname, "a");

  //     fprintf(fp, "%02d %02d %02d %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E\n",
  //             mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
  //             kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline, ach_thg);
  //     fclose(fp);

  //   }

  // #endif

// #if KERNEL_TEST==4
//   int NKernels = 6;

//   occa::kernel *testKernels = new occa::kernel[NKernels];
//   char kernelNames[NKernels][BUFSIZ];

//   for(int i=0; i<NKernels; i++)
//   {

        
//     sprintf(kernelNames[i], "ellipticPartialGradientTri2D_v%d", i);


//     testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/optGradientTri2D.okl",kernelNames[i], kernelInfo);
//     printf("Gradient Kernel #%02d\n", i);
//     printf("N: %d Np: %d Nfp: %d\n", mesh->N, mesh->Np, mesh->Nfaces*mesh->Nfp);

   
//     // sync processes
//     mesh->device.finish();
//     //    MPI_Barrier(MPI_COMM_WORLD);


//     if(i==5){
//       //occaTimerTic(mesh->device,"KernelTime");
//     tic = MPI_Wtime();  

//     occa::streamTag start = mesh->device.tagStream();
    
//       // assume 1 mpi process
//       for(int it=0;it<iterations;++it){
//         testKernels[i](mesh->NinternalElements,
//                        zero, // zero offset           
//                       mesh->o_vgeo,
//                       mesh->o_DrT,
//                       mesh->o_DsT,
//                       o_X,
//                       o_G);
//       }

//       occa::streamTag end = mesh->device.tagStream();
//       mesh->device.finish();  
//       toc = MPI_Wtime();
//       //      kernelElapsed    = toc-tic;
//       kernelElapsed = mesh->device.timeBetween(start,end);
//     }
//     else{
//         //occaTimerTic(mesh->device,"KernelTime");
//     tic = MPI_Wtime();  

//     occa::streamTag start = mesh->device.tagStream();
    
//        // assume 1 mpi process
//       for(int it=0;it<iterations;++it){
//         testKernels[i](mesh->NinternalElements,
//                        zero, // zero offset           
//                       mesh->o_vgeo,
//                       o_DrsT,
//                       mesh->o_DsT,
//                       o_X,
//                       o_G);
//       }


//       occa::streamTag end = mesh->device.tagStream();
//       mesh->device.finish();  
//       toc = MPI_Wtime();
//       //      kernelElapsed    = toc-tic;
//       kernelElapsed = mesh->device.timeBetween(start,end);

//    }

   
      
//       const dfloat alpha = 1.0;  // assumption on caching element/face values
//       // const dfloat alpha = 0.0;
//       // const dfloat alpha = -1.0/9.0 + 1.0/9.0*mesh->N;

//       // if(i<=4){ 
//       //   Nbytes       = (sizeof(dfloat)*(1*Np +4*alpha*Np + 4*(1.0-alpha)*1 + 4*Np)/2);
//       //   NbytesShared = (sizeof(dfloat)*(1*Np + Np*(2*Np))); 
//       //   flops        = Np*(Np*4 +6);  // All float ops only
//       // }


//       if(i==5){
//         Nbytes       = (sizeof(dfloat)*(1*Np +4*alpha*Np + 4*(1.0-alpha)*1 + 4*Np)/2);
//         NbytesShared = (sizeof(dfloat)*(1*Np + Np*(2*Np) + 4 + 4*Np + Np)); 
//         flops        = Np*(Np*4 +6);  // All float ops only
//       }

//       else{
//         Nbytes       = (sizeof(dfloat)*(1*Np +4*alpha*Np + 4*(1.0-alpha)*1 + 4*Np)/2);
//         NbytesShared = (sizeof(dfloat)*(1*Np + Np*(2*Np) + 1*Np)); 
//         flops        = Np*(Np*4 +6);  // All float ops only
//       }



//       occa::memory o_foo = mesh->device.malloc(Nbytes*mesh->Nelements);
//       occa::memory o_bah = mesh->device.malloc(Nbytes*mesh->Nelements);

//       mesh->device.finish(); 
//       tic = MPI_Wtime();

//       occa::streamTag startCopy = mesh->device.tagStream();
//       for(int it=0;it<iterations;++it){
//          o_bah.copyTo(o_foo);
//       }
//       occa::streamTag endCopy = mesh->device.tagStream();

//       mesh->device.finish();
//       toc = MPI_Wtime();
//       //      double copyElapsed = (toc-tic);
//       double copyElapsed = mesh->device.timeBetween(startCopy, endCopy);


//       // Compute Data
//       double copyBandwidth = mesh->Nelements*((Nbytes*iterations*2)/(1e9*copyElapsed));
//       double  bw           = mesh->Nelements*((Nbytes*iterations*2)/(1e9*kernelElapsed));

//       double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);
//       double d2dbound      = copyBandwidth*gflops/bw;

//       double smbound       = 7882*flops/( (double) NbytesShared); // 7751 for 1189 GHz
  
//       double intensity    = gflops/bw; 

//       double roofline     = mymin(d2dbound, smbound);

//       double max_thg_p100 = 4670; 
//       double ach_thg      = mymin(549*intensity, max_thg_p100);



//       printf("[ N\tK\tDOFS\tKernelTime\tCopyTime\tIntensity\tGFLOPS/s\t(d2dbound GFLOPS/s)\tBW(GB/s)\t(SMBOUND GFLOPS/s)\t(ROOFLINE GFLOPS/s)\tTH_peak]\n");
//       printf("%02d \t%02d\t%02d\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\n",
//               mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
//        kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline, ach_thg);

//       char fname[BUFSIZ];
//       sprintf(fname, "KernelData.dat");
//       FILE *fp;
//       fp = fopen(fname, "a");

//       fprintf(fp, "%02d %02d %02d %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E\n",
//               mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
//               kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline, ach_thg);
//       fclose(fp);

//     }

//   #endif

 


// #if KERNEL_TEST==5

  
//   int Nbl      = 20;
//   int Nmult    = 10;
//   int NKernels = Nbl*Nmult;

//   occa::kernel *testKernels = new occa::kernel[NKernels];
//   char kernelNames[NKernels][BUFSIZ];

//   dfloat mintime = 100.f;

//   int mintblock = 0; 
//   int mintmult  = 0; 

//   // for(int i=6; i<NKernels; i++)
//   // {
  
//   int i = 0; 
//   for (int b=1;b<=Nbl; b++){

//     for(int m =1; m<=Nmult; m++){

//     occa::kernelInfo kernelInfoT  = kernelInfo;
//     sprintf(kernelNames[i], "insSubCycleCubatureVolume2D_v5");

//     // for kernel 6
//     // kernelInfoT.addDefine("p_NbV", b);
//     // kernelInfoT.addDefine("p_Nmt", m);
    
//     // for kernel 3
//     //kernelInfoT.addDefine("p_cubNblockV",b);

//      kernelInfoT.addDefine("p_cNbV", b);
//      kernelInfoT.addDefine("p_cNmt", m);
   
   
//     testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",kernelNames[i], kernelInfoT);
//     printf("insSubCycleCubatureVolume Kernel #%02d\n", i);
//     printf("Nblock: %d cubNblock: %d N: %d Np: %d cubNp: %d\n", NblockV, cubNblockV, mesh->N, mesh->Np, mesh->cubNp);


//     // sync processes
//     mesh->device.finish();
//     //    MPI_Barrier(MPI_COMM_WORLD);

//     //occaTimerTic(mesh->device,"KernelTime");
//     tic = MPI_Wtime();  

//     occa::streamTag start = mesh->device.tagStream();

//       // assume 1 mpi process
//       for(int it=0;it<iterations;++it){
//         //printf("Cubature Points: %d", mesh->cubNp);
//         testKernels[i](mesh->Nelements,
//                   mesh->o_vgeo,
//                   mesh->o_cubDrWT,
//                   mesh->o_cubDsWT,
//                   mesh->o_cubInterpT,
//                   o_U,
//                   o_V,
//                   o_Ud,
//                   o_Vd,
//                   o_X,
//                   o_Y);
//       }

//       occa::streamTag end = mesh->device.tagStream();
//       mesh->device.finish();  
//       toc = MPI_Wtime();
//       //      kernelElapsed    = toc-tic;
//       kernelElapsed = mesh->device.timeBetween(start,end);

    
          
//         Nbytes       = (sizeof(dfloat)*(4*Np +4*Np + 2*Np)/2);

//         NbytesShared = (sizeof(dfloat)*(4*Nc + 4*Np*Nc + 4*Nc + 4*Np*Nc)); 
//         flops        = Np*6 + Np*Nc*8 + 4*Nc + 8*Np*Nc + 2*Np ;  // All float ops only

      
//       if(kernelElapsed<mintime){
//         mintime = kernelElapsed;
//         mintblock = b;
//         mintmult  = m; 

//       }


      
//       double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);

//       printf("[ N\tBlock\tNmult\tKernelTime\tGFLOPS/s\t mintime]\n");
//       printf("%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed, gflops,  mintime, mintblock, mintmult);
     
//      char fname[BUFSIZ];
//       sprintf(fname, "KernelOptimization.dat");
//       FILE *fp;
//       fp = fopen(fname, "a");

//       fprintf(fp,"%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed, gflops,  mintime, mintblock, mintmult);
//       fclose(fp);

//       i++;


//      }
//     }

//     // }

//   #endif


// #if KERNEL_TEST==6 // Optimize surface kernel

  
//   int Nbl      = 10;
//   int Nmult    = 5;
//   int NKernels = Nbl*Nmult;

//   occa::kernel *testKernels = new occa::kernel[NKernels];
//   char kernelNames[NKernels][BUFSIZ];

//   dfloat mintime = 100.f;

//   int mintblock = 0; 
//   int mintmult  = 0; 
 
//   int i = 0; 
//   for (int b=1;b<=Nbl; b++){

//     for(int m =1; m<=Nmult; m++){

//     occa::kernelInfo kernelInfoT  = kernelInfo;
//     sprintf(kernelNames[i], "insSubCycleCubatureSurface2D_v8");

    
//      // kernelInfoT.addDefine("p_cubNblockS", b);
//     kernelInfoT.addDefine("p_NnodesS5", m);
//     kernelInfoT.addDefine("p_NblockS5", b);
   
   
//     testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/optSubcycleSurface2D.okl",kernelNames[i], kernelInfoT);
//     printf("insSubCycleCubatureSurface Kernel #%02d\n", i);
//     printf("NblockS: %d cubNblockS: %d N: %d Np: %d Nfp: %d cubNfp: %d\n", NblockS, cubNblockS,  mesh->N, mesh->Np, mesh->Nfaces*mesh->Nfp,  mesh->Nfaces*mesh->intNfp );


//     // sync processes
//     mesh->device.finish();
//     // MPI_Barrier(MPI_COMM_WORLD);


//            //occaTimerTic(mesh->device,"KernelTime");
//     tic = MPI_Wtime();  

//     occa::streamTag start = mesh->device.tagStream();

//       // assume 1 mpi process
//       for(int it=0;it<iterations;++it){
//         //printf("Cubature Points: %d", mesh->cubNp);
//         testKernels[i](mesh->Nelements,
//                 mesh->o_sgeo,
//                 mesh->o_intInterpT,
//                 mesh->o_intLIFTT,
//                 mesh->o_vmapM,
//                 mesh->o_vmapP,
//                 mesh->o_EToB,
//                 time,
//                 mesh->o_intx,
//                 mesh->o_inty,
//                 o_U,
//                 o_V,
//                 o_Ud,
//                 o_Vd,
//                 o_X,
//                 o_Y);
//       }

//       occa::streamTag end = mesh->device.tagStream();
//       mesh->device.finish();  
//       toc = MPI_Wtime();

//       //      kernelElapsed    = toc-tic;
//       kernelElapsed = mesh->device.timeBetween(start,end);
    




  
//        Nbytes           = (sizeof(dfloat)*(8*Ntfp + 4*intNtfp + 4*mesh->Np) + sizeof(int)*(2*Ntfp))/2;

//        NbytesShared     = (sizeof(dfloat)*(8*Ntfp + 8*intNtfp*Nfp + 2*intNtfp + 2*Np*intNtfp)); 

//        flops            = intNtfp*( Nfp*16 + 6 + 1 + 28) + Np*intNtfp*4;

    
      
//       if(kernelElapsed<mintime){
//         mintime = kernelElapsed;
//         mintblock = b;
//         mintmult  = m; 

//       }


      
//       double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);

//       printf("[ N\tBlock\tNmult\tKernelTime\tGFLOPS/s\t mintime]\n");
//       printf("%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed/iterations, gflops,  mintime, mintblock, mintmult);
     
//      char fname[BUFSIZ];
//       sprintf(fname, "KernelOptimization.dat");
//       FILE *fp;
//       fp = fopen(fname, "a");

//       fprintf(fp,"%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed/iterations, gflops,  mintime, mintblock, mintmult);
//       fclose(fp);

//       i++;


//      }
//     }

//     // }

//   #endif




  // #if KERNEL_TEST==7 // Optimize surface kernel

  
  // int Nbl      = 10;
  // int Nmult    = 1;
  // int NKernels = Nbl*Nmult;

  // occa::kernel *testKernels = new occa::kernel[NKernels];
  // char kernelNames[NKernels][BUFSIZ];

  // dfloat mintime = 100.f;

  // int mintblock = 0; 
  // int mintmult  = 0; 


  // int i = 0; 
  // for (int b=1;b<=Nbl; b++){

  //   for(int m =1; m<=Nmult; m++){

  //   occa::kernelInfo kernelInfoT  = kernelInfo;
  //   sprintf(kernelNames[i], "ellipticPartialAxIpdgTri2D_v4");

  //    //   kernelInfoT.addDefine("p_Nblock", b);
  //   kernelInfoT.addDefine("p_Nmt", m);
  //   kernelInfoT.addDefine("p_NbV", b);
   
   
  //   testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/optAxIpdgTri2D.okl",kernelNames[i], kernelInfoT);
  //   printf("Ax Kernel #%02d\n", i);
  //   printf("N: %d Np: %d Nfp: %d\n", mesh->N, mesh->Np, mesh->Nfaces*mesh->Nfp);


  //   // sync processes
  //   mesh->device.finish();
  //   //    MPI_Barrier(MPI_COMM_WORLD);

  //   //occaTimerTic(mesh->device,"KernelTime");
  //   tic = MPI_Wtime();  

  //   occa::streamTag start = mesh->device.tagStream();

  //     // assume 1 mpi process
  //     for(int it=0;it<iterations;++it){
  //       testKernels[i](mesh->NinternalElements,
  //                                 mesh->o_internalElementIds,
  //                                 mesh->o_vmapM,
  //                                 mesh->o_vmapP,
  //                                 lambda, // lamda = 0.0
  //                                 tau,
  //                                 mesh->o_vgeo,
  //                                 mesh->o_sgeo,
  //                                 mesh->o_EToB,
  //                                 mesh->o_DrT,
  //                                 mesh->o_DsT,
  //                                 mesh->o_LIFTT,
  //                                 mesh->o_MM,
  //                                 o_G,
  //                                 o_X);
      
  //     }

  //     occa::streamTag end = mesh->device.tagStream();
  //     mesh->device.finish();  
  //     toc = MPI_Wtime();
  //     //      kernelElapsed    = toc-tic;
  //     kernelElapsed = mesh->device.timeBetween(start,end);


  //     const dfloat alpha = 1.0;  // assumption on caching element/face values
  //     // const dfloat alpha = 0.0;
  //     // const dfloat alpha = -1.0/9.0 + 1.0/9.0*mesh->N;

  //     Nbytes = (sizeof(dfloat)*(4*Np    // float4
  //                              +8*Ntfp // float4
  //                              +alpha*5*Nfaces + (1-alpha)*5*Ntfp   // alpha = 1 means cached face geometric factors
  //                              +alpha*5 +(1-alpha)*5*Np   // alpha 1 means drdx, J etc.. are all cached for an element
  //                              +1*mesh->Np))/2           // gloabl write of output
  //             + (sizeof(int)*(alpha*1 + (1-alpha)*Np*1  // element list
  //                             +2*Ntfp  // idM, idP
  //                             +alpha*Nfaces + (1-alpha)*Ntfp // BC's
  //                             ))/2;  

  //      NbytesShared = (sizeof(dfloat)*(3*Np + 3*Ntfp + Np*(2*Ntfp +2+4) + Ntfp*3 + Np*(2*Np + 1) + Np*(1*Ntfp + 1) + Np*(1*Np))); 

  //      flops = 1*Np + 21*Ntfp + Np*(4*Ntfp + 8) + 6*Ntfp+ Np*(4*Np + 2)+ Np*(1+2*Ntfp) + Np*(1+2*Np); 

      
  //     if(kernelElapsed<mintime){
  //       mintime = kernelElapsed;
  //       mintblock = b;
  //       mintmult  = m; 

  //     }


      
  //     double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);

  //     printf("[ N\tBlock\tNmult\tKernelTime\tGFLOPS/s\t mintime]\n");
  //     printf("%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed/iterations, gflops,  mintime, mintblock, mintmult);
     
  //    char fname[BUFSIZ];
  //     sprintf(fname, "KernelOptimization.dat");
  //     FILE *fp;
  //     fp = fopen(fname, "a");

  //     fprintf(fp,"%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed/iterations, gflops,  mintime, mintblock, mintmult);
  //     fclose(fp);

  //     i++;


  //    }
  //   }

  //   // }

  // #endif

  
 // #if KERNEL_TEST==8 // Optimize surface kernel

  
 //  int Nbl      = 1;
 //  int Nmult    = 1;
 //  int NKernels = Nbl*Nmult;

 //  occa::kernel *testKernels = new occa::kernel[NKernels];
 //  char kernelNames[NKernels][BUFSIZ];

 //  dfloat mintime = 100.f;

 //  int mintblock = 0; 
 //  int mintmult  = 0; 


 //  int i = 0; 
 //  for (int b=1;b<=Nbl; b++){

 //    for(int m =1; m<=Nmult; m++){

 //    occa::kernelInfo kernelInfoT  = kernelInfo;
 //    sprintf(kernelNames[i], "ellipticPartialGradientTri2D_v5");

 //     //   kernelInfoT.addDefine("p_Nblock", b);
 //    kernelInfoT.addDefine("p_Nmt2", m);
 //    kernelInfoT.addDefine("p_NbV2", b);
   
   
 //    testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/optGradientTri2D.okl",kernelNames[i], kernelInfoT);
 //    printf("Gradient Kernel #%02d\n", i);
 //    printf("N: %d Np: %d Nfp: %d\n", mesh->N, mesh->Np, mesh->Nfaces*mesh->Nfp);


 //    // sync processes
 //    mesh->device.finish();
 //    //    MPI_Barrier(MPI_COMM_WORLD);

 //    //occaTimerTic(mesh->device,"KernelTime");
 //    tic = MPI_Wtime();  

 //    occa::streamTag start = mesh->device.tagStream();


 //     #if 1
 //      // assume 1 mpi process
 //      for(int it=0;it<iterations;++it){
 //        testKernels[i](mesh->NinternalElements,
 //                       zero, // zero offset           
 //                      mesh->o_vgeo,
 //                      mesh->o_DrT,
 //                      mesh->o_DsT,
 //                      o_X,
 //                      o_G);

 //      }
    
 //     #else
 //       // assume 1 mpi process
 //      for(int it=0;it<iterations;++it){
 //        testKernels[i](mesh->NinternalElements,
 //                       zero, // zero offset           
 //                      mesh->o_vgeo,
 //                      o_DrsT,
 //                      mesh->o_DsT,
 //                      o_X,
 //                      o_G);
 //      }
 //      #endif





 //      occa::streamTag end = mesh->device.tagStream();
 //      mesh->device.finish();  
 //      toc = MPI_Wtime();
 //      //      kernelElapsed    = toc-tic;
 //      kernelElapsed = mesh->device.timeBetween(start,end);


 //      const dfloat alpha = 1.0;  // assumption on caching element/face values
 //      // const dfloat alpha = 0.0;
 //      // const dfloat alpha = -1.0/9.0 + 1.0/9.0*mesh->N;

      
 //      Nbytes       = (sizeof(dfloat)*(1*Np +4*alpha*Np + 4*(1.0-alpha)*1 + 4*Np)/2);
 //      NbytesShared = (sizeof(dfloat)*(1*Np + Np*(2*Np))); 
 //      flops        = Np*(Np*4 +6);  // All float ops only
      
 //      if(kernelElapsed<mintime){
 //        mintime = kernelElapsed;
 //        mintblock = b;
 //        mintmult  = m; 

 //      }


      
 //      double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);

 //      printf("[ N\tBlock\tNmult\tKernelTime\tGFLOPS/s\t mintime]\n");
 //      printf("%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed/iterations, gflops,  mintime, mintblock, mintmult);
     
 //      char fname[BUFSIZ];
 //      sprintf(fname, "KernelOptimization.dat");
 //      FILE *fp;
 //      fp = fopen(fname, "a");

 //      fprintf(fp,"%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed/iterations, gflops,  mintime, mintblock, mintmult);
 //      fclose(fp);

 //      i++;


 //     }
 //    }

 //    // }

 //  #endif
