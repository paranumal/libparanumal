#include <mpi.h>
#include "mesh2D.h"
#include "ellipticTri2D.h"

// 1 Advection Volume 2 Advection Surface 3 Ax  4 Gradient
// for optimization
// 5 advaction volume 6 advection surface 7 Ax kernel 8 Gradient
#define KERNEL_TEST 7

void insRunTimer2D(mesh2D *mesh, char *options, char *boundaryHeaderFileName){

  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  mesh->Nfields = 1;

  
  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", (rank)%2);

  occa::kernelInfo kernelInfo;
  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);
  
  kernelInfo.addInclude(boundaryHeaderFileName);


  iint index = 0, iterations = 100,  Nbytes=0,  zero = 0;  
  dfloat lambda = 0.0, tau = 1.0; 
  dfloat time   = 0.0; 
  iint  Ntotal     = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
  iint  cubNtotal  = (mesh->Nelements+mesh->totalHaloPairs)*mesh->cubNp;
  iint  Nftotal    = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfaces*mesh->Nfp;
  iint  intNftotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfaces*mesh->intNfp;

  dfloat *Z      = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  dfloat *ZM     = (dfloat*) calloc(Nftotal,sizeof(dfloat));
  dfloat *ZM1    = (dfloat*) calloc(intNftotal,sizeof(dfloat));
  dfloat *ZM4    = (dfloat*) calloc(4*Nftotal,sizeof(dfloat));
  dfloat *ZM2    = (dfloat*) calloc(2*Nftotal,sizeof(dfloat));


 // Test for dfloat2 in gradient kernel
  dfloat *DrsT  = (dfloat*) calloc(mesh->Np*mesh->Np*2, sizeof(dfloat));


  dfloat *G     = (dfloat*) calloc(Ntotal*4, sizeof(dfloat));
  

  occa::memory o_U, o_V, o_X, o_Y,  o_Ud, o_Vd, o_G, o_DrsT;
  occa::memory o_UM, o_UP, o_UdM, o_UdP;
  //occa::memory o_UT, o_VTd; 
  o_U   = mesh->device.malloc(Ntotal*sizeof(dfloat),Z);
  o_V   = mesh->device.malloc(Ntotal*sizeof(dfloat),Z);
  o_Ud  = mesh->device.malloc(Ntotal*sizeof(dfloat),Z);
  o_Vd  = mesh->device.malloc(Ntotal*sizeof(dfloat),Z);

  o_X   = mesh->device.malloc(Ntotal*sizeof(dfloat),Z);
  o_Y   = mesh->device.malloc(Ntotal*sizeof(dfloat),Z);

  o_G   = mesh->device.malloc(4*Ntotal*sizeof(dfloat),G);

  o_DrsT = mesh->device.malloc(2*mesh->Np*mesh->Np, DrsT);


  o_UM   = mesh->device.malloc(2*Nftotal*sizeof(dfloat),ZM2);
  o_UP   = mesh->device.malloc(2*Nftotal*sizeof(dfloat),ZM2);
  o_UdM  = mesh->device.malloc(2*Nftotal*sizeof(dfloat),ZM2);
  o_UdP  = mesh->device.malloc(2*Nftotal*sizeof(dfloat),ZM2);




  
  free(Z);  free(ZM);  free(ZM1); free(ZM2); free(ZM4); free(G); // free(ZM2);

 
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 512/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);
  
  // TW: why using 128 here ?
  iint maxNodesVolumeCub = mymax(mesh->cubNp,mesh->Np);  
  kernelInfo.addDefine("p_maxNodesVolumeCub", maxNodesVolumeCub);
  iint cubNblockV = mymax(1,1024/maxNodesVolumeCub); 
  //
  iint maxNodesSurfaceCub = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  
  kernelInfo.addDefine("p_maxNodesSurfaceCub",maxNodesSurfaceCub);

  iint cubNblockS = mymax(1,1024/maxNodesSurfaceCub);
    
  double flops;
 
  iint Np      = mesh->Np;
  iint Nc      = mesh->cubNp;
  iint Nfp     = mesh->Nfp; 
  iint Ntfp    = mesh->Nfaces*mesh->Nfp; 
  iint intNtfp = mesh->Nfaces*mesh->intNfp;
  iint Nfaces  = mesh->Nfaces;
  double tic   = 0.0, toc = 0.0, kernelElapsed=0.0;
  int NbytesShared = 0;
  int NbytesShared2 = 0;   



  occa::kernel TestKernel; 

  #if KERNEL_TEST==1

  int NKernels = 7;
  // TW: what do these variables mean ?
  dfloat  NMT1[10] = {2,2,2,2,2,2,3,2,2,1};
  dfloat  NMT2[10] = {1,2,2,2,2,3,3,2,2,2};

  occa::kernel *testKernels = new occa::kernel[NKernels];
  char kernelNames[NKernels][BUFSIZ];

  for(iint i=0; i<NKernels; i++)
  {
    
    sprintf(kernelNames[i], "insSubCycleCubatureVolume2D_v%d", i);

    testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/optSubcycleVolume2D.okl",kernelNames[i], kernelInfo);
    printf("insSubCycleCubatureVolume Kernel #%02d\n", i);
    printf("Nblock: %d cubNblock: %d N: %d Np: %d cubNp: %d\n", NblockV, cubNblockV, mesh->N, mesh->Np, mesh->cubNp);


    // sync processes
    mesh->device.finish();
    //    MPI_Barrier(MPI_COMM_WORLD);

    //occaTimerTic(mesh->device,"KernelTime");
    tic = MPI_Wtime();  

    occa::streamTag start = mesh->device.tagStream();

      // assume 1 mpi process
      for(int it=0;it<iterations;++it){
        //printf("Cubature Points: %d", mesh->cubNp);
        testKernels[i](mesh->Nelements,
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
      }

      occa::streamTag end = mesh->device.tagStream();
      mesh->device.finish();  
      toc = MPI_Wtime();
      //      kernelElapsed    = toc-tic;
      kernelElapsed = mesh->device.timeBetween(start,end);

      dfloat nmt =1; 

      if(i==0){
        Nbytes       = (sizeof(dfloat)*(4*Np*Nc*0 +4*Np + 2*Np)/2);

        NbytesShared = (sizeof(dfloat)*(4*Nc + 4*Np*Nc)); 

        //NbytesShared2 = (sizeof(dfloat)*(4*Nc + 4*Np*Nc + 1.*Np*Nc/nmt + 2.*Np*Nc/nmt ));  // Add operators

        flops = Nc*Np*8 + 4*Nc + Np*Nc*16 + Np*14 + Np*2;  // All float ops only
      }
       else
      {
          
        Nbytes        = (sizeof(dfloat)*(4*Np +4*Np + 2*Np)/2);
        NbytesShared  = (sizeof(dfloat)*(4*Np + 4*Np*Nc + 4*Nc + 4*Np*Nc)); 
        NbytesShared2 = (sizeof(dfloat)*(4*Np + 4*Np*Nc + 4*Nc + 4*Np*Nc + 1.*Np*Nc/nmt + 2.*Np*Nc/nmt )); // Add operators
        flops         = Np*6 + Np*Nc*8 + 4*Nc + 8*Np*Nc + 2*Np ;  // All float ops only

      }


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
      //      double copyElapsed = (toc-tic);
      double copyElapsed = mesh->device.timeBetween(startCopy, endCopy);

     
      // TW: please use 10^9 for GFLOPS and GB

     // Compute Data
      //double peakSharedStreaming = 
      double copyBandwidth = mesh->Nelements*((Nbytes*iterations*2)/(1e9*copyElapsed));
      double  bw           = mesh->Nelements*((Nbytes*iterations*2)/(1e9*kernelElapsed));

      double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);
      double d2dbound      = copyBandwidth*gflops/bw;

      // TW: avoid using hard coded numbers like this in formula, define a new variable peakSharedStreaming
      //     7882 is defined wrt to 10^9
      double smbound       = 7882*flops/( (double) NbytesShared);
      double intensity     = gflops/bw; 
      double l1bound      = 7882*flops/( (double) NbytesShared2);

      printf("l1bound :%.2e \n", l1bound);



      double roofline1     = mymin(d2dbound, smbound);
      double roofline2     = mymin(d2dbound, l1bound);
      double max_thg_p100  = 4670; 
      double ach_thg       = mymin(549*intensity, max_thg_p100);



      printf("[ N\tK\tDOFS\tKernelTime\tCopyTime\tIntensity\tGFLOPS/s\t(d2dbound GFLOPS/s)\tBW(GB/s)\t(SMBOUND GFLOPS/s)\t(ROOFLINE GFLOPS/s)\t(ROOFLINE2 GFLOPS/s)\tTH_peak]\n");
      printf("%02d \t%02d\t%02d\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\n",
              mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
       kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline1, roofline2, ach_thg);

      char fname[BUFSIZ];
      sprintf(fname, "KernelData.dat");
      FILE *fp;
      fp = fopen(fname, "a");

      fprintf(fp, "%02d %02d %02d %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E\n",
              mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
              kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline1, roofline2, ach_thg);
      fclose(fp);


    }

  #endif


  



  // SURFACE KERNEL
  #if KERNEL_TEST==2
  
  int NKernels = 8;

  dfloat  NMT[10] = {2,1,2,2,2,2,2,2,3,2};

  occa::kernel *testKernels = new occa::kernel[NKernels];

  char kernelNames[NKernels][BUFSIZ];

  for(iint i=0; i<NKernels; i++)
  {
        
    sprintf(kernelNames[i], "insSubCycleCubatureSurface2D_v%d", i);

    testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/optSubcycleSurface2D.okl",kernelNames[i], kernelInfo);

    printf("insSubCycleCubatureSurface Kernel #%02d\n", i);
    printf("Nblock: %d cubNblock: %d N: %d Np: %d cubNp: %d\n", NblockV, cubNblockV, mesh->N, mesh->Np, mesh->cubNp);

    
     // sync processes
    mesh->device.finish();
    //    MPI_Barrier(MPI_COMM_WORLD);

    // if(i==7){

    // //occaTimerTic(mesh->device,"KernelTime");
    // tic = MPI_Wtime();  

    // occa::streamTag start = mesh->device.tagStream();

    //   // assume 1 mpi process
    //   for(int it=0;it<iterations;++it){
    //     //printf("Cubature Points: %d", mesh->cubNp);
    //     testKernels[i](mesh->Nelements,
    //             mesh->o_sgeo,
    //             mesh->o_intInterpT,
    //             mesh->o_intLIFTT,
    //             mesh->o_vmapM,
    //             mesh->o_vmapP,
    //             mesh->o_EToB,
    //             time,
    //             mesh->o_intx,
    //             mesh->o_inty,
    //             o_UM,
    //             o_UP,
    //             o_UdM,
    //             o_UdP,
    //             o_X,
    //             o_Y);
    //   }

    //   occa::streamTag end = mesh->device.tagStream();
    //   mesh->device.finish();  
    //   toc = MPI_Wtime();
    //   //      kernelElapsed    = toc-tic;
    //   kernelElapsed = mesh->device.timeBetween(start,end);

    // }

    // else{
    //occaTimerTic(mesh->device,"KernelTime");
    tic = MPI_Wtime();  

    occa::streamTag start = mesh->device.tagStream();

      // assume 1 mpi process
      for(int it=0;it<iterations;++it){
        //printf("Cubature Points: %d", mesh->cubNp);
        testKernels[i](mesh->Nelements,
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
      }

      occa::streamTag end = mesh->device.tagStream();
      mesh->device.finish();  
      toc = MPI_Wtime();
      //      kernelElapsed    = toc-tic;
      kernelElapsed = mesh->device.timeBetween(start,end);
    // }



      
      dfloat nmt = 1; 

      // // TW - why ?
      //  if(i==4 || i==5){
      //   nmt = NMT[mesh->N-1];
      //   // nmt = cNmtS;
      //  }


      dfloat alpha =0.0; 


       // if(i==8){
       //  Nbytes        = (sizeof(dfloat)*(8*intNtfp + 4*intNtfp*alpha + 4*mesh->Nfaces*(1.0-alpha) + 4*mesh->Np) 
       //                  +sizeof(iint)*(0*Ntfp   +   intNtfp*alpha +   mesh->Nfaces*(1.0-alpha)))/2;

       //  NbytesShared  = (sizeof(dfloat)*(0*Ntfp + intNtfp*Nfp*alpha + 2*intNtfp + 2*Np*intNtfp)); 
       //  NbytesShared2 = (sizeof(dfloat)*(0*Ntfp + intNtfp*Nfp*alpha + 2*intNtfp + 2*Np*intNtfp 
       //                                      +intNtfp*Nfp/nmt + Np*intNtfp/nmt ));

       //  flops = intNtfp*( 0*Nfp*16 + 6 + 1 + 28) + Np*intNtfp*4; 

       // }
       // else{
          
       Nbytes           = (sizeof(dfloat)*(8*Ntfp + 4*intNtfp*alpha + 4*mesh->Nfaces*(1.0-alpha) + 4*mesh->Np) 
                          +sizeof(iint)*(2*Ntfp   +   intNtfp*alpha +   mesh->Nfaces*(1.0-alpha)))/2;

       NbytesShared     = (sizeof(dfloat)*(8*Ntfp + 8*intNtfp*Nfp + 2*intNtfp + 2*Np*intNtfp)); 

       NbytesShared2     = (sizeof(dfloat)*(8*Ntfp + 8*intNtfp*Nfp + 2*intNtfp + 2*Np*intNtfp 
                                                   +intNtfp*Nfp/nmt + Np*intNtfp/nmt ));

       // NbytesShared2     = (sizeof(dfloat)*(8*Ntfp + 4*mesh->Nfaces + 
					  //                                 8*intNtfp*Nfp + 
				   //                                   intNtfp*Nfp + 
					  //                                  2*intNtfp +
					  //                                  2*Np*intNtfp +
					  //                                    Np*intNtfp
					  //                                    ) );

       flops            = intNtfp*( Nfp*16 + 6 + 1 + 28) + Np*intNtfp*4;
      // }

      
            
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
      //      double copyElapsed = (toc-tic);
      double copyElapsed = mesh->device.timeBetween(startCopy, endCopy);

      // Compute Data
      double copyBandwidth = mesh->Nelements*((Nbytes*iterations*2)/(1e9*copyElapsed));
      double  bw           = mesh->Nelements*((Nbytes*iterations*2)/(1e9*kernelElapsed));

      double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);
      double d2dbound      = copyBandwidth*gflops/bw;

      double smbound       = 7882*flops/( (double) NbytesShared);
  
      double intensity     = gflops/bw; 

      double l1bound      = 7882*flops/( (double) NbytesShared2);

      double roofline1     = mymin(d2dbound, smbound);

      double roofline2     = mymin(mymin(d2dbound, smbound),l1bound);

      double max_thg_p100 = 4670; 
      double ach_thg      = mymin(549*intensity, max_thg_p100);



      printf("[ N\tK\tDOFS\tKernelTime\tCopyTime\tIntensity\tGFLOPS/s\t(d2dbound GFLOPS/s)\tBW(GB/s)\t(SMBOUND GFLOPS/s)\t(ROOFLINE GFLOPS/s)\t(ROOFLINE2 GFLOPS/s)\tTH_peak]\n");
      printf("%02d \t%02d\t%02d\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\n",
              mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
       kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline1, roofline2, ach_thg);

      char fname[BUFSIZ];
      sprintf(fname, "KernelData.dat");
      FILE *fp;
      fp = fopen(fname, "a");

      fprintf(fp, "%02d %02d %02d %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E\n",
              mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
              kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline1, roofline2, ach_thg);
      fclose(fp);


    }

  #endif


  #if KERNEL_TEST==3
  int NKernels = 5;

  occa::kernel *testKernels = new occa::kernel[NKernels];
  char kernelNames[NKernels][BUFSIZ];

  for(iint i=0; i<NKernels; i++)
  {
    
    sprintf(kernelNames[i], "ellipticPartialAxIpdgTri2D_v%d", i);


    testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/optAxIpdgTri2D.okl",kernelNames[i], kernelInfo);
    printf("Ax Kernel #%02d\n", i);
    printf("N: %d Np: %d Nfp: %d\n", mesh->N, mesh->Np, mesh->Nfaces*mesh->Nfp);


    // sync processes
    mesh->device.finish();
    //    MPI_Barrier(MPI_COMM_WORLD);

    //occaTimerTic(mesh->device,"KernelTime");
    tic = MPI_Wtime();  

    occa::streamTag start = mesh->device.tagStream();

      // assume 1 mpi process
      for(int it=0;it<iterations;++it){
        testKernels[i](mesh->NinternalElements,
                                  mesh->o_internalElementIds,
                                  mesh->o_vmapM,
                                  mesh->o_vmapP,
                                  lambda, // lamda = 0.0
                                  tau,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  mesh->o_EToB,
                                  mesh->o_DrT,
                                  mesh->o_DsT,
                                  mesh->o_LIFTT,
                                  mesh->o_MM,
                                  o_G,
                                  o_X);
      
      }

      occa::streamTag end = mesh->device.tagStream();
      mesh->device.finish();  
      toc = MPI_Wtime();
      //      kernelElapsed    = toc-tic;
      kernelElapsed = mesh->device.timeBetween(start,end);
      
      // const dfloat alpha = 1.0;  // assumption on caching element/face values
       const dfloat alpha = 0.0;
      // const dfloat alpha = -1.0/9.0 + 1.0/9.0*mesh->N;

       if(i==3){

         Nbytes = (sizeof(dfloat)*(4*Np    // float4
                                  +4*Ntfp  // float4
                                  +alpha*5*Nfaces + (1-alpha)*5*Ntfp   // alpha = 1 means cached face geometric factors
                                  +alpha*5 +(1-alpha)*5*Np   // alpha 1 means drdx, J etc.. are all cached for an element
                                  +1*mesh->Np))/2           // gloabl write of output
              + (sizeof(iint)*(alpha*1 + (1-alpha)*Np*1  // element list//  was 3
                              +2*Ntfp  // idM, idP
                              +alpha*Nfaces + (1-alpha)*Ntfp // BC's
                              ))/2;  



       }
       else{
       Nbytes = (sizeof(dfloat)*(4*Np    // float4
                               +8*Ntfp  // float4
                               +alpha*5*Nfaces + (1-alpha)*5*Ntfp   // alpha = 1 means cached face geometric factors
                               +alpha*5 +(1-alpha)*5*Np   // alpha 1 means drdx, J etc.. are all cached for an element
                               +1*mesh->Np))/2           // gloabl write of output
              + (sizeof(iint)*(alpha*1 + (1-alpha)*Np*1  // element list//  was 3
                              +2*Ntfp  // idM, idP
                              +alpha*Nfaces + (1-alpha)*Ntfp // BC's
                              ))/2;  

       }

     

       NbytesShared = (sizeof(dfloat)*(3*Np + 3*Ntfp + Np*(2*Ntfp +2+4) + Ntfp*3 + Np*(2*Np + 1) + Np*(1*Ntfp + 1) + Np*(1*Np))); 

       flops = 1*Np + 21*Ntfp + Np*(4*Ntfp + 8) + 6*Ntfp+ Np*(4*Np + 2)+ Np*(1+2*Ntfp) + Np*(1+2*Np); 



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
      //      double copyElapsed = (toc-tic);
      double copyElapsed = mesh->device.timeBetween(startCopy, endCopy);


      // Compute Data
      double copyBandwidth = mesh->Nelements*((Nbytes*iterations*2)/(1e9*copyElapsed));
      double  bw           = mesh->Nelements*((Nbytes*iterations*2)/(1e9*kernelElapsed));

      double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);
      double d2dbound      = copyBandwidth*gflops/bw;

      double smbound       = 7882*flops/( (double) NbytesShared);
  
      double intensity    = gflops/bw; 

      double roofline     = mymin(d2dbound, smbound);

      double max_thg_p100 = 4670; 
      double ach_thg      = mymin(549*intensity, max_thg_p100);



      printf("[ N\tK\tDOFS\tKernelTime\tCopyTime\tIntensity\tGFLOPS/s\t(d2dbound GFLOPS/s)\tBW(GB/s)\t(SMBOUND GFLOPS/s)\t(ROOFLINE GFLOPS/s)\tTH_peak]\n");
      printf("%02d \t%02d\t%02d\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\n",
              mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
       kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline, ach_thg);

      char fname[BUFSIZ];
      sprintf(fname, "KernelData.dat");
      FILE *fp;
      fp = fopen(fname, "a");

      fprintf(fp, "%02d %02d %02d %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E\n",
              mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
              kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline, ach_thg);
      fclose(fp);

    }

  #endif

#if KERNEL_TEST==4
  int NKernels = 6;

  occa::kernel *testKernels = new occa::kernel[NKernels];
  char kernelNames[NKernels][BUFSIZ];

  for(iint i=0; i<NKernels; i++)
  {

        
    sprintf(kernelNames[i], "ellipticPartialGradientTri2D_v%d", i);


    testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/optGradientTri2D.okl",kernelNames[i], kernelInfo);
    printf("Gradient Kernel #%02d\n", i);
    printf("N: %d Np: %d Nfp: %d\n", mesh->N, mesh->Np, mesh->Nfaces*mesh->Nfp);

   
    // sync processes
    mesh->device.finish();
    //    MPI_Barrier(MPI_COMM_WORLD);


    if(i==5){
      //occaTimerTic(mesh->device,"KernelTime");
    tic = MPI_Wtime();  

    occa::streamTag start = mesh->device.tagStream();
    
      // assume 1 mpi process
      for(int it=0;it<iterations;++it){
        testKernels[i](mesh->NinternalElements,
                       zero, // zero offset           
                      mesh->o_vgeo,
                      mesh->o_DrT,
                      mesh->o_DsT,
                      o_X,
                      o_G);
      }

      occa::streamTag end = mesh->device.tagStream();
      mesh->device.finish();  
      toc = MPI_Wtime();
      //      kernelElapsed    = toc-tic;
      kernelElapsed = mesh->device.timeBetween(start,end);
    }
    else{
        //occaTimerTic(mesh->device,"KernelTime");
    tic = MPI_Wtime();  

    occa::streamTag start = mesh->device.tagStream();
    
       // assume 1 mpi process
      for(int it=0;it<iterations;++it){
        testKernels[i](mesh->NinternalElements,
                       zero, // zero offset           
                      mesh->o_vgeo,
                      o_DrsT,
                      mesh->o_DsT,
                      o_X,
                      o_G);
      }


      occa::streamTag end = mesh->device.tagStream();
      mesh->device.finish();  
      toc = MPI_Wtime();
      //      kernelElapsed    = toc-tic;
      kernelElapsed = mesh->device.timeBetween(start,end);

   }

   
      
      const dfloat alpha = 1.0;  // assumption on caching element/face values
      // const dfloat alpha = 0.0;
      // const dfloat alpha = -1.0/9.0 + 1.0/9.0*mesh->N;

      // if(i<=4){ 
      //   Nbytes       = (sizeof(dfloat)*(1*Np +4*alpha*Np + 4*(1.0-alpha)*1 + 4*Np)/2);
      //   NbytesShared = (sizeof(dfloat)*(1*Np + Np*(2*Np))); 
      //   flops        = Np*(Np*4 +6);  // All float ops only
      // }


      if(i==5){
        Nbytes       = (sizeof(dfloat)*(1*Np +4*alpha*Np + 4*(1.0-alpha)*1 + 4*Np)/2);
        NbytesShared = (sizeof(dfloat)*(1*Np + Np*(2*Np) + 4 + 4*Np + Np)); 
        flops        = Np*(Np*4 +6);  // All float ops only
      }

      else{
        Nbytes       = (sizeof(dfloat)*(1*Np +4*alpha*Np + 4*(1.0-alpha)*1 + 4*Np)/2);
        NbytesShared = (sizeof(dfloat)*(1*Np + Np*(2*Np) + 1*Np)); 
        flops        = Np*(Np*4 +6);  // All float ops only
      }



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
      //      double copyElapsed = (toc-tic);
      double copyElapsed = mesh->device.timeBetween(startCopy, endCopy);


      // Compute Data
      double copyBandwidth = mesh->Nelements*((Nbytes*iterations*2)/(1e9*copyElapsed));
      double  bw           = mesh->Nelements*((Nbytes*iterations*2)/(1e9*kernelElapsed));

      double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);
      double d2dbound      = copyBandwidth*gflops/bw;

      double smbound       = 7882*flops/( (double) NbytesShared); // 7751 for 1189 GHz
  
      double intensity    = gflops/bw; 

      double roofline     = mymin(d2dbound, smbound);

      double max_thg_p100 = 4670; 
      double ach_thg      = mymin(549*intensity, max_thg_p100);



      printf("[ N\tK\tDOFS\tKernelTime\tCopyTime\tIntensity\tGFLOPS/s\t(d2dbound GFLOPS/s)\tBW(GB/s)\t(SMBOUND GFLOPS/s)\t(ROOFLINE GFLOPS/s)\tTH_peak]\n");
      printf("%02d \t%02d\t%02d\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\t%6.4E\n",
              mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
       kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline, ach_thg);

      char fname[BUFSIZ];
      sprintf(fname, "KernelData.dat");
      FILE *fp;
      fp = fopen(fname, "a");

      fprintf(fp, "%02d %02d %02d %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E %12.10E\n",
              mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
              kernelElapsed/iterations, copyElapsed/iterations, intensity, gflops, d2dbound, bw, smbound, roofline, ach_thg);
      fclose(fp);

    }

  #endif

 


#if KERNEL_TEST==5

  
  int Nbl      = 20;
  int Nmult    = 10;
  int NKernels = Nbl*Nmult;

  occa::kernel *testKernels = new occa::kernel[NKernels];
  char kernelNames[NKernels][BUFSIZ];

  dfloat mintime = 100.f;

  iint mintblock = 0; 
  iint mintmult  = 0; 

  // for(iint i=6; i<NKernels; i++)
  // {
  
  int i = 0; 
  for (iint b=1;b<=Nbl; b++){

    for(iint m =1; m<=Nmult; m++){

    occa::kernelInfo kernelInfoT  = kernelInfo;
    sprintf(kernelNames[i], "insSubCycleCubatureVolume2D_v5");

    // for kernel 6
    // kernelInfoT.addDefine("p_NbV", b);
    // kernelInfoT.addDefine("p_Nmt", m);
    
    // for kernel 3
    //kernelInfoT.addDefine("p_cubNblockV",b);

     kernelInfoT.addDefine("p_cNbV", b);
     kernelInfoT.addDefine("p_cNmt", m);
   
   
    testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/insSubCycle2D.okl",kernelNames[i], kernelInfoT);
    printf("insSubCycleCubatureVolume Kernel #%02d\n", i);
    printf("Nblock: %d cubNblock: %d N: %d Np: %d cubNp: %d\n", NblockV, cubNblockV, mesh->N, mesh->Np, mesh->cubNp);


    // sync processes
    mesh->device.finish();
    //    MPI_Barrier(MPI_COMM_WORLD);

    //occaTimerTic(mesh->device,"KernelTime");
    tic = MPI_Wtime();  

    occa::streamTag start = mesh->device.tagStream();

      // assume 1 mpi process
      for(int it=0;it<iterations;++it){
        //printf("Cubature Points: %d", mesh->cubNp);
        testKernels[i](mesh->Nelements,
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
      }

      occa::streamTag end = mesh->device.tagStream();
      mesh->device.finish();  
      toc = MPI_Wtime();
      //      kernelElapsed    = toc-tic;
      kernelElapsed = mesh->device.timeBetween(start,end);

    
          
        Nbytes       = (sizeof(dfloat)*(4*Np +4*Np + 2*Np)/2);

        NbytesShared = (sizeof(dfloat)*(4*Nc + 4*Np*Nc + 4*Nc + 4*Np*Nc)); 
        flops        = Np*6 + Np*Nc*8 + 4*Nc + 8*Np*Nc + 2*Np ;  // All float ops only

      
      if(kernelElapsed<mintime){
        mintime = kernelElapsed;
        mintblock = b;
        mintmult  = m; 

      }


      
      double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);

      printf("[ N\tBlock\tNmult\tKernelTime\tGFLOPS/s\t mintime]\n");
      printf("%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed, gflops,  mintime, mintblock, mintmult);
     
     char fname[BUFSIZ];
      sprintf(fname, "KernelOptimization.dat");
      FILE *fp;
      fp = fopen(fname, "a");

      fprintf(fp,"%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed, gflops,  mintime, mintblock, mintmult);
      fclose(fp);

      i++;


     }
    }

    // }

  #endif


#if KERNEL_TEST==6 // Optimize surface kernel

  
  int Nbl      = 10;
  int Nmult    = 5;
  int NKernels = Nbl*Nmult;

  occa::kernel *testKernels = new occa::kernel[NKernels];
  char kernelNames[NKernels][BUFSIZ];

  dfloat mintime = 100.f;

  iint mintblock = 0; 
  iint mintmult  = 0; 
 
  int i = 0; 
  for (iint b=1;b<=Nbl; b++){

    for(iint m =1; m<=Nmult; m++){

    occa::kernelInfo kernelInfoT  = kernelInfo;
    sprintf(kernelNames[i], "insSubCycleCubatureSurface2D_v8");

    
     // kernelInfoT.addDefine("p_cubNblockS", b);
    kernelInfoT.addDefine("p_NnodesS5", m);
    kernelInfoT.addDefine("p_NblockS5", b);
   
   
    testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/optSubcycleSurface2D.okl",kernelNames[i], kernelInfoT);
    printf("insSubCycleCubatureSurface Kernel #%02d\n", i);
    printf("NblockS: %d cubNblockS: %d N: %d Np: %d Nfp: %d cubNfp: %d\n", NblockS, cubNblockS,  mesh->N, mesh->Np, mesh->Nfaces*mesh->Nfp,  mesh->Nfaces*mesh->intNfp );


    // sync processes
    mesh->device.finish();
    // MPI_Barrier(MPI_COMM_WORLD);


           //occaTimerTic(mesh->device,"KernelTime");
    tic = MPI_Wtime();  

    occa::streamTag start = mesh->device.tagStream();

      // assume 1 mpi process
      for(int it=0;it<iterations;++it){
        //printf("Cubature Points: %d", mesh->cubNp);
        testKernels[i](mesh->Nelements,
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
      }

      occa::streamTag end = mesh->device.tagStream();
      mesh->device.finish();  
      toc = MPI_Wtime();

      //      kernelElapsed    = toc-tic;
      kernelElapsed = mesh->device.timeBetween(start,end);
    




  
       Nbytes           = (sizeof(dfloat)*(8*Ntfp + 4*intNtfp + 4*mesh->Np) + sizeof(iint)*(2*Ntfp))/2;

       NbytesShared     = (sizeof(dfloat)*(8*Ntfp + 8*intNtfp*Nfp + 2*intNtfp + 2*Np*intNtfp)); 

       flops            = intNtfp*( Nfp*16 + 6 + 1 + 28) + Np*intNtfp*4;

    
      
      if(kernelElapsed<mintime){
        mintime = kernelElapsed;
        mintblock = b;
        mintmult  = m; 

      }


      
      double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);

      printf("[ N\tBlock\tNmult\tKernelTime\tGFLOPS/s\t mintime]\n");
      printf("%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed/iterations, gflops,  mintime, mintblock, mintmult);
     
     char fname[BUFSIZ];
      sprintf(fname, "KernelOptimization.dat");
      FILE *fp;
      fp = fopen(fname, "a");

      fprintf(fp,"%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed/iterations, gflops,  mintime, mintblock, mintmult);
      fclose(fp);

      i++;


     }
    }

    // }

  #endif




  #if KERNEL_TEST==7 // Optimize surface kernel

  
  int Nbl      = 10;
  int Nmult    = 1;
  int NKernels = Nbl*Nmult;

  occa::kernel *testKernels = new occa::kernel[NKernels];
  char kernelNames[NKernels][BUFSIZ];

  dfloat mintime = 100.f;

  iint mintblock = 0; 
  iint mintmult  = 0; 


  int i = 0; 
  for (iint b=1;b<=Nbl; b++){

    for(iint m =1; m<=Nmult; m++){

    occa::kernelInfo kernelInfoT  = kernelInfo;
    sprintf(kernelNames[i], "ellipticPartialAxIpdgTri2D_v4");

     //   kernelInfoT.addDefine("p_Nblock", b);
    kernelInfoT.addDefine("p_Nmt", m);
    kernelInfoT.addDefine("p_NbV", b);
   
   
    testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/optAxIpdgTri2D.okl",kernelNames[i], kernelInfoT);
    printf("Ax Kernel #%02d\n", i);
    printf("N: %d Np: %d Nfp: %d\n", mesh->N, mesh->Np, mesh->Nfaces*mesh->Nfp);


    // sync processes
    mesh->device.finish();
    //    MPI_Barrier(MPI_COMM_WORLD);

    //occaTimerTic(mesh->device,"KernelTime");
    tic = MPI_Wtime();  

    occa::streamTag start = mesh->device.tagStream();

      // assume 1 mpi process
      for(int it=0;it<iterations;++it){
        testKernels[i](mesh->NinternalElements,
                                  mesh->o_internalElementIds,
                                  mesh->o_vmapM,
                                  mesh->o_vmapP,
                                  lambda, // lamda = 0.0
                                  tau,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  mesh->o_EToB,
                                  mesh->o_DrT,
                                  mesh->o_DsT,
                                  mesh->o_LIFTT,
                                  mesh->o_MM,
                                  o_G,
                                  o_X);
      
      }

      occa::streamTag end = mesh->device.tagStream();
      mesh->device.finish();  
      toc = MPI_Wtime();
      //      kernelElapsed    = toc-tic;
      kernelElapsed = mesh->device.timeBetween(start,end);


      const dfloat alpha = 1.0;  // assumption on caching element/face values
      // const dfloat alpha = 0.0;
      // const dfloat alpha = -1.0/9.0 + 1.0/9.0*mesh->N;

      Nbytes = (sizeof(dfloat)*(4*Np    // float4
                               +8*Ntfp // float4
                               +alpha*5*Nfaces + (1-alpha)*5*Ntfp   // alpha = 1 means cached face geometric factors
                               +alpha*5 +(1-alpha)*5*Np   // alpha 1 means drdx, J etc.. are all cached for an element
                               +1*mesh->Np))/2           // gloabl write of output
              + (sizeof(iint)*(alpha*1 + (1-alpha)*Np*1  // element list
                              +2*Ntfp  // idM, idP
                              +alpha*Nfaces + (1-alpha)*Ntfp // BC's
                              ))/2;  

       NbytesShared = (sizeof(dfloat)*(3*Np + 3*Ntfp + Np*(2*Ntfp +2+4) + Ntfp*3 + Np*(2*Np + 1) + Np*(1*Ntfp + 1) + Np*(1*Np))); 

       flops = 1*Np + 21*Ntfp + Np*(4*Ntfp + 8) + 6*Ntfp+ Np*(4*Np + 2)+ Np*(1+2*Ntfp) + Np*(1+2*Np); 

      
      if(kernelElapsed<mintime){
        mintime = kernelElapsed;
        mintblock = b;
        mintmult  = m; 

      }


      
      double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);

      printf("[ N\tBlock\tNmult\tKernelTime\tGFLOPS/s\t mintime]\n");
      printf("%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed/iterations, gflops,  mintime, mintblock, mintmult);
     
     char fname[BUFSIZ];
      sprintf(fname, "KernelOptimization.dat");
      FILE *fp;
      fp = fopen(fname, "a");

      fprintf(fp,"%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed/iterations, gflops,  mintime, mintblock, mintmult);
      fclose(fp);

      i++;


     }
    }

    // }

  #endif

  
 #if KERNEL_TEST==8 // Optimize surface kernel

  
  int Nbl      = 1;
  int Nmult    = 1;
  int NKernels = Nbl*Nmult;

  occa::kernel *testKernels = new occa::kernel[NKernels];
  char kernelNames[NKernels][BUFSIZ];

  dfloat mintime = 100.f;

  iint mintblock = 0; 
  iint mintmult  = 0; 


  int i = 0; 
  for (iint b=1;b<=Nbl; b++){

    for(iint m =1; m<=Nmult; m++){

    occa::kernelInfo kernelInfoT  = kernelInfo;
    sprintf(kernelNames[i], "ellipticPartialGradientTri2D_v5");

     //   kernelInfoT.addDefine("p_Nblock", b);
    kernelInfoT.addDefine("p_Nmt2", m);
    kernelInfoT.addDefine("p_NbV2", b);
   
   
    testKernels[i] = mesh->device.buildKernelFromSource(DHOLMES "/okl/optGradientTri2D.okl",kernelNames[i], kernelInfoT);
    printf("Gradient Kernel #%02d\n", i);
    printf("N: %d Np: %d Nfp: %d\n", mesh->N, mesh->Np, mesh->Nfaces*mesh->Nfp);


    // sync processes
    mesh->device.finish();
    //    MPI_Barrier(MPI_COMM_WORLD);

    //occaTimerTic(mesh->device,"KernelTime");
    tic = MPI_Wtime();  

    occa::streamTag start = mesh->device.tagStream();


     #if 1
      // assume 1 mpi process
      for(int it=0;it<iterations;++it){
        testKernels[i](mesh->NinternalElements,
                       zero, // zero offset           
                      mesh->o_vgeo,
                      mesh->o_DrT,
                      mesh->o_DsT,
                      o_X,
                      o_G);

      }
    
     #else
       // assume 1 mpi process
      for(int it=0;it<iterations;++it){
        testKernels[i](mesh->NinternalElements,
                       zero, // zero offset           
                      mesh->o_vgeo,
                      o_DrsT,
                      mesh->o_DsT,
                      o_X,
                      o_G);
      }

      #endif





      occa::streamTag end = mesh->device.tagStream();
      mesh->device.finish();  
      toc = MPI_Wtime();
      //      kernelElapsed    = toc-tic;
      kernelElapsed = mesh->device.timeBetween(start,end);


      const dfloat alpha = 1.0;  // assumption on caching element/face values
      // const dfloat alpha = 0.0;
      // const dfloat alpha = -1.0/9.0 + 1.0/9.0*mesh->N;

      
      Nbytes       = (sizeof(dfloat)*(1*Np +4*alpha*Np + 4*(1.0-alpha)*1 + 4*Np)/2);
      NbytesShared = (sizeof(dfloat)*(1*Np + Np*(2*Np))); 
      flops        = Np*(Np*4 +6);  // All float ops only
      
      if(kernelElapsed<mintime){
        mintime = kernelElapsed;
        mintblock = b;
        mintmult  = m; 

      }


      
      double gflops        = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);

      printf("[ N\tBlock\tNmult\tKernelTime\tGFLOPS/s\t mintime]\n");
      printf("%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed/iterations, gflops,  mintime, mintblock, mintmult);
     
      char fname[BUFSIZ];
      sprintf(fname, "KernelOptimization.dat");
      FILE *fp;
      fp = fopen(fname, "a");

      fprintf(fp,"%02d %02d %02d %12.10E %12.10E %12.10E %02d %02d\n", mesh->N, b, m, kernelElapsed/iterations, gflops,  mintime, mintblock, mintmult);
      fclose(fp);

      i++;


     }
    }

    // }

  #endif

   
  
  
  


  



 
  #if 0

   #if KERNEL_TEST==2

  TestKernel = mesh->device.buildKernelFromSource(DHOLMES 
                            "/okl/insSubCycle2D.okl","insSubCycleCubatureSurface2D"
                            ,kernelInfo);
  iint KernelId = 0; // Cubature Integration


  // TestKernel = mesh->device.buildKernelFromSource(DHOLMES 
  //                           "/okl/insAdvection2D.okl","insAdvectionSurface2D"
  //                           ,kernelInfo);
  // iint KernelId = 1; // Cubature Integration

  printf("Nblock: %d cubNblock: %d N: %d Nfp: %d cubNfp: %d\n", NblockS, cubNblockS, mesh->N, mesh->Nfp, mesh->intNfp);
  
  #endif


 
 
 double memory_used = mesh->device.memoryAllocated()/(1024.*1024.*1024);
 printf("Memory Allocated (GB) = %g\n", memory_used);


 
  

  double tic = 0.0, toc = 0.0, kernelElapsed=0.0;

 if(KernelId==0){
    // sync processes
    mesh->device.finish();
    MPI_Barrier(MPI_COMM_WORLD);

    //occaTimerTic(mesh->device,"KernelTime");
    tic = MPI_Wtime();  
    occa::streamTag start = mesh->device.tagStream();

    // assume 1 mpi process
    for(int it=0;it<iterations;++it){
      #if KERNEL_TEST==1
      //printf("Cubature Points: %d", mesh->cubNp);
      TestKernel(mesh->Nelements,
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
       TestKernel(mesh->Nelements,
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

    }

    occa::streamTag end = mesh->device.tagStream();
    mesh->device.finish();  
    toc = MPI_Wtime();
    kernelElapsed    = toc-tic;
    kernelElapsed    = mesh->device.timeBetween(start, end);
}
if(KernelId==1){
    // sync processes
    mesh->device.finish();
    MPI_Barrier(MPI_COMM_WORLD);

    //occaTimerTic(mesh->device,"KernelTime");
    tic = MPI_Wtime();  
    occa::streamTag start = mesh->device.tagStream();

    // assume 1 mpi process
    for(int it=0;it<iterations;++it){
      #if KERNEL_TEST==1
        TestKernel(mesh->Nelements,
                  mesh->o_vgeo,
                  mesh->o_DrT,
                  mesh->o_DsT,
                  zero,
                  o_U,
                  o_V,
                  o_X,
                  o_Y);

      #endif


      #if KERNEL_TEST==2
       TestKernel(mesh->Nelements,
        mesh->o_sgeo,
        mesh->o_LIFTT,
        mesh->o_vmapM,
        mesh->o_vmapP,
        mesh->o_EToB,
        time,
        mesh->o_x,
        mesh->o_y,
        index,
        o_U,
        o_V,
        o_X,
        o_Y);
      #endif

    }

    occa::streamTag end = mesh->device.tagStream();
    mesh->device.finish();  
    toc = MPI_Wtime();
    kernelElapsed    = toc-tic;
}

  
 
  

  #if KERNEL_TEST==1
    if(KernelId==0){ // Cubature Integration
     flops = Nc*Np*8 + 4*Nc + Np*Nc*16 + Np*14 + Np*2;  // All float ops only
     Nbytes =(sizeof(dfloat)*(6*mesh->Np +4*mesh->Np)/2);

    }
    else{
    flops =  3*Np + Np*Np*12 + Np*14;
    Nbytes =(sizeof(dfloat)*(2*mesh->Np  + 4*mesh->Np+ 2*mesh->Np)/2);

    }
  #endif

  #if KERNEL_TEST==2
    if(KernelId==0){ 
     flops = intNtfp*( Nfp*16 + 6 + 1 + 28) + Np*intNtfp*4;  // All float ops only  
      //Nbytes = (sizeof(dfloat)*(8*Ntfp + 4*intNtfp + 4*mesh->Np))/2;  
      Nbytes = (sizeof(dfloat)*(8*Ntfp + 4*mesh->Nfaces + 4*mesh->Np))/2;  //?????  
     
    }

    if(KernelId==1){ 
     flops = Ntfp*(6 + 28 + 1) + Np*Ntfp*4;  // All float ops only  
      Nbytes = (sizeof(dfloat)*(8*Ntfp + 4*mesh->Np))/2;   
     
    }

  #endif
  
   #if KERNEL_TEST==3
   flops = 1*Np + 21*Ntfp + Np*(4*Ntfp + 8) +6*Ntfp+ Np*(4*Np + 2)+ Np*(1+2*Ntfp) + Np*(1+2*Np);  // All float ops only      
  #endif

  #if KERNEL_TEST==4
   flops = Np*(Np*2 + Np*2 + 6);  // All float ops only      
  
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
   
  
  // Compute Data
 double copyBandwidth = mesh->Nelements*((Nbytes*iterations*2)/(1e9*copyElapsed));
 double  bw           = mesh->Nelements*((Nbytes*iterations*2)/(1e9*kernelElapsed));

  
  double gflops   = mesh->Nelements*flops*iterations/(1e9*kernelElapsed);
  double d2dbound = copyBandwidth*gflops/bw;


  // Intensity = flops/bytes
  double intensity = gflops/bw; 

  double max_thg_p100 = 4700; 
  double ach_thg     = mymin(549*intensity, max_thg_p100);

  printf("[ N\tK\tDOFS\tKernelTime\tCopyTime\tIntensity\tGFLOPS/s\t(d2dbound GFLOPS/s)\tBW(GB/s)\tTH_peak\t(GFLOPS/s)]\n");
  printf("%02d %02d\t%02d\t%17.15E\t%17.15E\t%17.15E\t%17.15E\t%17.15E\t%17.15E\t%17.15E\n",
              mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
              kernelElapsed, copyElapsed, intensity, gflops, d2dbound, bw, ach_thg);

    char fname[BUFSIZ];
    sprintf(fname, "KernelData.dat");
    FILE *fp;
    fp = fopen(fname, "a");

    fprintf(fp,"%02d %02d\t%02d\t%17.15E\t%17.15E\t%17.15E\t%17.15E\t%17.15E\t%17.15E\t%17.15E\n",
              mesh->N, mesh->Nelements,(mesh->Nelements*mesh->Np), 
              kernelElapsed, copyElapsed, intensity, gflops, d2dbound, bw, ach_thg);
    fclose(fp);
  // }

#endif
  




#if 0
  dfloat time = 0.0; 
  iint iterations = 10;
  // sync processes
  mesh->device.finish();
  MPI_Barrier(MPI_COMM_WORLD);

 //occaTimerTic(mesh->device,"KernelTime");
  double tic = MPI_Wtime();  
  occa::streamTag start = mesh->device.tagStream();

  // assume 1 mpi process
  for(int it=0;it<iterations;++it){
  
  #if KERNEL_TEST==1
    //printf("Cubature Points: %d", mesh->cubNp);
    testKernel1(mesh->Nelements,
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
  
  
  double copyBandwidth = mesh->Nelements*((Nbytes*iterations*2)/(1e9*copyElapsed));

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
  
  double gflops   = globalElements*flops*iterations/(1e9*globalElapsed);
  double d2dbound = copyBandwidth*flops/bw;
  double bwth     = 549*flops/bw;


  bw *= (globalElements*iterations)/(1024.*1024*1024*globalElapsed);



  if(rank==root){
    printf("[ RANKS N DOFS ELAPSEDTIME ITERATIONS (DOFS/RANKS) (DOFS/TIME/ITERATIONS/RANKS) (Kernel GFLOPS) (d2dbound GB/s) (copy GB/s) (achieved GP/s)]\n");
    printf("%02d %02d %02d %17.15lg %d %17.15E %17.15E %17.15E %17.15E %17.15E %17.15E\t\n",
           size, mesh->N, globalDofs, globalElapsed, iterations, globalDofs/(double)size,
           (globalDofs*iterations)/(globalElapsed*size), gflops, d2dbound, copyBandwidth, bw);

     printf("%02d\t%02d\t%17.15lg\t%17.15E\t%17.15E\n",
              mesh->N, globalDofs, globalElapsed, gflops, d2dbound);

    char fname[BUFSIZ];
    sprintf(fname, "KernelData.dat");
    FILE *fp;
    fp = fopen(fname, "a");

    fprintf(fp,"%02d\t%02d\t%17.15E\t%17.15E\t%17.15E\t%17.15E\t%17.15E\n",
              mesh->N, globalDofs, globalElapsed, gflops, d2dbound, bwth, bw);
    fclose(fp);




  }
  

  
 #endif 



  }



