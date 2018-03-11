#include "ellipticBenchmarkTet3D.h"

void ellipticRunBenchmark3D(solver_t *solver, char *options, occa::kernelInfo kernelInfo, char *kernelFileName, int Nblocks, int Nnodes){

  mesh3D *mesh = solver->mesh;

  int Ntrials = 10;

  size_t L1CacheSize = 24576; //L1 cache size of test device in bytes (24KB for P100)

  int NKernels;
  char kernelName[BUFSIZ];

  NKernels = 2;
  sprintf(kernelName, "ellipticPartialAxSparseTet3D");

  //  kernelInfo.addCompilerFlag("-G");

  dfloat time = 0.;
  printf("test 0 \n");

  char testkernelName[BUFSIZ];
  occa::kernel testKernel;

  int Nbytes = sizeof(dfloat)*mesh->Np*2; // load one field, save one filed (ignore geofacs)
  //add the matrices (?)
  //Nbytes += sizeof(dfloat)*(test[mesh->Np]-1)*3;
  Nbytes /= 2;
  occa::memory o_foo = mesh->device.malloc(Nbytes*mesh->Nelements);
  occa::memory o_bah = mesh->device.malloc(Nbytes*mesh->Nelements);
  mesh->device.finish();

  occa::streamTag startCopy = mesh->device.tagStream();
  for(int it=0;it<Ntrials;++it){
    o_bah.copyTo(o_foo);
  }
  occa::streamTag endCopy = mesh->device.tagStream();
  mesh->device.finish();

  double  copyElapsed = mesh->device.timeBetween(startCopy, endCopy);
  printf("copied \n");


  for(int i=0; i<NKernels; i++) {

    sprintf(testkernelName, "%s_v%d", kernelName,  i);
    printf("%s================= Kernel #%02d================================================\n\n", testkernelName, i);
    printf("Np = %d sizeof(dfloat) = %d Nelements = %d \n", mesh->Np, sizeof(dfloat), mesh->Nelements);

    testKernel = mesh->device.buildKernelFromSource(kernelFileName,testkernelName,kernelInfo);

    dfloat lambda = 0;
    // sync processes


    // time Ax
    double timeAx = 0.0f;
    double kernelElapsed =0.0f;;

    mesh->device.finish();

    occa::streamTag start[Ntrials], end[Ntrials];
    for(int it=0;it<Ntrials;++it){
      start[it] = mesh->device.tagStream();

#if 1
      // tet sparse kernel
      testKernel(mesh->Nelements, 
		 solver->o_localGatherElementList,
		 mesh->o_ggeo, 
		 mesh->o_IndTchar,
		 mesh->o_SrrT, 
		 mesh->o_SrsT, 
		 mesh->o_SrtT, 
		 mesh->o_SssT,		 
		 mesh->o_SstT,		 
		 mesh->o_SttT,		 
		 mesh->o_MM, 
		 lambda,
		 solver->o_p,
		 solver->o_grad);
#else
      testKernel(mesh->Nelements, 
		 solver->o_localGatherElementList,
		 mesh->o_ggeo, 
		 mesh->o_rowData,          
		 mesh->o_Srr, 
		 mesh->o_Srs, 
		 mesh->o_Srt, 
		 mesh->o_Sss,              
		 mesh->o_Sst,              
		 mesh->o_Stt,              
		 mesh->o_MM, 
		 lambda,
		 solver->o_p,
		 solver->o_grad);
#endif

      end[it] = mesh->device.tagStream();
    }
    mesh->device.finish();
    for(int it=0;it<Ntrials;++it){
      timeAx = mesh->device.timeBetween(start[it],end[it]);
      kernelElapsed +=timeAx;    
    }
    mesh->device.finish();
    kernelElapsed = kernelElapsed/(dfloat)Ntrials;

    //start
    //
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int   localDofs = mesh->Np*mesh->Nelements;
    int   localElements = mesh->Nelements;
    int   globalDofs;
    int   globalElements;
    double globalCopyElapsed;

    MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_INT,   MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce(&localElements,&globalElements,1, MPI_INT,   MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce(&copyElapsed,&globalCopyElapsed,1, MPI_DOUBLE,   MPI_MAX, 0, MPI_COMM_WORLD );


    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank==0){
      printf("Ntrials = %d global copy el %f local copy el %f size(dfloat) %d \n",Ntrials, globalCopyElapsed, copyElapsed, sizeof(dfloat));

      // count actual number of non-zeros

      int nnzs = 0;
      for(int n=0;n<mesh->maxNnzPerRow*mesh->Np;++n)
	nnzs += (mesh->Ind[n]!=0);
      
      printf("\n nnzs = %d matrix size %d padded row %d \n", nnzs, mesh->Np*mesh->maxNnzPerRow, mesh->maxNnzPerRow);


      // 6 flops per non-zero plus chain rule
      double flops = (double)(nnzs*12 + mesh->Np*11);
      double eqflops = mesh->Np*mesh->Np*12 + mesh->Np*11;

      double roofline = ((mesh->Nelements*flops*(double)Ntrials))/(1e9*globalCopyElapsed);
      printf("Nelements = %d flops = %d Ntrials = %d copy elapsed scaled = %f kernelElapsed %16.16f\n",mesh->Nelements, (int) flops, Ntrials,  1e9*globalCopyElapsed, kernelElapsed);

      double copyBandwidth   = mesh->Nelements*((Nbytes*Ntrials*2.)/(1e9*globalCopyElapsed));
      double kernelBandwidth = mesh->Nelements*((Nbytes*2.)/(1e9*kernelElapsed));
      double kernelGFLOPS = mesh->Nelements*flops/(1e9*kernelElapsed);

      double kernelEquivGFLOPS = mesh->Nelements*eqflops/(1e9*kernelElapsed);

      printf("Copy Bw %16.17g achieved BW %16.17g\n", copyBandwidth, kernelBandwidth);
      printf("ROOFLINE %16.17g \n", roofline);
      printf("GFLOPS %16.17f (%16.17f) \n", kernelGFLOPS, kernelEquivGFLOPS);
      printf("time per kernel %f \n",kernelElapsed);
      printf("PARAMETERS %d %d %16.17f \n ", Nblocks, Nnodes, kernelGFLOPS );

    }
  }
}
