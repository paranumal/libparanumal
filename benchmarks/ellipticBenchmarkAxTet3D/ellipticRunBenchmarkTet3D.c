/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "ellipticBenchmarkTet3D.h"

//#define BB_TESTS 1

void ellipticRunBenchmark3D(solver_t *solver, char *options, occa::kernelInfo kernelInfo, char *kernelFileName, int Nblocks, int Nnodes){

  mesh3D *mesh = solver->mesh;

  int Ntrials = 20;

  size_t L1CacheSize = 24576; //L1 cache size of test device in bytes (24KB for P100)

  int NKernels;
  char kernelName[BUFSIZ];

#ifndef BB_TESTS
  NKernels = 1;
  sprintf(kernelName, "ellipticPartialAxTet3D");
#else
  NKernels = 6;
  sprintf(kernelName, "ellipticPartialAxBBTet3D");
#endif
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

  // build faceNodes array on device
  occa::memory o_faceNodes = mesh->device.malloc(mesh->Nfp*mesh->Nfaces*sizeof(int), mesh->faceNodes);

  for(int i=0; i<NKernels; i++) {

    sprintf(testkernelName, "%s_v%d", kernelName,  i);
    
    printf("%s================= Kernel #%02d================================================\n\n", testkernelName, i);
    printf("Np = %d sizeof(dfloat) = %d Nelements = %d \n", mesh->Np, sizeof(dfloat), mesh->Nelements);

    //    kernelInfo.addDefine("p_Ne", Nnodes);
    //    kernelInfo.addDefine("p_Nb", Nblocks);

    int halfNp = (mesh->Np+1)/2;
    kernelInfo.addDefine("p_halfNp", halfNp);
    
    testKernel = mesh->device.buildKernelFromSource(kernelFileName,testkernelName,kernelInfo);

    dfloat lambda = 0;
    // sync processes

    // time Ax
    double timeAx = 0.0f;
    double kernelElapsed =0.0f;;

    mesh->device.finish();

    occa::streamTag start[Ntrials], end[Ntrials];
    start[0] = mesh->device.tagStream();
    for(int it=0;it<Ntrials+1;++it){

      if(it==1) // tag after warm up
	start[0] = mesh->device.tagStream();

#ifndef BB_TESTS
      // tet sparse kernel
      testKernel(mesh->Nelements, 
		 solver->o_localGatherElementList,
		 mesh->o_ggeo, 
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
      // v0: basic version
      int elementOffset = 0;

      testKernel(mesh->Nelements,
		 elementOffset,
		 mesh->o_sgeo,   // surface geofacs
		 o_faceNodes, // volume node indices of surface nodes
		 mesh->o_vgeo,   // volume geofacs (drdx, ...)
		 mesh->o_D0ids,     // sparse D matrices
		 mesh->o_D1ids,
		 mesh->o_D2ids,
		 mesh->o_D3ids,
		 mesh->o_Dvals,    
		 mesh->o_LIFTT,   // dense lift matrix
		 mesh->o_MM,     // dense mass matrix
		 lambda,
		 solver->o_p,
		 solver->o_grad);
#endif

    }
    end[0] = mesh->device.tagStream();
	  
    mesh->device.finish();
    kernelElapsed = mesh->device.timeBetween(start[0],end[0]);

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
#ifndef BB_TESTS
      double flops = (double)(nnzs*14 + mesh->Np*13);
#else
      double flops = (double)(mesh->Np*(39+15+1) +
			      mesh->Nfp*mesh->Nfaces*(6 + 2*mesh->Np) +
			      mesh->Np*15 +
			      (39+3)*mesh->Np + 2.*mesh->Np*mesh->Np);
#endif
      double eqflops = mesh->Np*mesh->Np*14 + mesh->Np*13;

      double roofline = ((mesh->Nelements*flops*(double)Ntrials))/(1e9*globalCopyElapsed);
      printf("Nelements = %d flops = %d Ntrials = %d "
	     "copy elapsed scaled = %f kernelElapsed %16.16f\n",
	     mesh->Nelements, (int) flops, Ntrials,  1e9*globalCopyElapsed, kernelElapsed);

      double copyBandwidth   = mesh->Nelements*((Nbytes*Ntrials*2.)/(1e9*globalCopyElapsed));
      double kernelBandwidth = mesh->Nelements*((Nbytes*2.)/(1e9*kernelElapsed));

      double kernelGFLOPS = mesh->Nelements*flops/(1e9*kernelElapsed);

      double kernelEquivGFLOPS = mesh->Nelements*eqflops/(1e9*kernelElapsed);

      printf("Copy Bw %16.17g achieved BW %16.17g\n", copyBandwidth, kernelBandwidth);
      printf("ROOFLINE %16.17g \n", roofline);
      printf("GFLOPS %16.17f (%16.17f) \n", kernelGFLOPS, kernelEquivGFLOPS);
      printf("time per kernel %f \n",kernelElapsed);
      printf("PARAMETERS %d %d %d %d %16.17f %16.17f %16.17f\n ",
	     mesh->N, i, Nblocks, Nnodes, kernelGFLOPS, kernelEquivGFLOPS, kernelGFLOPS/kernelEquivGFLOPS );

    }
  }
}
