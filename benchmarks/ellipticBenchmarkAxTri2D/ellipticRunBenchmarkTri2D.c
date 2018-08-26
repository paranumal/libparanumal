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

#include "ellipticBenchmarkTri2D.h"

void ellipticRunBenchmark2D(solver_t *solver, char *options, occa::kernelInfo kernelInfo, char *kernelFileName, int Nblocks, int Nnodes){

  mesh2D *mesh = solver->mesh;

  int Ntrials = 10;

  size_t L1CacheSize = 24576; //L1 cache size of test device in bytes (24KB for P100)

  int NKernels;
  char kernelName[BUFSIZ];

  NKernels = 4;
  sprintf(kernelName, "ellipticAxNEWTri2D");

  //  kernelInfo.addCompilerFlag("-G");

  dfloat time = 0.;
  printf("test 0 \n");

  /*int  * test = (int *) calloc(mesh->Np+1, sizeof(int));
    printf("test 1 \n");
    mesh->o_India.copyTo(test);
    printf("test 2\n");

    printf("test 1 \n");
    printf("\n");
    for (int nn=0; nn<mesh->Np+1; ++nn){
    printf(" %d ",  test[nn]);
    }
    */

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



  for(int i=1; i<NKernels; i++) {

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
#if 0
      testKernel(mesh->Nelements, 
          solver->o_localGatherElementList,
          mesh->o_ggeo, 
          mesh->o_sparseStackedNZ, 
          mesh->o_sparseSrrT, 
          mesh->o_sparseSrsT, 
          mesh->o_sparseSssT,		 
          mesh->o_MM, 
          lambda,
          solver->o_p,
          solver->o_grad);

      testKernel(mesh->Nelements, 
          solver->o_localGatherElementList,
          mesh->o_ggeo, 
          mesh->o_India,
          mesh->o_Indja, 
          mesh->o_Srr, 
          mesh->o_Srs, 
          mesh->o_Sss,              
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
          mesh->o_Sss,              
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
#if 1
      int nnzs = 0;
      printf("sparse nnz per row: %d \n", mesh->SparseNnzPerRow);
      for(int n=0;n<mesh->Np*mesh->SparseNnzPerRow;++n){
        //        nnzs += (fabs(mesh->sparseSrrT[n])>1e-13);
        //      nnzs += (fabs(mesh->sparseSrsT[n])>1e-13);
        //    nnzs += (fabs(mesh->sparseSssT[n])>1e-13);
        //printf("%16.16f %16.16f %16.16f \n", mesh->sparseSrrT[n],  mesh->sparseSrsT[n],  mesh->sparseSssT[n]);
        nnzs += (mesh->sparseStackedNZ[n]>0); 
      }
      //nnzs*=3;
#else
      nnzs = test[mesh->Np]-1;
#endif
      printf("\n nnzs = %d matrix size %d padded row %d \n", nnzs, mesh->Np*mesh->SparseNnzPerRow, mesh->SparseNnzPerRow);


      // 6 flops per non-zero plus chain rule
      double flops = (double)(nnzs*6 + mesh->Np*5);
      double eqflops = mesh->Np*mesh->Np*6 + mesh->Np*5;

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
