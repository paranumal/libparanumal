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

#include "elliptic.hpp"
#include "mesh/meshDefines3D.h"

#include <algorithm> 
#define nnz_t parAlmond::parCOO::nonZero_t

// compare on global indices
bool parallelCompareRowColumnV2(nnz_t &a, nnz_t &b){
  if(a.row < b.row) return +1;
  if(a.row > b.row) return  0;

  if(a.col < b.col) return +1;
  if(a.col > b.col) return  0;

  return 0;
}


int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  if(argc!=2)
    LIBP_ABORT(string("Usage: ./ellipticMain setupfile"));

  //create default settings
  platformSettings_t platformSettings(comm);
  meshSettings_t meshSettings(comm);
  ellipticSettings_t ellipticSettings(comm);
  ellipticAddRunSettings(ellipticSettings);

  //load settings from file
  ellipticSettings.parseFromFile(platformSettings, meshSettings,
                                 argv[1]);

  // set up platform
  platform_t platform(platformSettings);

  platformSettings.report();
  meshSettings.report();
  ellipticSettings.report();

  // set up mesh
  mesh_t& mesh = mesh_t::Setup(platform, meshSettings, comm);

  dfloat lambda = 0.0;
  ellipticSettings.getSetting("LAMBDA", lambda);

  // Boundary Type translation. Just defaults.
  int NBCTypes = 3;
  int BCType[NBCTypes] = {0,1,2};

  // set up elliptic solver
  elliptic_t& elliptic = elliptic_t::Setup(platform, mesh, ellipticSettings,
                                           lambda, NBCTypes, BCType);

  // run
  elliptic.Run();

#if 1
  int testHOST = 0;
  
  printf("Building Matrix\n");
  
  parAlmond::parCOO A(elliptic.platform, mesh.comm);

  if(testHOST)
    elliptic.BuildOperatorMatrixContinuous(A);

  occa::properties kernelInfo = mesh.props;

  if(sizeof(hlong)==8)
    kernelInfo["defines/hlong"]= "long long int";
  if(sizeof(hlong)==4)
    kernelInfo["defines/hlong"]= "int";

  kernelInfo["defines/" "p_G00ID"]= G00ID;
  kernelInfo["defines/" "p_G01ID"]= G01ID;
  kernelInfo["defines/" "p_G02ID"]= G02ID;
  kernelInfo["defines/" "p_G11ID"]= G11ID;
  kernelInfo["defines/" "p_G12ID"]= G12ID;
  kernelInfo["defines/" "p_G22ID"]= G22ID;
  kernelInfo["defines/" "p_GWJID"]= GWJID;

  char kernelName[BUFSIZ];
  switch(mesh.elementType){
  case TRIANGLES:
    sprintf(kernelName, "ellipticBuildOperatorMatrixContinuousTri2D");
    break;
  case QUADRILATERALS:
    sprintf(kernelName, "ellipticBuildOperatorMatrixContinuousQuad2D");
    break;
  case TETRAHEDRA:
    sprintf(kernelName, "ellipticBuildOperatorMatrixContinuousTet3D");
    break;
  case HEXAHEDRA:
    sprintf(kernelName, "ellipticBuildOperatorMatrixContinuousHex3D");
    break;
  }

  occa::kernel buildMatrixKernel =
    elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticBuildOperatorMatrixContinuous.okl",
				  kernelName,
				  kernelInfo);
  
  occa::memory o_maskedGlobalNumbering =
    platform.malloc(mesh.Np*mesh.Nelements*sizeof(hlong), elliptic.maskedGlobalNumbering);

  dlong Nnz = mesh.Nelements*mesh.Np*mesh.Np;
  // note - have to zero matrix before building because of unwritten boundary nodes
  nnz_t *d_AL = (nnz_t*) calloc(Nnz,sizeof(nnz_t));
  nnz_t *h_AL = (nnz_t*) calloc(Nnz,sizeof(nnz_t));
  nnz_t *h_AL2 = (nnz_t*) calloc(Nnz,sizeof(nnz_t));
  
  occa::memory o_AL =  platform.malloc(Nnz*sizeof(nnz_t), d_AL);

  platform.device.finish();
  double t0 = MPI_Wtime();
  
  switch(mesh.elementType){
  case TRIANGLES:
    buildMatrixKernel(mesh.Nelements, o_maskedGlobalNumbering,
		      mesh.o_S, mesh.o_MM, mesh.o_ggeo,
		      elliptic.lambda, o_AL);

    break;
  case QUADRILATERALS:
    buildMatrixKernel(mesh.Nelements, o_maskedGlobalNumbering,
		      mesh.o_D, mesh.o_ggeo,
		      elliptic.lambda, o_AL);

    break;
  case TETRAHEDRA:
    buildMatrixKernel(mesh.Nelements, o_maskedGlobalNumbering,
		      mesh.o_S, mesh.o_MM, mesh.o_ggeo, elliptic.lambda, o_AL);
    break;
  case HEXAHEDRA:
    buildMatrixKernel(mesh.Nelements, o_maskedGlobalNumbering,
		      mesh.o_D,  mesh.o_ggeo,
		      elliptic.lambda, o_AL);
    break;
  }

  platform.device.finish();
  double t1 = MPI_Wtime();

  if(testHOST){
    switch(mesh.elementType){
    case TRIANGLES:
      elliptic.BuildOperatorMatrixContinuousTri2D(h_AL); 
      break;
    case QUADRILATERALS:
      elliptic.BuildOperatorMatrixContinuousQuad2D(h_AL); 
      break;
    case TETRAHEDRA:
      elliptic.BuildOperatorMatrixContinuousTet3D(h_AL); 
      break;
    case HEXAHEDRA:
      elliptic.BuildOperatorMatrixContinuousHex3D(h_AL); 
      break;
    }
  }

  double t2 = MPI_Wtime();

  memcpy(h_AL2, h_AL, Nnz*sizeof(nnz_t));
  
  printf("DEVICE build took: %e\n", t1-t0);
  printf("HOST   build took: %e\n", t2-t1);

  o_AL.copyTo(d_AL);

#if 0
  // double check results
  double tol = 1e-10;
  for(int n=0;n<Nnz;++n){
    nnz_t tmp1 = h_AL[n], tmp2 = d_AL[n];
    if(tmp1.row != tmp2.row ||
       tmp1.col != tmp2.col ||
       fabs(tmp1.val-tmp2.val)>tol){

      printf("mismatch: %d,  (" hlongFormat "," hlongFormat ",%e) => (" hlongFormat"," hlongFormat ",%e)\n", 
	     n,
	     tmp1.row, tmp1.col, tmp1.val,
	     tmp2.row, tmp2.col, tmp2.val);
    }
  }
#endif

  deviceScan_t scanner(platform, DELLIPTIC "okl/nonZero.h", DELLIPTIC "okl/nonZeroCompare2.h", kernelInfo);
  deviceSort_t  sorter(platform, DELLIPTIC "okl/nonZero.h", DELLIPTIC "okl/nonZeroCompare.h", kernelInfo);

  platform.device.finish();
  double t3 = MPI_Wtime();
  
  // 1. sort based on row (fastest) then column in each row
  sorter.sort(Nnz, o_AL);

  platform.device.finish();
  double t4 = MPI_Wtime();

  int parallelCompareRowColumn(const void *a, const void *b);
  if(testHOST)
    qsort(h_AL, Nnz, sizeof(nnz_t), parallelCompareRowColumn);

  double t5 = MPI_Wtime();
  if(testHOST)
    std::sort(h_AL2, h_AL2+Nnz, parallelCompareRowColumnV2);
  double t6 = MPI_Wtime();

  // 2. perform scan  to find unique entries
  occa::memory o_tmp;
  dlong  *h_tmp;
  occa::memory o_scan = platform.device.malloc(Nnz*sizeof(dlong));
  
  scanner.mallocTemps(platform, Nnz, o_tmp, &h_tmp);

  platform.device.finish();
  double t7 = MPI_Wtime();

  scanner.scan(Nnz, o_AL, o_tmp, h_tmp, o_scan);
  
  platform.device.finish();
  double t8 = MPI_Wtime();

  occa::memory o_compactedAL;
  int includeLast = (elliptic.allNeumann);
  dlong compactedNnz = scanner.trashCompactor(platform, Nnz, sizeof(nnz_t), includeLast, o_AL, o_compactedAL);
  platform.device.finish();
  double t9 = MPI_Wtime();

  printf("compactedNnz=%d\n", compactedNnz);
  
  // 3. use block reduce with 32 threads to compress entries

  // 3.a extract starts
  // 3.b compress

  printf("DEVICE    sort took: %e\n", t4-t3);
  if(testHOST){
    printf("HOST     qsort took: %e\n", t5-t4);
    printf("HOST std::sort took: %e\n", t6-t5);
  }
  
  printf("DEVICE: sorted %e gdofs at a rate of %e gdofs/s\n", Nnz/1.e9, Nnz/(1.e9*(t4-t3)));
  if(testHOST){
    printf("HOST: qsorted %e gdofs at a rate of %e gdofs/s\n", Nnz/1.e9, Nnz/(1.e9*(t5-t4)));
    printf("HOST: std::sorted %e gdofs at a rate of %e gdofs/s\n", Nnz/1.e9, Nnz/(1.e9*(t6-t5)));
  }
  
  printf("DEVICE    scan took: %e\n", t8-t7);
  printf("DEVICE trash compactor: %e\n", t9-t8);
  

  if(testHOST){
    // check scan worked
    nnz_t *d_compactedAL = (nnz_t*) calloc(compactedNnz, sizeof(nnz_t));
    o_compactedAL.copyTo(d_compactedAL);
    
    dfloat tol = 1e-10;
    for(int n=0;n<A.nnz;++n){
      dfloat d = A.entries[n].val -  d_compactedAL[n].val;
      if(A.entries[n].row != d_compactedAL[n].row ||
	 A.entries[n].col != d_compactedAL[n].col ||
	 d*d>tol){
	printf("mismatch: %d,%d,%e => %d,%d,%e\n",
	       A.entries[n].row,
	       A.entries[n].col,
	       A.entries[n].val,
	       d_compactedAL[n].row,
	       d_compactedAL[n].col,
	       d_compactedAL[n].val);
	
      }
    }
  }
       
  
#if 0
  for(int n=0;n<compactedNnz;++n){
    nnz_t ent = d_compactedAL[n];
    printf("d_compactedAL[%d] = [%d,%d,%g]\n", n, ent.row, ent.col, ent.val);
  }
#endif


#endif
  
  // close down MPI
  MPI_Finalize();
  return LIBP_SUCCESS;
}
