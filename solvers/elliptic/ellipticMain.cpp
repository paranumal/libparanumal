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
#define nonZero_t parAlmond::parCOO::nonZero_t

// compare on global indices
bool parallelCompareRowColumnV2(nonZero_t &a, nonZero_t &b){
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

  // build scanner and sorter
  occa::properties kernelInfo = mesh.props;

  if(sizeof(hlong)==8)
    kernelInfo["defines/hlong"]= "long long int";
  if(sizeof(hlong)==4)
    kernelInfo["defines/hlong"]= "int";
  
  deviceScan_t scanner(platform, DELLIPTIC "okl/nonZero.h", DELLIPTIC "okl/nonZeroCompare2.h", kernelInfo);
  deviceSort_t  sorter(platform, DELLIPTIC "okl/nonZero.h", DELLIPTIC "okl/nonZeroCompare.h", kernelInfo);

  hlong BIG_NUM = 1 << (8*sizeof(hlong)-2);
  kernelInfo["defines/" "BIG_NUM"] = (hlong)BIG_NUM;
  
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
    platform.buildKernel(DELLIPTIC "/okl/ellipticBuildOperatorMatrixContinuous.okl",
			 kernelName,
			 kernelInfo);
  
  occa::memory o_maskedGlobalNumbering =
    platform.malloc(mesh.Np*mesh.Nelements*sizeof(hlong), elliptic.maskedGlobalNumbering);
  
  platform.device.finish();
  double t10 = MPI_Wtime();
  parAlmond::parCOO Ahost(elliptic.platform, mesh.comm);
  elliptic.BuildOperatorMatrixContinuous(Ahost);

  parAlmond::parCOO Adev(elliptic.platform, mesh.comm);
  occa::memory o_A;
  dlong devAnnz;

  platform.device.finish();
  double t11 = MPI_Wtime();
    
  elliptic.BuildOperatorMatrixContinuousDevice(buildMatrixKernel, o_maskedGlobalNumbering, BIG_NUM, sorter, scanner, Adev, o_A, devAnnz);
  
  platform.device.finish();
  double t12 = MPI_Wtime();

  double elapsedHost = t11-t10;
  double elapsedDevice = t12-t11;

  double globalElapsedHost, globalElapsedDevice;
  int root = 0;
  
  MPI_Reduce(&elapsedHost, &globalElapsedHost, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
  MPI_Reduce(&elapsedDevice, &globalElapsedDevice, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);

  if(mesh.rank==root){
    printf("DEVICE build whole matrix (host): %e\n", elapsedHost);
    printf("DEVICE build whole matrix (dev): %e\n",  elapsedDevice);
  }
  
  
  Adev.nnz = devAnnz;
  Adev.entries = (nonZero_t*) calloc(Adev.nnz, sizeof(nonZero_t));
  o_A.copyTo(Adev.entries);
  
  int testHOST = 1;
  if(testHOST){

    if(Ahost.nnz!=Adev.nnz){ printf("mismatch in HOST and DEVICE non-zero count: %d to %d\n",
				    Ahost.nnz, Adev.nnz); }
    
    dfloat tol = 1e-15;
    for(int n=0;n<mymin(10,mymin(Ahost.nnz,Adev.nnz));++n){
      nonZero_t Ahostn = Ahost.entries[n];
      nonZero_t Adevn  = Adev.entries[n];
      
      dfloat d = Ahostn.val -  Adevn.val;
      if(Ahostn.row != Adevn.row ||
	 Ahostn.col  != Adevn.col ||
	 d*d>tol){

	printf("mismatch: %d,%d,%e => %d,%d,%e\n",
	       Ahostn.row, Ahostn.col,  Ahostn.val,
	       Adevn.row,  Adevn.col,   Adevn.val);
      }
    }
  }
  
  // close down MPI
  MPI_Finalize();
  return LIBP_SUCCESS;
}
