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

#if 0
typedef struct{
  hlong gnum;
  int   lnum;
}globalNode_t;

int compareGlobalNodes2(const void *a, const void *b){
  
  globalNode_t *ea = (globalNode_t*) a;
  globalNode_t *eb = (globalNode_t*) b;
  
  if(ea->gnum < eb->gnum) return -1;
  if(ea->gnum > eb->gnum) return +1;
  
  return 0;
}
#endif

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

  platform.device.finish();

  parAlmond::parCOO Ahost(elliptic.platform, mesh.comm);

  int testHOST = 0;

  if(testHOST){
    double t1 = MPI_Wtime();
    elliptic.BuildOperatorMatrixContinuous(Ahost);
    double t2 = MPI_Wtime();
    printf("HOST took %g secs to build matrix\n", t2-t1);
  }
  occa::memory o_A;
  dlong devAnnz;

  platform.device.finish();
  double t3 = MPI_Wtime();
  elliptic.BuildOperatorMatrixContinuousDevice(o_A, devAnnz);

  platform.device.finish();
  double t4 = MPI_Wtime();
  printf("DEVICE took %g secs to build matrix\n", t4-t3);

  if(testHOST){
    nonZero_t *h_A = (nonZero_t*) calloc(devAnnz, sizeof(nonZero_t));
    o_A.copyTo(h_A);
    
    dfloat tol = 1e-15;
    for(int n=0;n<mymin(Ahost.nnz,devAnnz);++n){
      nonZero_t Ahostn = Ahost.entries[n];
      nonZero_t Adevn  = h_A[n];
      
      dfloat d = Ahostn.val -  Adevn.val;
      if(Ahostn.row != Adevn.row ||
	 Ahostn.col  != Adevn.col ||
	 d*d>tol){

	printf("mismatch: (host) %d,%d,%e => %d,%d,%e (dev)\n",
	       Ahostn.row, Ahostn.col,  Ahostn.val,
	       Adevn.row,  Adevn.col,   Adevn.val);
      }
    }

  }
  
  // close down MPI
  MPI_Finalize();
  return LIBP_SUCCESS;
}
