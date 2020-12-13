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

  hlong BIG_NUM = ((hlong)1) << (8*sizeof(hlong)-2);
  kernelInfo["defines/" "BIG_NUM"] = (const int64_t)BIG_NUM;
  
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

  int testHOST = 0;

  if(testHOST)
    elliptic.BuildOperatorMatrixContinuous(Ahost);

  parAlmond::parCOO Adev(elliptic.platform, mesh.comm);
  occa::memory o_A;
  dlong devAnnz;

  // START ONE-RING VERSION
  platform.device.finish();
  double ticA = MPI_Wtime();

  // A. build one ring mesh
  mesh.HaloRingSetup();

  // B. transfer into halo
  dlong Nel = mesh.Nelements;
  dlong ringNel = mesh.totalRingElements;
  dlong allNel = Nel+ringNel;

  // i. import ggeo
  size_t NggeoBlk;
  if(mesh.elementType == TRIANGLES ||
     mesh.elementType == TETRAHEDRA)
    NggeoBlk = mesh.Nggeo;
  else
    NggeoBlk = mesh.Nggeo*mesh.Np;
  
  occa::memory o_ringGgeo =
    platform.device.malloc(allNel*NggeoBlk*sizeof(dfloat));

  o_ringGgeo.copyFrom(mesh.o_ggeo, Nel*NggeoBlk*sizeof(dfloat), 0, 0);
					   
  mesh.ringHalo->Exchange(o_ringGgeo, NggeoBlk, ogs_dfloat);

  // ii. import globalNumbers
  size_t NnumBlk = mesh.Np;
  occa::memory o_ringMaskedGlobalNumbering =
    platform.device.malloc(allNel*NnumBlk*sizeof(hlong));

  o_ringMaskedGlobalNumbering.copyFrom(o_maskedGlobalNumbering, Nel*NnumBlk*sizeof(hlong), 0, 0);

  mesh.ringHalo->Exchange(o_ringMaskedGlobalNumbering, NnumBlk, ogs_hlong);

  hlong *ringMaskedGlobalNumbering = (hlong*) calloc(allNel*NnumBlk,sizeof(hlong));
  o_ringMaskedGlobalNumbering.copyTo(ringMaskedGlobalNumbering);
  
  // iii. build output map
  hlong Ngather = elliptic.ogsMasked->Ngather;
  hlong *globalRowStarts = (hlong*) calloc(mesh.size+1,sizeof(hlong));
  MPI_Allgather(&Ngather, 1, MPI_HLONG, globalRowStarts+1, 1, MPI_HLONG, mesh.comm);
  for(int r=0;r<mesh.size;++r) {
    globalRowStarts[r+1] = globalRowStarts[r]+globalRowStarts[r+1];
  }

  hlong gatherStart = globalRowStarts[mesh.rank];
  hlong gatherEnd   = globalRowStarts[mesh.rank+1];

  // for each element build map to shuffle output
  globalNode_t *gnodes = (globalNode_t*) calloc(mesh.Np, sizeof(globalNode_t));
  dlong    *colMap = (dlong*) calloc(mesh.Np*allNel, sizeof(dlong));
  for(dlong e=0;e<allNel;++e){
    for(int n=0;n<mesh.Np;++n){
      hlong id = ringMaskedGlobalNumbering[e*mesh.Np+n];
      gnodes[n].gnum = (id<0) ? BIG_NUM:id;
      gnodes[n].lnum = n;
    }
    qsort(gnodes, mesh.Np, sizeof(globalNode_t), compareGlobalNodes2);
    for(int n=0;n<mesh.Np;++n){
      colMap[e*mesh.Np+gnodes[n].lnum] = n; // order to map write index into
    }
#if 0
    for(int n=0;n<mesh.Np;++n){
      printf("colMap[%d] = %d\n", n, colMap[e*mesh.Np+n]);
    }
#endif
  }

  occa::memory o_colMap = platform.device.malloc(allNel*mesh.Np*sizeof(dlong), colMap);


  // find where to put the row segments
  dlong *rowCount = (dlong*) calloc(mesh.Np*allNel, sizeof(dlong));
  for(dlong e=0;e<allNel;++e){
    for(int n=0;n<mesh.Np;++n){
      hlong id = ringMaskedGlobalNumbering[e*mesh.Np+n];
      if(gatherStart<=id && id<gatherEnd){
	++(rowCount[id-gatherStart]);
      }
    }
  }
  dlong *rowStarts = (dlong*) calloc(Ngather+1, sizeof(dlong));
  for(dlong n=0;n<Ngather;++n){
    rowStarts[n+1] = rowStarts[n] + rowCount[n];
  }
				     
  dlong *rowMap = (dlong*) calloc(mesh.Np*allNel, sizeof(dlong));
  memset(rowMap, -1, mesh.Np*allNel*sizeof(dlong));
  for(dlong e=0;e<allNel;++e){
    for(int n=0;n<mesh.Np;++n){

      hlong id = ringMaskedGlobalNumbering[e*mesh.Np+n];
      //      printf("globalNum=%lld\n", id);
      if(gatherStart<=id && id<gatherEnd){
	rowMap[e*mesh.Np+n] = (rowStarts[id-gatherStart]++);
	//	printf("rowMap = %d\n", rowMap[e*mesh.Np+n]);
      }
    }
  }

  dlong maxDegree = 0;
  for(dlong g=0;g<Ngather;++g){
    dlong deg = (g==0) ? rowStarts[g] : rowStarts[g]-rowStarts[g-1];
    ///    printf("rowCount[%d] = %d\n", g,deg);
    maxDegree = mymax(maxDegree, deg);
  }
  printf("maxDegree=%d\n", maxDegree);

  kernelInfo["defines/" "MAX_DEGREE"]= (int) maxDegree;
  
  occa::kernel assembleMatrixKernel = 
    platform.buildKernel(DELLIPTIC "/okl/ellipticBuildOperatorMatrixContinuous.okl",
			 "ellipticAssembleEntries",
			 kernelInfo);
  
  occa::memory o_rowStarts = platform.device.malloc(Ngather*sizeof(dlong), rowStarts);

  dlong allNgather = rowStarts[Ngather-1];
  printf("gatherStart = %lld, gatherEnd = %lld, Ngather=%d, allNgather=%d\n",
	 gatherStart, gatherEnd, Ngather, allNgather*mesh.Np);
  
  occa::memory o_rowMap = platform.device.malloc(allNel*mesh.Np*sizeof(dlong), rowMap);
  
  occa::memory o_AL = platform.device.malloc(allNgather*mesh.Np*sizeof(nonZero_t));
  occa::memory o_rowCounts = platform.malloc(Ngather*sizeof(dlong));
  occa::memory o_AL2 = platform.malloc(allNgather*mesh.Np*sizeof(nonZero_t)); // wasteful until fixed
  
  switch(mesh.elementType){
  case TRIANGLES:
    buildMatrixKernel(allNel, o_ringMaskedGlobalNumbering,  o_rowMap, o_colMap,
		      mesh.o_S, mesh.o_MM, o_ringGgeo, 
		      lambda, o_AL);

    break;
  case QUADRILATERALS:
    buildMatrixKernel(allNel, o_ringMaskedGlobalNumbering, o_rowMap, o_colMap,
		      mesh.o_D, o_ringGgeo,
		      lambda, o_AL);

    break;
  case TETRAHEDRA:
    buildMatrixKernel(allNel, o_ringMaskedGlobalNumbering, o_rowMap, o_colMap,
		      mesh.o_S, mesh.o_MM, o_ringGgeo, 
		      lambda, o_AL);
    break;
  case HEXAHEDRA:
    buildMatrixKernel(allNel, o_ringMaskedGlobalNumbering, o_rowMap, o_colMap,
		      mesh.o_D,  o_ringGgeo,
		      lambda, o_AL);
    break;
  }

  assembleMatrixKernel(Ngather, o_rowStarts, o_AL, o_AL2, o_rowCounts);
  
  platform.device.finish();
  double ticB = MPI_Wtime();

#if 1
  if(testHOST){
    // just for testing
    dlong *rowCounts = (dlong*) calloc(Ngather, sizeof(dlong));
    o_rowCounts.copyTo(rowCounts);

    nonZero_t *h_AL2 = (nonZero_t*) calloc(allNgather*mesh.Np, sizeof(nonZero_t));
    nonZero_t *h_AL3 = (nonZero_t*) calloc(allNgather*mesh.Np, sizeof(nonZero_t));
    o_AL2.copyTo(h_AL2);

    //  dlong cnt = 0;
    dlong skip = 0;
    for(dlong n=0;n<Ngather;++n){
      dlong start = (n==0) ? 0:rowStarts[n-1]*mesh.Np;
      dlong end = rowStarts[n]*mesh.Np;
      dlong cnt = rowCounts[n];
      //    printf("Node %d (orign = %d, rowCount = %d) \n", n, end-start, cnt);
      for(dlong m=start;m<start+cnt;++m){
	//      if(!(m%mesh.Np)) printf("---\n");
	//      printf("%lld %lld %g\n", h_AL[m].row, h_AL[m].col, h_AL[m].val);
	h_AL3[skip++] = h_AL2[m];
      }
    }

    for(dlong n=0;n<mymin(skip, Ahost.nnz);++n){
      nonZero_t Ahostn = Ahost.entries[n];
      nonZero_t Adevn  = h_AL3[n];
      dfloat tol = 1e-15;
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
#endif
  printf("elapsed: %f\n", ticB-ticA);
  exit(-1);


  // END 
  


  
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
  
  if(testHOST){

    if(Ahost.nnz!=Adev.nnz){ printf("mismatch in HOST and DEVICE non-zero count: %d to %d\n",
				    Ahost.nnz, Adev.nnz); exit(-1); }
    
    dfloat tol = 1e-15;
    for(int n=0;n<mymin(Ahost.nnz,Adev.nnz);++n){
      nonZero_t Ahostn = Ahost.entries[n];
      nonZero_t Adevn  = Adev.entries[n];
      
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
