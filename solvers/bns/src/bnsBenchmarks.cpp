
/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "bns.hpp"

void bns_t::Benchmarks(){

  deviceMemory<dfloat> o_Q   = platform.device.malloc(Nfields*mesh.Nelements*mesh.Np*sizeof(dfloat));
  deviceMemory<dfloat> o_RHS = platform.device.malloc(Nfields*mesh.Nelements*mesh.Np*sizeof(dfloat));


  printf("Starting volumBenchmark\n");
  volumeBenchmark(o_Q, o_RHS);

#if 0
  printf("Starting surfaceBernsteinBenchmark\n");
  surfaceBernsteinBenchmark(o_Q, o_RHS);
  
  printf("Starting volumeBernsteinBenchmark\n");
  volumeBernsteinBenchmark(o_Q, o_RHS);

  printf("Starting surfaceBenchmark\n");
  surfaceBenchmark(o_Q, o_RHS);
#endif


}


void bns_t::volumeBenchmark(deviceMemory<dfloat> &o_Q, deviceMemory<dfloat> &o_RHS){

  
  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  const dfloat p_half = 1./2.;
  int blockMax = 768; // 512;
  
  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;
  kernelInfo["defines/" "p_Nfields"]= Nfields;
  kernelInfo["defines/" "p_half"]= p_half;
  kernelInfo["defines/" "p_Npmlfields"]= Npmlfields;
  kernelInfo["defines/" "p_sqrt2"]= (dfloat)sqrt(2.0);
  
  // set kernel name suffix
  std::string suffix;
  if(mesh.elementType==Mesh::TRIANGLES)
    suffix = "Tri2D";
  if(mesh.elementType==Mesh::QUADRILATERALS)
    suffix = "Quad2D";
  if(mesh.elementType==Mesh::TETRAHEDRA)
    suffix = "Tet3D";
  if(mesh.elementType==Mesh::HEXAHEDRA)
    suffix = "Hex3D";

  std::string oklFilePrefix = DBNS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string volumeFileName, volumeKernelName;

  occa::streamTag start, stop;

  // kernels from volume file
  volumeFileName   = oklFilePrefix + "bnsVolume" + suffix + oklFileSuffix;
  volumeKernelName = "bnsVolume" + suffix;

  // set up flop counts for tets
  printf("Nelements = %d, NnonPmlElements = %d\n", mesh.Nelements, mesh.NnonPmlElements);
  
  long long int NFLOP   = mesh.NnonPmlElements*mesh.Np*((Nfields*mesh.dim*2*(long long int)mesh.Np) + 32);
  long long int NBYTES  = mesh.NnonPmlElements*(mesh.Nvgeo + mesh.Np*(Nfields*2 + mesh.dim) )*sizeof(dfloat);

  int bestNvol = 0, bestNblockV = 0, bestKnl = 0;
  double bestElapsed = 1e9;;

  dfloat T = 1;
  int maxNvol = 1;
  for(int Nvol=1;Nvol<=maxNvol;++Nvol){
    for(int NblockV=1;NblockV<=blockMax/mesh.Np;++NblockV){
      int LDS = (mesh.Np*Nvol*Nfields*NblockV+Nvol*mesh.Nvgeo*NblockV)*sizeof(dfloat);
      // limit case based on shared array storage
      if(LDS<36*1024){
	int Nkernels = 2;
	for(int knl=0;knl<Nkernels;++knl){
	  properties_t volumeKernelInfo = kernelInfo;
	  
	  volumeKernelInfo["defines/" "p_NblockV"]= NblockV;
	  volumeKernelInfo["defines/" "p_Nvol"]= Nvol;
	  volumeKernelInfo["defines/" "p_knl"]= knl;
	  
	  platform.device.finish();
	
	  volumeKernel =  platform.buildKernel(volumeFileName, volumeKernelName,
					       volumeKernelInfo);
	

	  //  rhsVolume(mesh.NnonPmlElements, mesh.o_nonPmlElements, o_Q, o_RHS, T);
	  
	  int Nwarm = 1;
	  for(int w=0;w<Nwarm;++w){
	    volumeKernel(mesh.NnonPmlElements,
			 mesh.o_nonPmlElements,
			 mesh.o_vgeo,
			 mesh.o_D,
			 mesh.o_x,
			 mesh.o_y,
			 mesh.o_z,
			 T,
			 c,
			 nu,
			 o_Q,
			 o_RHS);
	  }
	
	  platform.device.finish();
	
	  start = platform.device.tagStream();
	
	  int Nrun = 2;
	  for(int r=0;r<Nrun;++r){
	    volumeKernel(mesh.NnonPmlElements,
			 mesh.o_nonPmlElements,
			 mesh.o_vgeo,
			 mesh.o_D,
			 mesh.o_x,
			 mesh.o_y,
			 mesh.o_z,
			 T,
			 c,
			 nu,
			 o_Q,
			 o_RHS);
	  }
	
	  stop = platform.device.tagStream();
	
	  platform.device.finish();
	
	  double elapsed = platform.device.timeBetween(start, stop);
	  elapsed /= Nrun;
	
	  double GFLOPS = NFLOP/(1.e9*elapsed);
	  double GBS = NBYTES/(1.e9*elapsed);
	
	  printf("%02d, %02d, %02d, %02d, %5.4e, %5.4e, %5.4e, %lld, %d %%%% VOL: N, Knl, Nvol, NblockV, elapsed, GFLOP/s, GB/s, NFLOP, LDS usage\n",
		 mesh.N, knl, Nvol, NblockV, elapsed, GFLOPS, GBS, NFLOP, LDS);
	
	  if(elapsed<bestElapsed){
	    bestElapsed = elapsed;
	    bestNvol = Nvol;;
	    bestNblockV = NblockV;
	    bestKnl = knl;
	  }
	}
      }
    }  
  }

  double bestGFLOPS = NFLOP/(1.e9*bestElapsed);
  double bestGBS = NBYTES/(1.e9*bestElapsed);
  
  printf("%02d, %02d, %02d, %02d, %5.4e, %5.4e, %5.4e %%%% VOL, BEST: N, Knl, Nvol, NblockV, elapsed, GFLOP/s, GB/s\n",
	 mesh.N, bestKnl, bestNvol, bestNblockV, bestElapsed, bestGFLOPS, bestGBS);

}



void bns_t::surfaceBenchmark(deviceMemory<dfloat> &o_Q, deviceMemory<dfloat> &o_RHS){

  
  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties
  
  const dfloat p_half = 1./2.;
  int blockMax = 256;
  
  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;
  kernelInfo["defines/" "p_Nfields"]= Nfields;
  kernelInfo["defines/" "p_half"]= p_half;
  kernelInfo["defines/" "p_NfacesNfp"]= mesh.Nfaces*mesh.Nfp;

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;
  
  // set kernel name suffix
  std::string suffix;
  if(mesh.elementType==Mesh::TRIANGLES)
    suffix = "Tri2D";
  if(mesh.elementType==Mesh::QUADRILATERALS)
    suffix = "Quad2D";
  if(mesh.elementType==Mesh::TETRAHEDRA)
    suffix = "Tet3D";
  if(mesh.elementType==Mesh::HEXAHEDRA)
    suffix = "Hex3D";

  std::string oklFilePrefix = DBNS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string surfaceFileName, surfaceKernelName;

  occa::streamTag start, stop;

  // kernels from surface file
  surfaceFileName   = oklFilePrefix + "bnsSurface" + suffix + oklFileSuffix;
  surfaceKernelName = "bnsSurface" + suffix;

  // set up flop counts for tets
  long long int NFLOP = ((long long int)mesh.Nelements)*(mesh.Nfp*mesh.Nfaces*(10 + 18 + 4) + mesh.Np*((1+mesh.Nfaces*mesh.Nfp)*Nfields*2));

  long long int NBYTES   = mesh.Nelements*(
					   (mesh.Nsgeo*mesh.Nfaces  +
					    mesh.Nfp*mesh.Nfaces*Nfields*2 +
					    mesh.Np*(Nfields*2) )*sizeof(dfloat)
					   + (mesh.Nfp*mesh.Nfaces*2*sizeof(dlong)));
												       
  int bestNblockS = 0;
  int bestNsur = 0;
  double bestElapsed = 1e9;;
  dfloat T = 0;

  for(int Nsur=1;Nsur<=5;++Nsur){
    for(int NblockS=1;NblockS<=blockMax/maxNodes;++NblockS){
      properties_t surfaceKernelInfo = kernelInfo;
      
      surfaceKernelInfo["defines/" "p_NblockS"]= NblockS;
      surfaceKernelInfo["defines/" "p_Nsur"]= Nsur;
      
      platform.device.finish();
      
      surfaceKernel =  platform.buildKernel(surfaceFileName, surfaceKernelName,
					    surfaceKernelInfo);
      
      
      int Nwarm = 100;
      for(int w=0;w<Nwarm;++w){
	
	if (mesh.NinternalElements)
	  surfaceKernel(mesh.NinternalElements,
			mesh.o_internalElementIds,
			mesh.o_sgeo,
			mesh.o_LIFT,
			mesh.o_vmapM,
			mesh.o_vmapP,
			mesh.o_EToB,
			T,
			mesh.o_x,
			mesh.o_y,
			mesh.o_z,
			o_Q,
			o_RHS);
	
      }
      
      platform.device.finish();
      
      start = platform.device.tagStream();
      
      
      int Nrun = 10;
      for(int r=0;r<Nrun;++r){
	if (mesh.NinternalElements)
	  surfaceKernel(mesh.NinternalElements,
			mesh.o_internalElementIds,
			mesh.o_sgeo,
			mesh.o_LIFT,
			mesh.o_vmapM,
			mesh.o_vmapP,
			mesh.o_EToB,
			T,
			mesh.o_x,
			mesh.o_y,
			mesh.o_z,
			o_Q,
			o_RHS);
      }
      
      stop = platform.device.tagStream();
      
      platform.device.finish();
      
      double elapsed = platform.device.timeBetween(start, stop);
      elapsed /= Nrun;
      
      double GFLOPS = NFLOP/(1.e9*elapsed);
      double GBS = NBYTES/(1.e9*elapsed);
      
      printf("%02d, %02d, %5.4e, %5.4e, %5.4e %%%% SURF, Nsur, NblockS, elapsed, GFLOP/s, GB/s\n",
	     Nsur, NblockS, elapsed, GFLOPS, GBS);
      
      if(elapsed<bestElapsed){
	bestElapsed = elapsed;
	bestNblockS = NblockS;
	bestNsur = Nsur;
      }
    }
  }

  
  double bestGFLOPS = NFLOP/(1.e9*bestElapsed);
  double bestGBS = NBYTES/(1.e9*bestElapsed);
  
  printf("%02d, %02d, %02d, %5.4e, %5.4e, %5.4e %%%% SURF BEST - N, Nsur, NblockS, elapsed, GFLOP/s, GB/s\n",
	 mesh.N, bestNsur, bestNblockS, bestElapsed, bestGFLOPS, bestGBS);
  
}


// do something smarter here later
#if 0
#include "data3dN04.h"
#include "data3dN05.h"
#else
//#include "data2dN04.h"
//#include "data2dN06.h"
#include "data3dN04.h"
#endif

typedef struct {
  unsigned char x;
  unsigned char y;
  unsigned char z;
  unsigned char w;
}uchar4;

void bns_t::volumeBernsteinBenchmark(deviceMemory<dfloat> &o_Q, deviceMemory<dfloat> &o_RHS){

  // load Bernstein stuff
  deviceMemory<int> o_D1_ids = platform.device.malloc(mesh.Np*sizeof(int)*4, p_D1_ids[0]);
  deviceMemory<int> o_D2_ids = platform.device.malloc(mesh.Np*sizeof(int)*4, p_D2_ids[0]);
  deviceMemory<int> o_D3_ids = platform.device.malloc(mesh.Np*sizeof(int)*4, p_D3_ids[0]);
  deviceMemory<int> o_D4_ids = platform.device.malloc(mesh.Np*sizeof(int)*4, p_D4_ids[0]);
  deviceMemory<dfloat> o_D_vals = platform.device.malloc(mesh.Np*sizeof(dfloat)*4, p_D_vals);

  uchar4 *c_D1_ids = (uchar4*) calloc(mesh.Np, sizeof(uchar4));
  uchar4 *c_D2_ids = (uchar4*) calloc(mesh.Np, sizeof(uchar4));
  uchar4 *c_D3_ids = (uchar4*) calloc(mesh.Np, sizeof(uchar4));
  uchar4 *c_D4_ids = (uchar4*) calloc(mesh.Np, sizeof(uchar4));

  for(int n=0;n<mesh.Np;++n){
    c_D1_ids[n].x = p_D1_ids[n][0];
    c_D1_ids[n].y = p_D1_ids[n][1];
    c_D1_ids[n].z = p_D1_ids[n][2];
    c_D1_ids[n].w = (mesh.dim==3) ? p_D1_ids[n][3]: 0;

    c_D2_ids[n].x = p_D2_ids[n][0];
    c_D2_ids[n].y = p_D2_ids[n][1];
    c_D2_ids[n].z = p_D2_ids[n][2];
    c_D2_ids[n].w = (mesh.dim==3) ? p_D2_ids[n][3]: 0;
   
    c_D3_ids[n].x = p_D3_ids[n][0];
    c_D3_ids[n].y = p_D3_ids[n][1];
    c_D3_ids[n].z = p_D3_ids[n][2];
    c_D3_ids[n].w = (mesh.dim==3) ?p_D3_ids[n][3]: 0;
    
    c_D4_ids[n].x = p_D4_ids[n][0];
    c_D4_ids[n].y = p_D4_ids[n][1];
    c_D4_ids[n].z = p_D4_ids[n][2];
    c_D4_ids[n].w = (mesh.dim==3) ?p_D4_ids[n][3]:0;
  }
  
  deviceMemory<uchar4> o_cD1_ids = platform.device.malloc(mesh.Np*sizeof(uchar4), c_D1_ids);
  deviceMemory<uchar4> o_cD2_ids = platform.device.malloc(mesh.Np*sizeof(uchar4), c_D2_ids);
  deviceMemory<uchar4> o_cD3_ids = platform.device.malloc(mesh.Np*sizeof(uchar4), c_D3_ids);
  deviceMemory<uchar4> o_cD4_ids = platform.device.malloc(mesh.Np*sizeof(uchar4), c_D4_ids);

  
  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  const dfloat p_half = 1./2.;
  int blockMax = 768; // 512;
  
  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;
  kernelInfo["defines/" "p_Nfields"]= Nfields;
  kernelInfo["defines/" "p_half"]= p_half;

  // set kernel name suffix
  std::string suffix;
  if(mesh.elementType==Mesh::TRIANGLES)
    suffix = "TriBB2D";
  if(mesh.elementType==Mesh::QUADRILATERALS)
    suffix = "QuadBB2D";
  if(mesh.elementType==Mesh::TETRAHEDRA)
    suffix = "TetBB3D";
  if(mesh.elementType==Mesh::HEXAHEDRA)
    suffix = "HexBB3D";

  std::string oklFilePrefix = DBNS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string volumeFileName, volumeKernelName;

  occa::streamTag start, stop;

  // kernels from volume file
  volumeFileName   = oklFilePrefix + "bnsVolume" + suffix + oklFileSuffix;
  volumeKernelName = "bnsVolume" + suffix;

  // set up flop counts for tets
  long long int NFLOP = mesh.Nelements*mesh.Np*((4*3*2*(long long int)mesh.Np) + 32);
  long long int NBYTES   = mesh.Nelements*(mesh.Nvgeo + mesh.Np*(Nfields*2 + 0*mesh.dim*mesh.Np) )*sizeof(dfloat);

  int bestNvol = 0, bestNblockV = 0, bestKnl = 0;
  double bestElapsed = 1e9;;
  
  for(int Nvol=1;Nvol<=8;++Nvol){
    for(int NblockV=1;NblockV<=blockMax/mesh.Np;++NblockV){
      int LDS = (mesh.Np*Nvol*Nfields*NblockV+Nvol*mesh.Nvgeo*NblockV)*sizeof(dfloat);
      // limit case based on shared array storage
      if(LDS<48*1024){

	int Nkernels = 3;
	for(int knl=2;knl<Nkernels;++knl){
	  
	  properties_t volumeKernelInfo = kernelInfo;
	
	  volumeKernelInfo["defines/" "p_NblockV"]= NblockV;
	  volumeKernelInfo["defines/" "p_Nvol"]= Nvol;
	  volumeKernelInfo["defines/" "p_knl"]= knl;
	
	  platform.device.finish();
	
	  volumeKernel =  platform.buildKernel(volumeFileName, volumeKernelName,
					       volumeKernelInfo);
	
	
	  int Nwarm = 20;
	  for(int w=0;w<Nwarm;++w){
	    volumeKernel(mesh.Nelements,
			 mesh.o_vgeo,
			 o_cD1_ids, o_cD2_ids, o_cD3_ids, o_cD4_ids,
			 o_D_vals,
			 o_Q,
			 o_RHS);
	  }
	
	  platform.device.finish();
	
	  start = platform.device.tagStream();
	
	
	  int Nrun = 20;
	  for(int r=0;r<Nrun;++r){
	    volumeKernel(mesh.Nelements,
			 mesh.o_vgeo,
			 o_cD1_ids, o_cD2_ids, o_cD3_ids, o_cD4_ids,
			 o_D_vals,
			 o_Q,		 
			 o_RHS);
	  }
	
	  stop = platform.device.tagStream();
	
	  platform.device.finish();
	
	  double elapsed = platform.device.timeBetween(start, stop);
	  elapsed /= Nrun;
	
	  double GFLOPS = NFLOP/(1.e9*elapsed);
	  double GBS = NBYTES/(1.e9*elapsed);
	
	  printf("%02d, %02d, %02d, %02d, %5.4e, %5.4e, %5.4e, %lld, %d %%%% BB-VOL, Nvol, NblockV, elapsed, GFLOP/s, GB/s, NFLOP, LDS usage\n",
		 mesh.N, knl, Nvol, NblockV, elapsed, GFLOPS, GBS, NFLOP, LDS);
	
	  if(elapsed<bestElapsed){
	    bestElapsed = elapsed;
	    bestNvol = Nvol;;
	    bestNblockV = NblockV;
	    bestKnl = knl;
	  }
	}
      }
    }  
  }
  double bestGFLOPS = NFLOP/(1.e9*bestElapsed);
  double bestGBS = NBYTES/(1.e9*bestElapsed);
  
  printf("%02d, %02d, %02d, %02d, %5.4e, %5.4e, %5.4e %%%% BB-VOL, BESTL N, Knl, Nvol, NblockV, elapsed, GFLOP/s, GB/s\n",
	 mesh.N, bestKnl, bestNvol, bestNblockV, bestElapsed, bestGFLOPS, bestGBS);

}


void bns_t::surfaceBernsteinBenchmark(deviceMemory<dfloat> &o_Q, deviceMemory<dfloat> &o_RHS){

  int minNfp7 = std::min(7,mesh.Nfp);
  
  dfloat *tmp_EEL_vals = (dfloat*) calloc(mesh.Np*p_EEL_nnz, sizeof(dfloat));
  int    *tmp_EEL_ids  = (int*)    calloc(mesh.Np*p_EEL_nnz, sizeof(int));
  dfloat *tmp_L0_vals  = (dfloat*) calloc(mesh.Nfp*minNfp7, sizeof(dfloat));
  int    *tmp_L0_ids   = (int*)    calloc(mesh.Nfp*minNfp7, sizeof(int));

  for(int n=0;n<mesh.Np;++n){
    for(int m=0;m<p_EEL_nnz;++m){
      tmp_EEL_vals[n+m*mesh.Np] = p_EEL_vals[n][m];
      tmp_EEL_ids[n+m*mesh.Np] = p_EEL_ids[n][m];
    }
  }

  for(int n=0;n<mesh.Nfp;++n){
    for(int m=0;m<minNfp7;++m){
      tmp_L0_vals[n+m*mesh.Nfp] = p_L0_vals[n][m];
      tmp_L0_ids[n+m*mesh.Nfp] = p_L0_ids[n][m];
    }
  }
  
  // BB stuff
  deviceMemory<dfloat> o_EEL_vals = platform.device.malloc(mesh.Np*p_EEL_nnz*sizeof(dfloat), tmp_EEL_vals);
  deviceMemory<int>    o_EEL_ids  = platform.device.malloc(mesh.Np*p_EEL_nnz*sizeof(int),    tmp_EEL_ids);
  deviceMemory<dfloat> o_L0_vals  = platform.device.malloc(mesh.Nfp*minNfp7*sizeof(dfloat),        tmp_L0_vals);
  deviceMemory<int>    o_L0_ids   = platform.device.malloc(mesh.Nfp*minNfp7*sizeof(int),           tmp_L0_ids);

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties
  
  const dfloat p_half = 1./2.;
  int blockMax = 256;
  
  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;
  kernelInfo["defines/" "p_Nfields"]= Nfields;
  kernelInfo["defines/" "p_half"]= p_half;
  kernelInfo["defines/" "p_NfacesNfp"]= mesh.Nfaces*mesh.Nfp;

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;
  
  // set kernel name suffix
  std::string suffix;
  if(mesh.elementType==Mesh::TRIANGLES)
    suffix = "TriBB2D";
  if(mesh.elementType==Mesh::QUADRILATERALS)
    suffix = "QuadBB2D";
  if(mesh.elementType==Mesh::TETRAHEDRA)
    suffix = "TetBB3D";
  if(mesh.elementType==Mesh::HEXAHEDRA)
    suffix = "HexBB3D";

  std::string oklFilePrefix = DBNS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string surfaceFileName, surfaceKernelName;

  occa::streamTag start, stop;

  // kernels from surface file
  surfaceFileName   = oklFilePrefix + "bnsSurface" + suffix + oklFileSuffix;
  surfaceKernelName = "bnsSurface" + suffix;

  // set up flop counts for tets
  long long int NFLOP = ((long long int)mesh.Nelements)*(mesh.Nfp*mesh.Nfaces*(10 + 18 + 4) + mesh.Np*((1+mesh.Nfaces*mesh.Nfp)*Nfields*2));

  long long int NBYTES   = mesh.Nelements*(
					   (mesh.Nsgeo*mesh.Nfaces  +
					    mesh.Nfp*mesh.Nfaces*Nfields*2 +
					    mesh.Np*(Nfields*2) )*sizeof(dfloat)
					   + (mesh.Nfp*mesh.Nfaces*2*sizeof(dlong)));
												       
  int bestNblockS = 0;
  int bestNsur = 0;
  double bestElapsed = 1e9;;
  dfloat T = 0;

  for(int Nsur=1;Nsur<=5;++Nsur){
    for(int NblockS=1;NblockS<=blockMax/maxNodes;++NblockS){

      int LDS = (mesh.Nfp*mesh.Nfaces*Nsur*NblockS*4 + Nsur*NblockS*mesh.dim*mesh.Nfaces)*sizeof(dfloat);
      if(LDS<36*1024){
      properties_t surfaceKernelInfo = kernelInfo;
      
      surfaceKernelInfo["defines/" "p_NblockS"]= NblockS;
      surfaceKernelInfo["defines/" "p_Nsur"]= Nsur;
      surfaceKernelInfo["defines/" "p_NfpNfaces"]= mesh.Nfp*mesh.Nfaces;
      surfaceKernelInfo["defines/" "p_L0_nnz"]= minNfp7; // need to transpose data for L0
      surfaceKernelInfo["defines/" "p_EEL_nnz"]= p_EEL_nnz; // need to transpose data for L0
      
      platform.device.finish();
      
      surfaceKernel =  platform.buildKernel(surfaceFileName, surfaceKernelName,
					    surfaceKernelInfo);
      
      
      int Nwarm = 100;
      for(int w=0;w<Nwarm;++w){
	
	surfaceKernel(mesh.Nelements,
		      mesh.o_sgeo,
		      mesh.o_vmapM,
		      mesh.o_vmapP,
		      o_EEL_ids,
		      o_EEL_vals,
		      o_L0_ids,
		      o_L0_vals,
		      o_Q,
		      o_RHS);
	
      }
      
      platform.device.finish();
      
      start = platform.device.tagStream();
      
      
      int Nrun = 10;
      for(int r=0;r<Nrun;++r){

	surfaceKernel(mesh.Nelements,
		      mesh.o_sgeo,
		      mesh.o_vmapM,
		      mesh.o_vmapP,
		      o_EEL_ids,
		      o_EEL_vals,
		      o_L0_ids,
		      o_L0_vals,
		      o_Q,
		      o_RHS);
      }
      
      stop = platform.device.tagStream();
      
      platform.device.finish();
      
      double elapsed = platform.device.timeBetween(start, stop);
      elapsed /= Nrun;
      
      double GFLOPS = NFLOP/(1.e9*elapsed);
      double GBS = NBYTES/(1.e9*elapsed);
      
      printf("%02d, %02d, %5.4e, %5.4e, %5.4e, %d %%%% SURF-BB, Nsur, NblockS, elapsed, GFLOP/s, GB/s, LDS usage\n",
	     Nsur, NblockS, elapsed, GFLOPS, GBS, LDS);
      
      if(elapsed<bestElapsed){
	bestElapsed = elapsed;
	bestNblockS = NblockS;
	bestNsur = Nsur;
      }
    }
  }
  }
  
  double bestGFLOPS = NFLOP/(1.e9*bestElapsed);
  double bestGBS = NBYTES/(1.e9*bestElapsed);
  
  printf("%02d, %02d, %02d, %5.4e, %5.4e, %5.4e %%%% BB SURF BEST - N, Nsur, NblockS, elapsed, GFLOP/s, GB/s\n",
	 mesh.N, bestNsur, bestNblockS, bestElapsed, bestGFLOPS, bestGBS);
  
}

