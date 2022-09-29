
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

#include "acoustics.hpp"

void acoustics_t::Benchmarks(){

  deviceMemory<dfloat> o_Q   = platform.device.malloc(Nfields*mesh.Nelements*mesh.Np*sizeof(dfloat));
  deviceMemory<dfloat> o_RHS = platform.device.malloc(Nfields*mesh.Nelements*mesh.Np*sizeof(dfloat));

  volumeBenchmark(o_Q, o_RHS);
}


void acoustics_t::volumeBenchmark(deviceMemory<dfloat> &o_Q, deviceMemory<dfloat> &o_RHS){

  
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

  std::string oklFilePrefix = DACOUSTICS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string volumeFileName, volumeKernelName;

  occa::streamTag start, stop;

  // kernels from volume file
  volumeFileName   = oklFilePrefix + "acousticsVolume" + suffix + oklFileSuffix;
  volumeKernelName = "acousticsVolume" + suffix;

  // set up flop counts for tets
  long long int NFLOP = mesh.Nelements*mesh.Np*(4*3*2*(long long int)mesh.Np + 32);
  long long int NBYTES   = mesh.Nelements*(mesh.Nvgeo + mesh.Np*(Nfields*2 + 0*mesh.dim*mesh.Np) )*sizeof(dfloat);

int bestNvol = 0, bestNblockV = 0;
  double bestElapsed = 1e9;;
  
  for(int Nvol=1;Nvol<=7;++Nvol){
    for(int NblockV=1;NblockV<=blockMax/mesh.Np;++NblockV){
      properties_t volumeKernelInfo = kernelInfo;
      
      kernelInfo["defines/" "p_NblockV"]= NblockV;
      kernelInfo["defines/" "p_Nvol"]= Nvol;

      platform.device.finish();
      
      volumeKernel =  platform.buildKernel(volumeFileName, volumeKernelName,
					   volumeKernelInfo);
      

      int Nwarm = 100;
      for(int w=0;w<Nwarm;++w){
	volumeKernel(mesh.Nelements,
		     mesh.o_vgeo,
		     mesh.o_D,
		     o_Q,
		     o_RHS);
      }

      platform.device.finish();
      
      start = platform.device.tagStream();
      

      int Nrun = 10;
      for(int r=0;r<Nrun;++r){
	volumeKernel(mesh.Nelements,
		     mesh.o_vgeo,
		     mesh.o_D,
		     o_Q,		 
		     o_RHS);
      }
      
      stop = platform.device.tagStream();
      
      platform.device.finish();
      
      double elapsed = platform.device.timeBetween(start, stop);
      elapsed /= Nrun;
      
      double GFLOPS = NFLOP/(1.e9*elapsed);
      double GBS = NBYTES/(1.e9*elapsed);
      
      printf("%02d, %02d, %5.4e, %5.4e, %5.4e %%%% Nvol, NblockV, elapsed, GFLOPS, GB/s\n",
	     Nvol, NblockV, elapsed, GFLOPS, GBS);

      if(elapsed<bestElapsed){
	bestElapsed = elapsed;
	bestNvol = Nvol;;
	bestNblockV = NblockV;
      }
    }
  }
  
  double bestGFLOPS = NFLOP/(1.e9*bestElapsed);
  double bestGBS = NBYTES/(1.e9*bestElapsed);
  
  printf("%02d, %02d, %5.4e, %5.4e, %5.4e %%%% BEST - Nvol, NblockV, elapsed, GFLOPS, GB/S\n",
	 bestNvol, bestNblockV, bestElapsed, bestGFLOPS, bestGBS);

}

