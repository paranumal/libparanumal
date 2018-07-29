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

void ellipticRunBenchmark3D(solver_t *solver, char *options, occa::kernelInfo kernelInfo, char *kernelFileName, int Nblocks, int Nnodes) {

  mesh3D *mesh = solver->mesh;
	
  int Ntrials = 4;
	
  size_t L1CacheSize = 24576; //L1 cache size of test device in bytes (24KB for P100)
	
  int NKernels;
  char kernelName[BUFSIZ];

  if (strstr(kernelFileName,"ellipticBRGradientVolume3D.okl")) {
    NKernels = 1;
    sprintf(kernelName, "ellipticBRGradientVolume3D");
  } else if (strstr(kernelFileName,"ellipticBRGradientSurface3D.okl")) {
    NKernels = 1;
    sprintf(kernelName, "ellipticBRGradientSurface3D");
  } else if (strstr(kernelFileName,"ellipticBRDivergenceVolume3D.okl")) {
    NKernels = 1;
    sprintf(kernelName, "ellipticBRDivergenceVolume3D");
  } else if (strstr(kernelFileName,"ellipticBRDivergenceSurface3D.okl")) {
    NKernels = 1;
    sprintf(kernelName, "ellipticBRBBDivergenceSurface3D");
  } else if (strstr(kernelFileName,"ellipticBRBBGradientVolume3D.okl")) {
    NKernels = 5;
    sprintf(kernelName, "ellipticBRBBGradientVolume3D");
  } else if (strstr(kernelFileName,"ellipticBRBBGradientPackedVolume3D.okl")) {
    NKernels = 2;
    sprintf(kernelName, "ellipticBRBBGradientPackedVolume3D");

  } else if (strstr(kernelFileName,"ellipticBRBBGradientSurface3D.okl")) {
    NKernels = 1;
    sprintf(kernelName, "ellipticBRBBGradientSurface3D");
  } else if (strstr(kernelFileName,"ellipticBRBBDivergenceVolume3D.okl")) {
    NKernels = 1;
    sprintf(kernelName, "ellipticBRBBDivergenceVolume3D");
  } else if (strstr(kernelFileName,"ellipticBRBBDivergenceSurface3D.okl")) {
    NKernels = 1;
    sprintf(kernelName, "ellipticBRBBDivergenceSurface3D");
  }
		
  dfloat time = 0.;
	
  char testkernelName[BUFSIZ];
  occa::kernel testKernel;
  for(int i=0; i<NKernels; i++) {
	
    sprintf(testkernelName, "%s_v%d", kernelName,  i);
    printf("%s Kernel #%03d\n", kernelFileName, i);
			
    testKernel = mesh->device.buildKernelFromSource(kernelFileName,testkernelName,kernelInfo);
			
    // sync processes
    mesh->device.finish();
    occa::streamTag start = mesh->device.tagStream();
    for(int it=0;it<Ntrials;++it) {
      if (strstr(kernelFileName,"ellipticBRGradientVolume3D.okl")) {
	testKernel(mesh->Nelements,
		   mesh->o_vgeo,
		   mesh->o_DrT,
		   mesh->o_DsT,
		   mesh->o_DtT,
		   solver->o_p,
		   solver->o_grad);
      } else if (strstr(kernelFileName,"ellipticBRGradientSurface3D.okl")) {
	testKernel(mesh->Nelements,
		   mesh->o_vmapM,
		   mesh->o_vmapP,
		   mesh->o_sgeo,
		   mesh->o_EToB,
		   mesh->o_LIFTT,
		   solver->o_p,
		   solver->o_grad);
      } else if (strstr(kernelFileName,"ellipticBRDivergenceVolume3D.okl")) {
	testKernel(mesh->Nelements,
		   mesh->o_vgeo,
		   mesh->o_DrT,
		   mesh->o_DsT,
		   mesh->o_DtT,
		   solver->o_grad,
		   solver->o_Ap);
      } else if (strstr(kernelFileName,"ellipticBRDivergenceSurface3D.okl")) {
	testKernel(mesh->Nelements,
		   mesh->o_vmapM,
		   mesh->o_vmapP,
		   solver->lambda,
		   solver->tau,
		   mesh->o_vgeo,
		   mesh->o_sgeo,
		   mesh->o_EToB,
		   mesh->o_LIFTT,
		   mesh->o_MM,
		   solver->o_p,
		   solver->o_grad,
		   solver->o_Ap);
      } else if (strstr(kernelFileName,"ellipticBRBBGradientVolume3D.okl")) {
	testKernel(mesh->Nelements,
		   mesh->o_vgeo,
		   mesh->o_D0ids,
		   mesh->o_D1ids,
		   mesh->o_D2ids,
		   mesh->o_D3ids,
		   mesh->o_Dvals,
		   solver->o_p,
		   solver->o_grad);
      } else if (strstr(kernelFileName,"ellipticBRBBGradientPackedVolume3D.okl")) {
	testKernel(mesh->Nelements,
		   mesh->o_vgeo,
		   mesh->o_D0ids,
		   mesh->o_packedDids,
		   mesh->o_Dvals,
		   solver->o_p,
		   solver->o_grad);
      } else if (strstr(kernelFileName,"ellipticBRBBGradientSurface3D.okl")) {
	testKernel(mesh->Nelements,
		   mesh->o_vmapM,
		   mesh->o_vmapP,
		   mesh->o_sgeo,
		   mesh->o_EToB,
		   mesh->o_L0ids,
		   mesh->o_L0vals,
		   mesh->o_ELids,
		   mesh->o_ELvals,
		   solver->o_p,
		   solver->o_grad);
      } else if (strstr(kernelFileName,"ellipticBRBBDivergenceVolume3D.okl")) {
	testKernel(mesh->Nelements,
		   mesh->o_vgeo,
		   mesh->o_D0ids,
		   mesh->o_D1ids,
		   mesh->o_D2ids,
		   mesh->o_D3ids,
		   mesh->o_Dvals,
		   solver->o_grad,
		   solver->o_Ap);
      } else if (strstr(kernelFileName,"ellipticBRBBDivergenceSurface3D.okl")) {
	testKernel(mesh->Nelements,
		   mesh->o_vmapM,
		   mesh->o_vmapP,
		   solver->lambda,
		   solver->tau,
		   mesh->o_vgeo,
		   mesh->o_sgeo,
		   mesh->o_EToB,
		   mesh->o_L0ids,
		   mesh->o_L0vals,
		   mesh->o_ELids,
		   mesh->o_ELvals,
		   mesh->o_MM,
		   solver->o_p,
		   solver->o_grad,
		   solver->o_Ap);
      }
    }
    occa::streamTag end = mesh->device.tagStream();
    mesh->device.finish();
    double kernelElapsed = mesh->device.timeBetween(start,end);
			
    size_t Nflops;          //flops per element performed
    size_t NbytesGlobal;    //bytes of global memory R+W per element
    size_t NbytesShared;    //bytes of shared memory R+W per element
    long long int NbytesCacheMiss; //bytes of streamed operators which overflow L1 cache
			
			
    if (strstr(kernelFileName,"ellipticBRGradientVolume3D.okl")) {
      NbytesGlobal  = sizeof(dfloat)*1*mesh->Np;   //volume data
      NbytesGlobal += sizeof(dfloat)*9;            //vgeo factors
      NbytesGlobal += sizeof(dfloat)*3*mesh->Np;   //write grad
					
      NbytesCacheMiss  = sizeof(dfloat)*3*mesh->Np*mesh->Np; //DrT, DsT, and DtT
      NbytesCacheMiss -= L1CacheSize; //L1 cache size
					
      NbytesShared  = sizeof(dfloat)*1*mesh->Np;             //load of volume data to shared
      NbytesShared += sizeof(dfloat)*3*mesh->Np*mesh->Np;   //Dr Ds Dt matvec
					
      Nflops  = 2*3*mesh->Np*mesh->Np;  //Dr Ds matvec
      Nflops += mesh->Np*3*5;          //vgeo factors
    } else if (strstr(kernelFileName,"ellipticBRGradientSurface3D.okl")) {
      NbytesGlobal  = sizeof(dfloat)*2*mesh->Nfaces*mesh->Nfp;  //trace data
      NbytesGlobal += sizeof(dfloat)*5*mesh->Nfaces;            //sgeo factors
      NbytesGlobal += sizeof(dfloat)*(3+3)*mesh->Np;            //read + write rhs
      NbytesGlobal += sizeof(int)*2*mesh->Nfaces*mesh->Nfp;    //vmapM and vmapP
      NbytesGlobal += sizeof(int)*mesh->Nfaces;                //EToB flag
					
      NbytesCacheMiss  = sizeof(dfloat)*mesh->Np*mesh->Nfaces*mesh->Nfp; //LIFT
      NbytesCacheMiss -= L1CacheSize; //L1 cache size
					
      NbytesShared  = sizeof(dfloat)*3*mesh->Nfp*mesh->Nfaces;            //store of fluxes
      NbytesShared += sizeof(dfloat)*3*mesh->Np*mesh->Nfp*mesh->Nfaces;   //LIFT matvec
					
      Nflops  = mesh->Nfp*mesh->Nfaces*(2+3*3);       //eval fluxes
      Nflops += mesh->Np*mesh->Nfp*mesh->Nfaces*2*3;  // LIFT matvec
    } else if (strstr(kernelFileName,"ellipticBRDivergenceVolume3D.okl")) {
      NbytesGlobal  = sizeof(dfloat)*3*mesh->Np;   //volume data
      NbytesGlobal += sizeof(dfloat)*9;            //vgeo factors
      NbytesGlobal += sizeof(dfloat)*1*mesh->Np;   //write Aq
					
      NbytesCacheMiss  = sizeof(dfloat)*3*mesh->Np*mesh->Np; //DrT, DsT, and DtT
      NbytesCacheMiss -= L1CacheSize; //L1 cache size
					
      NbytesShared  = sizeof(dfloat)*3*mesh->Np;             //load of volume data to shared
      NbytesShared += sizeof(dfloat)*9*mesh->Np*mesh->Np;   //Dr Ds matvec
					
      Nflops  = 2*9*mesh->Np*mesh->Np;  //Dr Ds matvec
      Nflops += mesh->Np*17;             //vgeo factors
    } else if (strstr(kernelFileName,"ellipticBRDivergenceSurface3D.okl")) {
      NbytesGlobal  = sizeof(dfloat)*8*mesh->Nfaces*mesh->Nfp;  //trace data
      NbytesGlobal += sizeof(dfloat)*5*mesh->Nfaces;            //sgeo factors
      NbytesGlobal += sizeof(dfloat)*(2+1)*mesh->Np;            //read q + Aq and write Aq
      NbytesGlobal += sizeof(int)*2*mesh->Nfaces*mesh->Nfp;    //vmapM and vmapP
      NbytesGlobal += sizeof(int)*mesh->Nfaces;                //EToB flag
					
      NbytesCacheMiss  = sizeof(dfloat)*mesh->Np*mesh->Nfaces*mesh->Nfp; //LIFT
      NbytesCacheMiss += sizeof(dfloat)*mesh->Np*mesh->Np; // MM
      NbytesCacheMiss -= L1CacheSize; //L1 cache size
					
      NbytesShared  = sizeof(dfloat)*1*mesh->Nfp*mesh->Nfaces;  //store of fluxes
      NbytesShared += sizeof(dfloat)*1*mesh->Np*mesh->Nfp*mesh->Nfaces;   //LIFT matvec
      NbytesShared += sizeof(dfloat)*1*mesh->Np;                //store of Aq
      NbytesShared += sizeof(dfloat)*1*mesh->Np*mesh->Np;       //MM matvec
					
      Nflops  = mesh->Nfp*mesh->Nfaces*(2+12);       //eval fluxes
      Nflops += mesh->Np*mesh->Nfp*mesh->Nfaces*2*1;  // LIFT matvec
      Nflops += mesh->Np*2;  // store Aq + lambda*q
      Nflops += mesh->Np*mesh->Np*2;  // MM matvec
    } else if (strstr(kernelFileName,"ellipticBRBBGradientVolume3D.okl")) {
      NbytesGlobal  = sizeof(dfloat)*1*mesh->Np;   //volume data
      NbytesGlobal += sizeof(dfloat)*9;            //vgeo factors
      NbytesGlobal += sizeof(dfloat)*3*mesh->Np;   //write grad
					
      NbytesCacheMiss  = sizeof(dfloat)*4*mesh->Np; //Dvals
      NbytesCacheMiss += sizeof(int)*4*4*mesh->Np; //D0ids D1ids D2ids D3ids
      NbytesCacheMiss -= L1CacheSize; //L1 cache size
					
      NbytesShared  = sizeof(dfloat)*1*mesh->Np;       //load of volume data to shared
      NbytesShared += sizeof(dfloat)*3*2*4*mesh->Np;   //D0ids D1ids D2ids D3ids
					
      Nflops  = 3*12*mesh->Np;  //D0ids D1ids D2ids D3ids
      Nflops += mesh->Np*3*5;   //vgeo factors
    } else if (strstr(kernelFileName,"ellipticBRBBGradientPackedVolume3D.okl")) {
      NbytesGlobal  = sizeof(dfloat)*1*mesh->Np;   //volume data
      NbytesGlobal += sizeof(dfloat)*9;            //vgeo factors
      NbytesGlobal += sizeof(dfloat)*3*mesh->Np;   //write grad

      // TW need to revisit
      NbytesCacheMiss  = sizeof(dfloat)*4*mesh->Np; //Dvals
      NbytesCacheMiss += sizeof(int)*4*mesh->Np;   //D0ids
      NbytesCacheMiss += sizeof(char)*4*3*mesh->Np;   // D1ids D2ids D3ids
      NbytesCacheMiss -= L1CacheSize; //L1 cache size
      NbytesCacheMiss = 0;

      // TW: this is wrong (need to add all shared loads and stores)
      NbytesShared  = sizeof(dfloat)*1*mesh->Np;       //load of volume data to shared
      NbytesShared += sizeof(int)*4*mesh->Np;   //D0ids
      NbytesShared += sizeof(dfloat)*4*mesh->Np;   //Dvals
      NbytesShared += sizeof(char)*4*3*mesh->Np;   //D0ids D1ids D2ids D3ids
					
      Nflops  = (7+9*3)*mesh->Np;  //D0ids D1ids D2ids D3ids
      Nflops += mesh->Np*3*5;   //vgeo factors

    } else if (strstr(kernelFileName,"ellipticBRBBGradientSurface3D.okl")) {
      NbytesGlobal  = sizeof(dfloat)*2*mesh->Nfaces*mesh->Nfp;  //trace data
      NbytesGlobal += sizeof(dfloat)*5*mesh->Nfaces;            //sgeo factors
      NbytesGlobal += sizeof(dfloat)*(3+3)*mesh->Np;            //read + write rhs
      NbytesGlobal += sizeof(int)*2*mesh->Nfaces*mesh->Nfp;    //vmapM and vmapP
      NbytesGlobal += sizeof(int)*mesh->Nfaces;                //EToB flag
					
      NbytesCacheMiss  = sizeof(dfloat)*mesh->Nfp*7;               //L0vals
      NbytesCacheMiss += sizeof(int)*mesh->Nfp*7;                 //L0ids
      NbytesCacheMiss += sizeof(dfloat)*mesh->Np*mesh->max_EL_nnz; //ELvals
      NbytesCacheMiss += sizeof(int)*mesh->Np*mesh->max_EL_nnz;   //ELids
      NbytesCacheMiss -= L1CacheSize; //L1 cache size
					
      NbytesShared  = sizeof(dfloat)*3*mesh->Nfp*mesh->Nfaces;     //store of fluxes
      NbytesShared += sizeof(dfloat)*3*7*mesh->Nfaces*mesh->Nfp;   //L0
      NbytesShared += sizeof(dfloat)*3*mesh->Nfaces*mesh->Nfp;     //store L0flux
      NbytesShared += sizeof(dfloat)*3*mesh->Np*mesh->max_EL_nnz;  //EL
					
      Nflops  = mesh->Nfp*mesh->Nfaces*(2+3*3); //eval fluxes
      Nflops += mesh->Nfp*mesh->Nfaces*3*2*7;     //L0
      Nflops += mesh->Np*3*2*mesh->max_EL_nnz;    //EL
    } else if (strstr(kernelFileName,"ellipticBRBBDivergenceVolume3D.okl")) {
      NbytesGlobal  = sizeof(dfloat)*3*mesh->Np;   //volume data
      NbytesGlobal += sizeof(dfloat)*9;            //vgeo factors
      NbytesGlobal += sizeof(dfloat)*1*mesh->Np;   //write Aq
					
      NbytesCacheMiss  = sizeof(dfloat)*4*mesh->Np; //Dvals
      NbytesCacheMiss += sizeof(int)*4*4*mesh->Np; //D0ids D1ids D2ids D3ids
      NbytesCacheMiss -= L1CacheSize; //L1 cache size
					
      NbytesShared  = sizeof(dfloat)*3*mesh->Np;       //load of volume data to shared
      NbytesShared += sizeof(dfloat)*3*3*2*4*mesh->Np;   //D0ids D1ids D2ids D3ids
					
      Nflops  = 3*3*12*mesh->Np;  //D0ids D1ids D2ids D3ids
      Nflops += mesh->Np*17;             //vgeo factors
    } else if (strstr(kernelFileName,"ellipticBRBBDivergenceSurface3D.okl")) {
      NbytesGlobal  = sizeof(dfloat)*8*mesh->Nfaces*mesh->Nfp;  //trace data
      NbytesGlobal += sizeof(dfloat)*5*mesh->Nfaces;            //sgeo factors
      NbytesGlobal += sizeof(dfloat)*(2+1)*mesh->Np;            //read q + Aq and write Aq
      NbytesGlobal += sizeof(int)*2*mesh->Nfaces*mesh->Nfp;    //vmapM and vmapP
      NbytesGlobal += sizeof(int)*mesh->Nfaces;                //EToB flag
					
      NbytesCacheMiss  = sizeof(dfloat)*mesh->Nfp*7;               //L0vals
      NbytesCacheMiss += sizeof(int)*mesh->Nfp*7;                 //L0ids
      NbytesCacheMiss += sizeof(dfloat)*mesh->Np*mesh->max_EL_nnz; //ELvals
      NbytesCacheMiss += sizeof(int)*mesh->Np*mesh->max_EL_nnz;   //ELids
      NbytesCacheMiss += sizeof(dfloat)*mesh->Np*mesh->Np; // MM
      NbytesCacheMiss -= L1CacheSize; //L1 cache size
					
      NbytesShared  = sizeof(dfloat)*1*mesh->Nfp*mesh->Nfaces;  //store of fluxes
      NbytesShared += sizeof(dfloat)*1*7*mesh->Nfaces*mesh->Nfp;   //L0
      NbytesShared += sizeof(dfloat)*1*mesh->Nfaces*mesh->Nfp;     //store L0flux
      NbytesShared += sizeof(dfloat)*1*mesh->Np*mesh->max_EL_nnz;  //EL
      NbytesShared += sizeof(dfloat)*1*mesh->Np;                //store of Aq
      NbytesShared += sizeof(dfloat)*1*mesh->Np*mesh->Np;       //MM matvec
					
      Nflops  = mesh->Nfp*mesh->Nfaces*(2+12);       //eval fluxes
      Nflops += mesh->Nfp*mesh->Nfaces*2*7;     //L0
      Nflops += mesh->Np*2*mesh->max_EL_nnz;    //EL
      Nflops += mesh->Np*2;  // store Aq + lambda*q
      Nflops += mesh->Np*mesh->Np*2;  // MM matvec
    }
				
    //global memory bandwidth benchmark
    NbytesGlobal /= 2; //bi-directional memory bus
    occa::memory o_foo = mesh->device.malloc(NbytesGlobal*mesh->Nelements);
    occa::memory o_bah = mesh->device.malloc(NbytesGlobal*mesh->Nelements);
			
    mesh->device.finish();
    occa::streamTag startCopy = mesh->device.tagStream();
    for(int it=0;it<Ntrials;++it) {
      o_bah.copyTo(o_foo);
    }
    occa::streamTag endCopy = mesh->device.tagStream();
    mesh->device.finish();
    double copyElapsed = mesh->device.timeBetween(startCopy, endCopy);
			
    if (NbytesCacheMiss > 0) {
      //caching bandwidth benchmark (estimate with 1000 elements)
      NbytesCacheMiss /= 2; //bi-directional memory bus
      occa::memory o_fooCache = mesh->device.malloc(NbytesCacheMiss*1000);
      occa::memory o_bahCache = mesh->device.malloc(NbytesCacheMiss*1000);
					
      mesh->device.finish();
      startCopy = mesh->device.tagStream();
      for(int it=0;it<Ntrials;++it) {
	o_bahCache.copyTo(o_fooCache);
      }
      endCopy = mesh->device.tagStream();
      mesh->device.finish();
      //copyElapsed += mesh->device.timeBetween(startCopy, endCopy)*(((dfloat) mesh->Nelements)/(Nnodes*1000));
    }
				
    // Compute Data
    double copyBandwidth   = mesh->Nelements*((NbytesGlobal*Ntrials*2)/(1e9*copyElapsed));
    double kernelBandwidth = mesh->Nelements*((NbytesGlobal*Ntrials*2)/(1e9*kernelElapsed));
			
    double copyGFLOPS   = mesh->Nelements*Nflops*Ntrials/(1e9*copyElapsed);
    double kernelGFLOPS = mesh->Nelements*Nflops*Ntrials/(1e9*kernelElapsed);
			
    double shmemBound = 7882*Nflops/( (double) NbytesShared);
    double intensity  = kernelGFLOPS*Ntrials/kernelBandwidth;
    double roofline   = mymin(copyGFLOPS, shmemBound);
			
    printf("%s\t %d \t %d \t %6.4E\t %6.4E\t %6.4E\t %6.4E\t %6.4E\t %d \t %d\t %6.4E "
	   "[ N K KernelTime CopyTime (Kernel GFLOPS/s) (Copy GFLOPS/s) (KernelTime/copyTime) Nblocks Nnodes (GNODES/s)\n ]",
	   testkernelName, mesh->N, mesh->Nelements, kernelElapsed/Ntrials, copyElapsed/Ntrials, kernelGFLOPS, copyGFLOPS, kernelElapsed/copyElapsed, Nblocks, Nnodes,
	   mesh->Nelements*Ntrials*mesh->Np/(1e9*kernelElapsed));
  }
}
