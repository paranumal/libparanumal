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

#include "advection.hpp"

int main(int argc, char **argv){


  // start up MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  if(argc!=2)
    LIBP_ABORT(string("Usage: ./advectionMain setupfile"));

  advectionSettings_t settings(comm); //sets default settings
  settings.readSettingsFromFile(argv[1]);
  if (!rank) settings.report();

  // set up occa device
  occa::device device;
  occa::properties props;
  occaDeviceConfig(device, comm, settings, props);

  // set up mesh
  mesh_t& mesh = mesh_t::Setup(device, comm, settings, props);

  // set up linear algebra module
  linAlg_t& linAlg = linAlg_t::Setup(device, settings, props);

  // set up advection solver
  advection_t& advection = advection_t::Setup(mesh, linAlg);

  // run
  advection.Run();

  // close down MPI
  MPI_Finalize();
  return LIBP_SUCCESS;


  // test mass matrix inversion
  dfloat *diagInvMassMatrix = (dfloat*) calloc(mesh->Np*mesh->Nelements, sizeof(dfloat));

  dfloat *cubInterpRowSums = (dfloat*) calloc(mesh->cubNq, sizeof(dfloat));
  for(int a=0;a<mesh->cubNq;++a){
    for(int i=0;i<mesh->Nq;++i){
      cubInterpRowSums[a] += mesh->cubInterp[a*mesh->Nq+i];
    }
  }

#pragma omp parallel for
  for(int e=0;e<mesh->Nelements;++e){

#if 0
    for(int n=0;n<mesh->Np;++n){
      diagInvMassMatrix[n+e*mesh->Np] = 1./mesh->vgeo[n + e*mesh->Np*mesh->Nvgeo + JWID*mesh->Np];
    }
#else
    for(int k=0;k<mesh->Nq;++k){
      for(int j=0;j<mesh->Nq;++j){
	for(int i=0;i<mesh->Nq;++i){
	  int id = i + j*mesh->Nq + k*mesh->Nq*mesh->Nq;
	  dfloat res = 0;

	  for(int c=0;c<mesh->cubNq;++c){
	    for(int b=0;b<mesh->cubNq;++b){
	      for(int a=0;a<mesh->cubNq;++a){
		hlong cid = a + b*mesh->cubNq + c*mesh->cubNq*mesh->cubNq;
		dfloat JW = mesh->cubvgeo[mesh->Nvgeo*mesh->cubNp*e + cid + JWID*mesh->cubNp];

		res +=
		  JW*mesh->cubInterp[i+mesh->Nq*a]*cubInterpRowSums[a]
		  *mesh->cubInterp[j+mesh->Nq*b]*cubInterpRowSums[b]
		  *mesh->cubInterp[k+mesh->Nq*c]*cubInterpRowSums[c];
	      }
	    }
	  }
	  diagInvMassMatrix[id + e*mesh->Np] = 1./res;
	}
      }
    }
#endif
  }


  int *iterations = (int*) calloc(mesh->Nelements, sizeof(int));
  dfloat *residuals = (dfloat*) calloc(mesh->Nelements, sizeof(dfloat));
  advection->o_diagInvMassMatrix = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), diagInvMassMatrix);

  mesh->device.finish();


  MPI_Barrier(MPI_COMM_WORLD);

  mesh->device.finish();


  // run
  advectionRun(advection, newOptions);

  mesh->device.finish();

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
