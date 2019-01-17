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

#include "advection.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2){
    printf("usage2: ./advectionMain setupfile\n");
    exit(-1);
  }

  // if argv > 2 then should load input data from argv
  setupAide newOptions(argv[1]);
  
  // set up mesh stuff
  string fileName;
  int N, dim, elementType;

  newOptions.getArgs("MESH FILE", fileName);
  newOptions.getArgs("POLYNOMIAL DEGREE", N);
  newOptions.getArgs("ELEMENT TYPE", elementType);
  newOptions.getArgs("MESH DIMENSION", dim);

  dfloat tol = 1e-8;

  newOptions.getArgs("MASS MATRIX TOLERANCE", tol);

  
  // set up mesh
  mesh_t *mesh;
  switch(elementType){
  case TRIANGLES:
  case TETRAHEDRA:
    printf("Triangles and tetrahedra are not currently supported for this code, exiting ...\n");
    exit(-1);
  case QUADRILATERALS:
    mesh = meshSetupQuad2D((char*)fileName.c_str(), N); break;
  case HEXAHEDRA:
    mesh = meshSetupHex3D((char*)fileName.c_str(), N); break;
  }

  if(elementType==HEXAHEDRA){
    
    /* rescale to unit box and transform */
    hlong allNelements = mesh->Nelements+mesh->totalHaloPairs;
    for(int n=0;n<allNelements*mesh->Np;++n){
      mesh->x[n] = 0.5*(mesh->x[n]+1);
      mesh->y[n] = 0.5*(mesh->y[n]+1);
      mesh->z[n] = 0.5*(mesh->z[n]+1);
    }
    
    // compute geometric factors
    meshGeometricFactorsHex3D(mesh);
    meshSurfaceGeometricFactorsHex3D(mesh);
  }
  
  char *boundaryHeaderFileName; // could sprintf
  if(dim==2)
    boundaryHeaderFileName = strdup(DADVECTION "/advectionBox2D.h"); // default
  if(dim==3)
    boundaryHeaderFileName = strdup(DADVECTION "/advectionBox3D.h"); // default

  // set up advection stuff
  advection_t *advection = advectionSetup(mesh, newOptions, boundaryHeaderFileName);

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

#if 0

  occa::memory o_qnew = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), iterations);
  occa::memory o_iterations = mesh->device.malloc(mesh->Nelements*sizeof(int), iterations);
  occa::memory o_residuals  = mesh->device.malloc(mesh->Nelements*sizeof(dfloat), residuals);
  for(int n=0;n<mesh->Np*mesh->Nelements;++n){
    advection->rhsq[n] = sin(M_PI*mesh->x[n])*sin(M_PI*mesh->y[n])*sin(M_PI*mesh->z[n]); // drand48();
  }
  advection->o_rhsq.copyFrom(advection->rhsq);

  int maxIterations = 3;

  mesh->device.finish();

#if 0
  for(int test=0;test<10;++test)
  advection->invertMassMatrixKernel(mesh->Nelements,
				    tol,
				    maxIterations,
				    mesh->o_cubvgeo,
				    mesh->o_cubInterpT,
				    advection->o_diagInvMassMatrix,
				    advection->o_rhsq,
				    advection->o_q,
				    o_iterations);
#else
  for(int test=0;test<10;++test)
    advection->invertMassMatrixKernel(mesh->NinternalElements,
				      mesh->o_internalElementIds,
				      mesh->dt,
				      mesh->rka[0],
				      mesh->rkb[0],
				      tol,
				      maxIterations,
				      mesh->o_cubvgeo,
				      mesh->o_cubInterpT,
				      advection->o_diagInvMassMatrix,
				      advection->o_rhsq,
				      advection->o_resq,
				      advection->o_q,
				      o_qnew,
				      o_iterations);
#endif
				    
  mesh->device.finish();
  
  o_iterations.copyTo(iterations);
  o_residuals.copyTo(residuals);

  for(hlong e=0;e<50;++e){
    printf("element %d took %d iterations to reach residual^2 of %17.15lg\n",
	   e, iterations[e], residuals[e]);
  }
  mesh->device.finish();
  MPI_Finalize();
  exit(-1);
#endif
  
  mesh->device.finish();

  
  MPI_Barrier(MPI_COMM_WORLD);

  mesh->device.finish();

#if 0    
  occa::streamTag start = mesh->device.tagStream();

  int rk = 0, Nsteps = 50;

  for(int test=0;test<Nsteps;++test)
    advection->invertMassMatrixCombinedKernel(mesh->NinternalElements,
					      mesh->o_internalElementIds,
					      mesh->dt,
					      mesh->rka[rk], 
					      mesh->rkb[rk], 
					      mesh->o_vgeo, 
					      mesh->o_Dmatrices,
					      advection->o_advectionVelocityJW,
					      mesh->o_vmapM,
					      mesh->o_vmapP,
					      advection->o_advectionVelocityM,
					      advection->o_advectionVelocityP,
					      mesh->o_cubvgeo,
					      mesh->o_cubInterpT,
					      advection->o_diagInvMassMatrix,
					      tol,
					      maxIterations,
					      advection->o_resq,
					      advection->o_q,
					      advection->o_qtmp0);
  
  occa::streamTag end = mesh->device.tagStream();
  
  mesh->device.finish();

  double elapsed = mesh->device.timeBetween(start, end);

  printf("%d %d %d %lg %lg %lg %lg %lg \%\%[ MMDG: N, Nel, Nodes, elapsed, time/step, nodes/time, gnodes*Nsteps/time, gnodes*Nstages*Nsteps/time, %s ]\n",
	 mesh->N,
	 mesh->Nelements,
	 mesh->Np*mesh->Nelements,
	 elapsed,
	 elapsed/Nsteps,
	 (1.*mesh->Np*mesh->Nelements)/(elapsed),
	 (1.*mesh->Np*mesh->Nelements*Nsteps)/(1.e9*elapsed),
	 (1.*mesh->Np*mesh->Nelements*Nsteps)/(1.e9*elapsed),
	 newOptions.getArgs("ADVECTION FORMULATION").c_str()
	 );

#endif

  // run
  advectionRun(advection, newOptions);

  mesh->device.finish();
  
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
