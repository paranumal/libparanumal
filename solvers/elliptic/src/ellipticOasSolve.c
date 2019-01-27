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

#include "elliptic.h"

void ellipticOasSolve(elliptic_t *elliptic, dfloat lambda,
		      occa::memory &o_r, occa::memory &o_z){

  // TW: do not interleave these tasks yet
  
  // 1. restrict to coarse grid
  precon_t *precon = elliptic->precon;
  mesh_t *mesh = elliptic->mesh;

  elliptic_t *elliptic1 = (elliptic_t*) precon->ellipticOneRing; // should rename
  mesh_t *mesh1 = elliptic1->mesh;

  // hack to zero initial guess
  dfloat *h_x = (dfloat*) calloc(mesh1->Np*mesh1->Nelements, sizeof(dfloat));
  elliptic1->o_x.copyFrom(h_x);
  //  elliptic1->o_r.copyFrom(h_x);

  // TW: possibility these device have difference queues
  mesh1->device.finish();
  mesh->device.finish();

  // TW: do I need to mask here
  if (elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_r);
  
  // 3. collect patch rhs    
  ellipticOneRingExchange(elliptic, elliptic1, mesh1->Np*sizeof(dfloat), o_r, elliptic1->o_r);

  // TW: possibility these device have difference queues
  mesh1->device.finish();
  mesh->device.finish();

  // TW: what tolerance to use ?
  dfloat tol = 1e-1;
  
  // patch solve
  if(mesh->rank==0) printf("Starting extended partition iterations:\n");

  ellipticSolve(elliptic1, lambda, tol, elliptic1->o_r, elliptic1->o_x); // may need to zero o_x

  // sum up patches
  ogsGatherScatter(elliptic1->o_x, ogsDfloat, ogsAdd, elliptic->precon->oasOgs);

  // do we need to scale by 1/overlapDegree ?
  
  // just retain core [ actually need to gs all the element contributions]
  elliptic->dotMultiplyKernel(mesh->Nelements*mesh->Np, elliptic->precon->oasOgs->o_invDegree, elliptic1->o_x, o_z);

#if 0
  // 2. solve coarse problem
  //   a. call solver
  //   b. prolongate (watch out for +=)

  // TW: QUESTIONABLE FROM HERE ---->
  MPI_Barrier(mesh->comm);
  if(mesh->rank==0) printf("Starting coarse grid iterations:\n");
  
  elliptic_t *ellipticOasCoarse = (elliptic_t*) (precon->ellipticOasCoarse);
  mesh_t   *meshCoarse   = ellipticOasCoarse->mesh;

  mesh1->device.finish();
  mesh->device.finish();
  meshCoarse->device.finish();
  
  precon->oasRestrictionKernel(meshCoarse->Nelements,
			       precon->o_oasRestrictionMatrix,
			       o_r, ellipticOasCoarse->o_r);

  mesh1->device.finish();
  mesh->device.finish();
  meshCoarse->device.finish();

  if (ellipticOasCoarse->Nmasked) meshCoarse->maskKernel(ellipticOasCoarse->Nmasked, ellipticOasCoarse->o_maskIds,
							 ellipticOasCoarse->o_r);
  
  ogsGatherScatter(ellipticOasCoarse->o_r, ogsDfloat, ogsAdd, ellipticOasCoarse->ogs);
  
  ellipticOasCoarse->dotMultiplyKernel(meshCoarse->Nelements*meshCoarse->Np,
				       meshCoarse->ogs->o_invDegree,
				       ellipticOasCoarse->o_r,
				       ellipticOasCoarse->o_r);

  dfloat *h_xCoarse = (dfloat*) calloc(meshCoarse->Np*meshCoarse->Nelements, sizeof(dfloat));
  ellipticOasCoarse->o_x.copyFrom(h_xCoarse);
  free(h_xCoarse);
  
  ellipticSolve(ellipticOasCoarse, lambda, tol, ellipticOasCoarse->o_r, ellipticOasCoarse->o_x);
  
  mesh1->device.finish();
  mesh->device.finish();
  meshCoarse->device.finish();
  
  // prolongate to QN (note kernel expects restriction matrix)
  // do we need to weight the sum against patches?
  precon->oasProlongationKernel(mesh->Nelements, precon->o_oasRestrictionMatrix,
				ellipticOasCoarse->o_x, o_z);

  mesh1->device.finish();
  mesh->device.finish();
  meshCoarse->device.finish();

  MPI_Barrier(mesh->comm);
  if(mesh->rank==0) printf("Ending coarse grid iterations:\n Outer ");  
  // TW: QUESTIONABLE TO HERE <----
#endif
  
  // TW: is this needed ?
  if (elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_z);
  
  free(h_x);
}
