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
  free(h_x);

  // TW: possibility these device have difference queues
  mesh1->device.finish();
  mesh->device.finish();

  // 3. collect patch rhs    
  ellipticOneRingExchange(elliptic, elliptic1, mesh1->Np*sizeof(dfloat), o_r, elliptic1->o_r);

  // TW: possibility these device have difference queues
  mesh1->device.finish();
  mesh->device.finish();

  // TW: what tolerance to use ?
  dfloat tol1 = 1e-2;
  
  // patch solve
  if(mesh->rank==0) printf("Starting extended partition iterations:\n");
  ellipticSolve(elliptic1, lambda, tol1, elliptic1->o_r, elliptic1->o_x); // may need to zero o_x

  // will gather over all patches - so have to remove local multiplicity
  elliptic1->dotMultiplyKernel(mesh1->Nelements*mesh1->Np, elliptic1->ogs->o_invDegree, elliptic1->o_x, elliptic1->o_z);
  
  // sum up overlapping patches
  ogsGatherScatter(elliptic1->o_z, ogsDfloat, ogsAdd, elliptic->precon->oasOgs);

  o_z.copyFrom(elliptic1->o_z, mesh->Nelements*mesh->Np*sizeof(dfloat), 0);
  
  // 2. solve coarse problem
  //   a. call solver
  //   b. prolongate (watch out for +=)

  // TW: QUESTIONABLE FROM HERE ---->
  if(mesh->rank==0) printf("Starting coarse grid iterations:\n");
  
  elliptic_t *ellipticOasCoarse = (elliptic_t*) (precon->ellipticOasCoarse);
  mesh_t   *meshCoarse   = ellipticOasCoarse->mesh;

  mesh1->device.finish();
  mesh->device.finish();
  meshCoarse->device.finish();
  
  occa::memory o_tmpCoarse = mesh->device.malloc(meshCoarse->Np*meshCoarse->Nelements*sizeof(dfloat));
  
#if 1
  precon->oasRestrictionKernel(meshCoarse->Nelements,
			       precon->o_oasRestrictionMatrix,
			       o_r, ellipticOasCoarse->o_r);
#else
  // do this on mesh device
  precon->oasRestrictionKernel(meshCoarse->Nelements,
			       precon->o_oasRestrictionMatrix,
			       o_r, o_tmpCoarse);

  // copy from mesh device to coarse device
  o_tmpCoarse.copyTo(ellipticOasCoarse->o_r.ptr());
#endif
  
  mesh1->device.finish();
  mesh->device.finish();
  meshCoarse->device.finish();
  
  // why do I Have to do (1/deg)*S*G*o_rCoarse here ? ---------->
  ogsGatherScatter(ellipticOasCoarse->o_r, ogsDfloat, ogsAdd, ellipticOasCoarse->ogs);

  ellipticOasCoarse->dotMultiplyKernel(meshCoarse->Nelements*meshCoarse->Np,
				       meshCoarse->ogs->o_invDegree,
				       ellipticOasCoarse->o_r,
				       ellipticOasCoarse->o_r);
  // <----------
  
  dfloat *h_xCoarse = (dfloat*) calloc(meshCoarse->Np*meshCoarse->Nelements, sizeof(dfloat));
  ellipticOasCoarse->o_x.copyFrom(h_xCoarse);
  free(h_xCoarse);

  dfloat tolCoarse = 1e-2;
  
  ellipticSolve(ellipticOasCoarse, lambda, tolCoarse, ellipticOasCoarse->o_r, ellipticOasCoarse->o_x);
  
  mesh1->device.finish();
  mesh->device.finish();
  meshCoarse->device.finish();

  // prolongate to QN (note kernel expects restriction matrix)
  // do we need to weight the sum against patches?
#if 1
  precon->oasProlongationKernel(mesh->Nelements, precon->o_oasRestrictionMatrix,
				ellipticOasCoarse->o_x, o_z);
#else
  // copy from coarse device to mesh device
  o_tmpCoarse.copyFrom(ellipticOasCoarse->o_x.ptr());

  // this happens on the mesh device
  precon->oasProlongationKernel(mesh->Nelements, precon->o_oasRestrictionMatrix,
				o_tmpCoarse, o_z);
#endif

  mesh1->device.finish();
  mesh->device.finish();
  meshCoarse->device.finish();

  o_tmpCoarse.free();
  
  if(mesh->rank==0) printf("Ending coarse grid iterations:\n Outer ");  
  
}
