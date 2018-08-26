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

#include "massHex3D.h"


void massOperator3D(solver_t *solver, dfloat lambda,
		    occa::memory &o_q, occa::memory &o_Aq, const char *options){

  mesh_t *mesh = solver->mesh;

  //  occaTimerTic(mesh->device,"AxKernel");

  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;

  // compute local element operations and store result in o_Aq
  if(strstr(options, "CONTINUOUS")){
    ogs_t *nonHalo = solver->nonHalo;
    ogs_t *halo = solver->halo;

    // Ax for C0 halo elements  (on default stream - otherwise local Ax swamps)
    mesh->device.setStream(solver->dataStream);
    mesh->device.finish();
    mesh->device.setStream(solver->defaultStream);
    mesh->device.finish();

    dfloat zero = 0;
    solver->o_pAp.copyFrom(&zero);
    {
      if(solver->NglobalGatherElements){
	solver->partialAxKernel(solver->NglobalGatherElements, solver->o_globalGatherElementList,
				solver->o_gjGeo, solver->o_gjI, o_q, o_Aq, solver->o_grad, solver->o_Aw);
      }
      
      if(halo->Ngather){
	mesh->gatherKernel(halo->Ngather, halo->o_gatherOffsets, halo->o_gatherLocalIds, o_Aq, halo->o_gatherTmp);
      }
      
      if(halo->Ngather){
	halo->o_gatherTmp.copyTo(halo->gatherTmp);
      }
      
      // Ax for C0 internal elements
      if(solver->NlocalGatherElements){
	solver->partialAxKernel(solver->NlocalGatherElements, solver->o_localGatherElementList,
				solver->o_gjGeo, solver->o_gjI, o_q, o_Aq, solver->o_grad, solver->o_Aw);
      }
    }

    // C0 halo gather-scatter (on data stream)
    if(halo->Ngather){
      occa::streamTag tag;
      
      // MPI based gather scatter using libgs
      gsParallelGatherScatter(halo->gatherGsh, halo->gatherTmp, dfloatString, "add");
      
      // copy totally gather halo data back from HOST to DEVICE
      halo->o_gatherTmp.copyFrom(halo->gatherTmp);
      
      // wait for async copy
      tag = mesh->device.tagStream();
      mesh->device.waitFor(tag);
      
      // do scatter back to local nodes
      mesh->scatterKernel(halo->Ngather, halo->o_gatherOffsets, halo->o_gatherLocalIds, halo->o_gatherTmp, o_Aq);
      
      // make sure the scatter has finished on the data stream
      tag = mesh->device.tagStream();
      mesh->device.waitFor(tag);
    }
    
    // finalize gather using local and global contributions
    mesh->device.setStream(solver->defaultStream);
#if 0
    if(nonHalo->Ngather)
      mesh->gatherScatterKernel(nonHalo->Ngather, nonHalo->o_gatherOffsets, nonHalo->o_gatherLocalIds, o_Aq);
#endif
  }
  else{


  }

  //  occaTimerToc(mesh->device,"AxKernel");
}


dfloat massScaledAdd(solver_t *solver, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b){

  mesh_t *mesh = solver->mesh;

  int Ntotal = mesh->Nelements*mesh->Np;

  occaTimerTic(mesh->device,"scaledAddKernel");

  // b[n] = alpha*a[n] + beta*b[n] n\in [0,Ntotal)
  mesh->scaledAddKernel(Ntotal, alpha, o_a, beta, o_b);

  occaTimerToc(mesh->device,"scaledAddKernel");

}

dfloat massWeightedInnerProduct(solver_t *solver,
				occa::memory &o_w,
				occa::memory &o_a,
				occa::memory &o_b,
				const char *options){


  mesh_t *mesh = solver->mesh;
  dfloat *tmp = solver->tmp;
  int Nblock = solver->Nblock;
  int Ntotal = mesh->Nelements*mesh->Np;

  occa::memory &o_tmp = solver->o_tmp;

  occaTimerTic(mesh->device,"weighted inner product2");
  //  printf("Nblock = %d, Ntotal = %d, ratio = %lf\n", Nblock, Ntotal, ((double)Ntotal)/Nblock);
  if(strstr(options,"CONTINUOUS")||strstr(options, "PROJECT"))
    mesh->weightedInnerProduct2Kernel(Ntotal, o_w, o_a, o_b, o_tmp);
  else
    mesh->innerProductKernel(Ntotal, o_a, o_b, o_tmp);

  occaTimerToc(mesh->device,"weighted inner product2");

  o_tmp.copyTo(tmp);

  dfloat wab = 0;
  for(int n=0;n<Nblock;++n){
    wab += tmp[n];
  }

  dfloat globalwab = 0;
  MPI_Allreduce(&wab, &globalwab, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);


  return globalwab;
}


void massPreconditioner3D(solver_t *solver,
			  occa::memory &o_r,
			  occa::memory &o_zP,
			  occa::memory &o_z,
			  const char *options){

  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  ogs_t    *ogs = solver->ogs; // C0 Gather ScatterTri info

  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;

  if(strstr(options, "JACOBI")){

    int Ntotal = mesh->Np*mesh->Nelements;
    // Jacobi preconditioner
    occaTimerTic(mesh->device,"dotDivideKernel");
    mesh->dotDivideKernel(Ntotal, o_r, precon->o_diagA, o_z);
    occaTimerToc(mesh->device,"dotDivideKernel");
  }
  else // turn off preconditioner
    o_z.copyFrom(o_r);

}

int massSolveHex3D(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const int maxIterations, const char *options){

  mesh_t *mesh = solver->mesh;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // convergence tolerance (currently absolute)
  const dfloat tol = 1e-10;

  occa::memory &o_p  = solver->o_p;
  occa::memory &o_z  = solver->o_z;
  occa::memory &o_zP = solver->o_zP;
  occa::memory &o_Ap = solver->o_Ap;
  occa::memory &o_Ax = solver->o_Ax;

  occa::streamTag startTag = mesh->device.tagStream();

  occaTimerTic(mesh->device,"PCG");

  mesh->device.setStream(solver->defaultStream);

  // gather-scatter
  if(strstr(options,"CONTINUOUS")||strstr(options, "PROJECT"))
    massParallelGatherScatter(mesh, solver->ogs, o_r, o_r, dfloatString, "add");

  // compute A*x
  massOperator3D(solver, lambda, o_x, o_Ax, options);

  // subtract r = b - A*x
  massScaledAdd(solver, -1.f, o_Ax, 1.f, o_r);

  occaTimerTic(mesh->device,"Preconditioner");
  if(strstr(options,"PCG")){
    // Precon^{-1} (b-A*x)
    massPreconditioner3D(solver, o_r, o_zP, o_z, options); // r => rP => zP => z

    // p = z
    o_p.copyFrom(o_z); // PCG
  }
  else{
    // p = r
    o_p.copyFrom(o_r); // CG
  }
  occaTimerToc(mesh->device,"Preconditioner");

  // dot(r,r)
  dfloat rdotr0 = massWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_r, options);
  dfloat rdotz0 = massWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_z, options);
  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;
  int Niter = 0;
  dfloat alpha, beta, pAp;

  while(Niter<maxIterations && rdotr0>(tol*tol)){
    // A*p
    massOperator3D(solver, lambda, o_p, o_Ap, options);

    // dot(p,A*p)
    pAp = massWeightedInnerProduct(solver, solver->o_invDegree, o_p, o_Ap, options);

    if(strstr(options,"PCG"))
      // alpha = dot(r,z)/dot(p,A*p)
      alpha = rdotz0/pAp;
    else
      // alpha = dot(r,r)/dot(p,A*p)
      alpha = rdotr0/pAp;

    // x <= x + alpha*p
    massScaledAdd(solver,  alpha, o_p,  1.f, o_x);

    // r <= r - alpha*A*p
    massScaledAdd(solver, -alpha, o_Ap, 1.f, o_r);

    // dot(r,r)
    rdotr1 = massWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_r, options);

    occaTimerTic(mesh->device,"Preconditioner");
    if(strstr(options,"PCG")){

      // z = Precon^{-1} r
      massPreconditioner3D(solver, o_r, o_zP, o_z, options);

      // dot(r,z)
      rdotz1 = massWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_z, options);

      // flexible pcg beta = (z.(-alpha*Ap))/zdotz0
      if(strstr(options,"FLEXIBLE")){
	dfloat zdotAp = massWeightedInnerProduct(solver, solver->o_invDegree, o_z, o_Ap, options);
	beta = -alpha*zdotAp/rdotz0;
      }
      else{
	beta = rdotz1/rdotz0;
      }

      // p = z + beta*p
      massScaledAdd(solver, 1.f, o_z, beta, o_p);

      // switch rdotz0 <= rdotz1
      rdotz0 = rdotz1;
    }
    else{
      beta = rdotr1/rdotr0;

      // p = r + beta*p
      massScaledAdd(solver, 1.f, o_r, beta, o_p);
    }
    occaTimerToc(mesh->device,"Preconditioner");

    // switch rdotr0 <= rdotr1
    rdotr0 = rdotr1;

#if 0
    if(rank==0)
      printf("iter=%05d pAp = %g norm(r) = %g\n", Niter, pAp, sqrt(rdotr0));
#endif
    ++Niter;
  };

  occaTimerToc(mesh->device,"PCG");

  occa::streamTag stopTag = mesh->device.tagStream();

  double elapsed = mesh->device.timeBetween(startTag, stopTag);
  double gElapsed;
  MPI_Allreduce(&elapsed, &gElapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  occa::printTimer();

  return Niter;
}
