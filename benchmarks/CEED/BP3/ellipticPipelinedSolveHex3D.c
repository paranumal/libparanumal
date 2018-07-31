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

#include "ellipticHex3D.h"

#define COMMS 1


void ellipticCombinedUpdateStart(solver_t *solver, dfloat alpha, dfloat beta, 
				 dfloat *localResults, occa::memory &o_results, dfloat *globalResults,
				 occa::memory &o_Aw, occa::memory &o_p, occa::memory &o_r, 
				 occa::memory &o_s,  occa::memory &o_w, occa::memory &o_x, occa::memory &o_z,
				 MPI_Request *request, const char *options){


  int Ntotal = solver->mesh->Nelements*solver->mesh->Np;
  int degreeWeighted = (strstr(options,"CONTINUOUS")!=NULL);
  
  localResults[0] = 0;
  localResults[1] = 0;
  o_results.copyFrom(localResults);
  
  solver->combinedUpdateKernel(Ntotal, degreeWeighted, solver->o_invDegree, alpha, beta,
			       o_Aw, o_p, o_r, o_s, o_w, o_x, o_z, o_results);
  
  o_results.copyTo(localResults);
  //MPI_Iallreduce(localResults, globalResults, 2, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD, request);
  MPI_Allreduce(localResults, globalResults, 2, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
  
}

void ellipticCombinedUpdateEnd(solver_t *solver, dfloat *globalResults, dfloat *gamm, dfloat *delta, MPI_Request *request){
  
  MPI_Status status;
  //  MPI_Wait(request, &status);
  
  *gamm = globalResults[0];
  *delta = globalResults[1];
}



void ellipticOperator3D(solver_t *solver, dfloat lambda,
      occa::memory &o_q, occa::memory &o_Aq, const char *options){

  mesh_t *mesh = solver->mesh;

  occaTimerTic(mesh->device,"AxKernel");
  
  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;

  // compute local element operations and store result in o_Aq
  if(strstr(options, "CONTINUOUS")){
    //    mesh->AxKernel(mesh->Nelements, mesh->o_ggeo, mesh->o_D, lambda, o_q, o_Aq);
#if 0
    solver->AxKernel(mesh->Nelements, solver->o_gjGeo, solver->o_gjD, solver->o_gjI, lambda, o_q, o_Aq);

    ellipticParallelGatherScatter(mesh, solver->ogs, o_Aq, o_Aq, dfloatString, "add");
#else

    //    solver->AxKernel(mesh->Nelements, solver->o_gjGeo, solver->o_gjD, solver->o_gjI, lambda, o_q, o_Aq);
    ogs_t *nonHalo = solver->nonHalo;
    ogs_t *halo = solver->halo;

    // Ax for C0 halo elements  (on default stream - otherwise local Ax swamps)
#if 0
    mesh->device.setStream(solver->dataStream);
    mesh->device.finish();
    mesh->device.setStream(solver->defaultStream);
    mesh->device.finish();
#endif
    dfloat zero = 0;
    solver->o_pAp.copyFrom(&zero);
    {
      if(solver->NglobalGatherElements){
	//	mesh->device.setStream(solver->dataStream);
      
	solver->partialAxKernel(solver->NglobalGatherElements, solver->o_globalGatherElementList,
				solver->o_gjGeo, solver->o_gjD, solver->o_gjI, lambda, o_q, o_Aq, 
				solver->o_pAp);
      }

      if(halo->Ngather){
	//	mesh->device.setStream(solver->dataStream);
	
	mesh->gatherKernel(halo->Ngather, halo->o_gatherOffsets, halo->o_gatherLocalIds, o_Aq, halo->o_gatherTmp);
      }

#if COMMS==1
      if(halo->Ngather){
	//	mesh->device.setStream(solver->dataStream);
	// avoid async copy [ otherwise we compete with the local Ax ]
	halo->o_gatherTmp.copyTo(halo->gatherTmp);
      }      
#endif  
      // Ax for C0 internal elements
      if(solver->NlocalGatherElements){
	//	mesh->device.setStream(solver->defaultStream);

	solver->partialAxKernel(solver->NlocalGatherElements, solver->o_localGatherElementList,
				solver->o_gjGeo, solver->o_gjD, solver->o_gjI, lambda, o_q, o_Aq, 
				solver->o_pAp);
      } 
    }

    // C0 halo gather-scatter (on data stream)
    if(halo->Ngather){
      occa::streamTag tag;   

      //      mesh->device.setStream(solver->dataStream);
      //      tag = mesh->device.tagStream();
      //      mesh->device.waitFor(tag);

#if COMMS==1
      // MPI based gather scatter using libgs
      gsParallelGatherScatter(halo->gatherGsh, halo->gatherTmp, dfloatString, "add"); 
      
      // copy totally gather halo data back from HOST to DEVICE
      //      mesh->device.setStream(solver->dataStream);
      halo->o_gatherTmp.copyFrom(halo->gatherTmp); 
#endif 
      // wait for async copy
      //      occa::streamTag tag = mesh->device.tagStream();
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
    if(nonHalo->Ngather)
      mesh->gatherScatterKernel(nonHalo->Ngather, nonHalo->o_gatherOffsets, nonHalo->o_gatherLocalIds, o_Aq);
    
#endif    
  }
  else{
    // should not be hard coded
    dfloat tau = 2.f*(mesh->Nq)*(mesh->Nq+2)/3.;

    int offset = 0;

#if COMMS==1    
    ellipticStartHaloExchange3D(solver, o_q, sendBuffer, recvBuffer);
    ellipticInterimHaloExchange3D(solver, o_q, sendBuffer, recvBuffer);
#endif

    solver->partialGradientKernel(mesh->Nelements, offset, mesh->o_vgeo, mesh->o_D, o_q, solver->o_grad);

    if(mesh->NinternalElements)
      solver->partialIpdgKernel(mesh->NinternalElements,
				mesh->o_internalElementIds,
				mesh->o_vmapM,
				mesh->o_vmapP,
				lambda,
				tau,
				mesh->o_vgeo,
				mesh->o_sgeo,
				mesh->o_D,
				solver->o_grad,
				o_Aq);

#if COMMS==1    
    ellipticEndHaloExchange3D(solver, o_q, recvBuffer);
#endif
    
    if(mesh->totalHaloPairs){
      offset = mesh->Nelements;
      solver->partialGradientKernel(mesh->totalHaloPairs, offset, mesh->o_vgeo, mesh->o_D, o_q, solver->o_grad);
    }

    if(mesh->NnotInternalElements)
      solver->partialIpdgKernel(mesh->NnotInternalElements,
				mesh->o_notInternalElementIds,
				mesh->o_vmapM,
				mesh->o_vmapP,
				lambda,
				tau,
				mesh->o_vgeo,
				mesh->o_sgeo,
				mesh->o_D,
				solver->o_grad,
				o_Aq);
    
  }

  occaTimerToc(mesh->device,"AxKernel");
}


void ellipticMatrixFreeAx(void **args, occa::memory o_q, occa::memory o_Aq, const char* options) {

  solver_t *solver = (solver_t *) args[0];
  dfloat  *lambda  = (dfloat *)  args[1];

  mesh_t *mesh = solver->mesh;
  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;

  // compute local element operations and store result in o_Aq
  if(strstr(options, "CONTINUOUS")){
    solver->AxKernel(mesh->Nelements, solver->o_gjGeo, solver->o_gjD, solver->o_gjI, lambda, o_q, o_Aq,
		     solver->o_invDegree, solver->o_pAp);
  }
  else{
    // tau should not be hard coded
    dfloat tau = 2.f*(mesh->Nq)*(mesh->Nq+2)/3.;
    
    int offset = 0;
    
    ellipticStartHaloExchange3D(solver, o_q, sendBuffer, recvBuffer);

    solver->partialGradientKernel(mesh->Nelements, offset, mesh->o_vgeo, mesh->o_D, o_q, solver->o_grad);

    ellipticInterimHaloExchange3D(solver, o_q, sendBuffer, recvBuffer);

    if(mesh->NinternalElements)
      solver->partialIpdgKernel(mesh->NinternalElements,
				mesh->o_internalElementIds,
				mesh->o_vmapM,
				mesh->o_vmapP,
				lambda,
				tau,
				mesh->o_vgeo,
				mesh->o_sgeo,
				mesh->o_D,
				solver->o_grad,
				o_Aq);
   
    ellipticEndHaloExchange3D(solver, o_q, recvBuffer);

    if(mesh->totalHaloPairs){
      offset = mesh->Nelements;      
      solver->partialGradientKernel(mesh->totalHaloPairs, offset, mesh->o_vgeo, mesh->o_D, o_q, solver->o_grad);
    }

#if 1
    if(mesh->NnotInternalElements)
      solver->partialIpdgKernel(mesh->NnotInternalElements,
				mesh->o_notInternalElementIds,
				mesh->o_vmapM,
				mesh->o_vmapP,
				lambda,
				tau,
				mesh->o_vgeo,
				mesh->o_sgeo,
				mesh->o_D,
				solver->o_grad,
				o_Aq);
#else
    solver->ipdgKernel(mesh->Nelements,
		       mesh->o_vmapM,
		       mesh->o_vmapP,
		       lambda,
		       tau,
		       mesh->o_vgeo,
		       mesh->o_sgeo,
		       mesh->o_D,
		       solver->o_grad,
		       o_Aq);
#endif
  }
}


dfloat ellipticScaledAdd(solver_t *solver, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b){

  mesh_t *mesh = solver->mesh;

  int Ntotal = mesh->Nelements*mesh->Np;

  occaTimerTic(mesh->device,"scaledAddKernel");
  
  // b[n] = alpha*a[n] + beta*b[n] n\in [0,Ntotal)
  mesh->scaledAddKernel(Ntotal, alpha, o_a, beta, o_b);

  occaTimerToc(mesh->device,"scaledAddKernel");
  
}

dfloat ellipticWeightedInnerProduct(solver_t *solver,
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
  if(strstr(options,"CONTINUOUS"))
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


dfloat ellipticStartCombinedInnerProduct(solver_t *solver,
					 occa::memory &o_r,
					 occa::memory &o_w,
					 occa::memory &o_results,
					 const char *options){
  
  mesh_t *mesh = solver->mesh;
  dfloat *tmp = solver->tmp;
  int Ntotal = mesh->Nelements*mesh->Np;
  
  occaTimerTic(mesh->device,"weighted inner product2");
  
  int degreeWeighted = (strstr(options,"CONTINUOUS")!=NULL);
  dfloat results[2];
  results[0] = 0;
  results[1] = 0;
  o_results.copyFrom(results);
  
  solver->combinedInnerProductKernel(Ntotal, solver->o_invDegree, o_r, o_w, o_results, degreeWeighted);
  
  occaTimerToc(mesh->device,"weighted inner product2");
}


void ellipticIntermediateCombinedInnerProduct(solver_t *solver,
					      occa::memory &o_results,
					      dfloat *localResults,
					      dfloat *globalResults,
					      MPI_Request *request){

  
  o_results.copyTo(localResults);
  
#if 0
  MPI_Allreduce(localResults, globalResults, 2, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
#else
  MPI_Iallreduce(localResults, globalResults, 2, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD, request);
#endif
}


void ellipticFinishCombinedInnerProduct(solver_t *solver,
					dfloat *globalResults,
					dfloat *rdotr,
					dfloat *rdotw,
					MPI_Request *request){


  MPI_Status status;
  MPI_Wait(request, &status);
  
  *rdotr = globalResults[0];
  *rdotw = globalResults[1];
}


void ellipticPreconditioner3D(solver_t *solver,
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

int ellipticSolveHex3D(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const int maxIterations, const char *options){

  mesh_t *mesh = solver->mesh;
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0)
    printf("Pipelined Conjugate Gradient: \n");
  
  // convergence tolerance (currently absolute)
  const dfloat tol = 1e-6;
  const dfloat one = 1.;
  
  occa::memory &o_p  = solver->o_p;
  occa::memory &o_w  = solver->o_w;
  occa::memory &o_s  = solver->o_s;
  occa::memory &o_z  = solver->o_z;
  occa::memory &o_zP = solver->o_zP;
  occa::memory &o_Aw = solver->o_Aw;
  occa::memory &o_Ax = solver->o_Ax;

  occa::memory o_results = mesh->device.malloc(2*sizeof(dfloat));
  
  occa::streamTag startTag = mesh->device.tagStream();
  
  occaTimerTic(mesh->device,"PCG");

  mesh->device.setStream(solver->defaultStream);
  
  // gather-scatter
  if(strstr(options,"CONTINUOUS"))
    ellipticParallelGatherScatter(mesh, solver->ogs, o_r, o_r, dfloatString, "add");

  // compute A*x
  ellipticOperator3D(solver, lambda, o_x, o_Ax, options);
  
  // subtract r = b - A*x
  ellipticScaledAdd(solver, -one, o_Ax, one, o_r);

  // w = A*r
  ellipticOperator3D(solver, lambda, o_r, o_w, options);

  dfloat gam = 1, delta = 1, alpha, beta, oldgam = 1;
  int Niter = 0;

  MPI_Request *request = (MPI_Request*) calloc(1, sizeof(MPI_Request));

  occa::memory o_localIpBuffer = mesh->device.mappedAlloc(2*sizeof(dfloat), NULL);
  occa::memory o_globalIpBuffer = mesh->device.mappedAlloc(2*sizeof(dfloat), NULL);
  
  dfloat *localResults = (dfloat*) o_localIpBuffer.getMappedPointer();
  dfloat *globalResults = (dfloat*) o_globalIpBuffer.getMappedPointer();
  
  while(Niter<maxIterations){
    // save last gamma
    oldgam = gam;
    
    // q = A*w 
    ellipticOperator3D(solver, lambda, o_w, o_Aw, options); 

#if COMMS==1

    if(Niter>0){

      ellipticCombinedUpdateEnd(solver, globalResults, &gam, &delta, request);

      beta = gam/oldgam;
      alpha = gam/(delta - beta*gam/alpha);
    }
    else{
      
      gam   = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_r, options);
      delta = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_w, options);
      
      beta = 0;
      alpha = gam/delta;
    }
#else
    beta = .1;
    alpha = .2;
#endif

    ellipticCombinedUpdateStart(solver, alpha, beta, localResults, o_results, globalResults,
				o_Aw, o_p, o_r, o_s, o_w, o_x, o_z, request, options);

    if(rank==0)
      printf("iter = %05d alpha = %g, beta = %g, gam = %g\n", Niter, alpha, beta, gam);
    
    ++Niter;
  };

  occaTimerToc(mesh->device,"PCG");

  occa::streamTag stopTag = mesh->device.tagStream();

  double elapsed = mesh->device.timeBetween(startTag, stopTag);
  double gElapsed;
  MPI_Allreduce(&elapsed, &gElapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
  //  if(rank==0)
  //    printf("elapsed = %g iter=%05d pAp = %g norm(r) = %g\n",
  //	   gElapsed, Niter, pAp, sqrt(rdotr0));

  occa::printTimer();

  //  printf("total number of nodes: %d\n", mesh->Np*mesh->Nelements);

  return Niter;
}
