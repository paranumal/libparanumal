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


void ellipticOperator3D(solver_t *solver, dfloat lambda,
                        occa::memory &o_q, occa::memory &o_Aq, const char *options){
                        
	mesh_t *mesh = solver->mesh;
	
	//  occaTimerTic(mesh->device,"AxKernel");
	
	dfloat *sendBuffer = solver->sendBuffer;
	dfloat *recvBuffer = solver->recvBuffer;
	
	// compute local element operations and store result in o_Aq
	if(strstr(options, "CONTINUOUS")){
		//    mesh->AxKernel(mesh->Nelements, mesh->o_ggeo, mesh->o_D, lambda, o_q, o_Aq);
#if 0
		solver->AxKernel(mesh->Nelements, solver->o_gjGeo, solver->o_gjD2, solver->o_gjI, lambda, o_q, o_Aq);
		
		ellipticParallelGatherScatter(mesh, solver->ogs, o_Aq, o_Aq, dfloatString, "add");
#endif
		
		//    solver->AxKernel(mesh->Nelements, solver->o_gjGeo, solver->o_gjD2, solver->o_gjI, lambda, o_q, o_Aq);
		ogs_t *nonHalo = solver->nonHalo;
		ogs_t *halo = solver->halo;
		
		// Ax for C0 halo elements  (on default stream - otherwise local Ax swamps)
#if 1
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
				if(strstr(options, "COLLOCATION")){
				
					solver->partialAxKernel(solver->NglobalGatherElements, solver->o_globalGatherElementList,
					                        solver->o_gjGeo, solver->o_gjD, lambda, o_q, o_Aq, solver->o_grad);
					                        
				}
				else{
					//e9-e12, changed for e6 to e8
					//  printf("no collocation\n");
					solver->partialAxKernel(solver->NglobalGatherElements, solver->o_globalGatherElementList,
					                        solver->o_gjGeo, solver->o_gjD2, solver->o_gjI, lambda, o_q, o_Aq,
					                        solver->o_grad);
				}
			}
			
			if(halo->Ngather){
				//	mesh->device.setStream(solver->dataStream);
				
				mesh->gatherKernel(halo->Ngather, halo->o_gatherOffsets, halo->o_gatherLocalIds, o_Aq, halo->o_gatherTmp);
			}
			
			if(halo->Ngather){
				//	mesh->device.setStream(solver->dataStream);
				// avoid async copy [ otherwise we compete with the local Ax ]
				halo->o_gatherTmp.copyTo(halo->gatherTmp);
			}
			
			// Ax for C0 internal elements
			if(solver->NlocalGatherElements){
				//	mesh->device.setStream(solver->defaultStream);
				if(strstr(options, "COLLOCATION")){
					//	printf("collocation 2\n");
					solver->partialAxKernel(solver->NlocalGatherElements, solver->o_localGatherElementList,
					                        solver->o_gjGeo, solver->o_gjD, lambda, o_q, o_Aq, solver->o_grad);
					                        
				}
				else{
					//		printf("no collocation\n");
					//e9-e12, changed for e6 to e8
					solver->partialAxKernel(solver->NlocalGatherElements, solver->o_localGatherElementList,
					                        solver->o_gjGeo, solver->o_gjD2, solver->o_gjI, lambda, o_q,
					                        solver->o_grad, o_Aq);
				}
			}
		}
		
		// C0 halo gather-scatter (on data stream)
		if(halo->Ngather){
			occa::streamTag tag;
			
			//      mesh->device.setStream(solver->dataStream);
			//      tag = mesh->device.tagStream();
			//      mesh->device.waitFor(tag);
			
			// MPI based gather scatter using libgs
			gsParallelGatherScatter(halo->gatherGsh, halo->gatherTmp, dfloatString, "add");
			
			// copy totally gather halo data back from HOST to DEVICE
			//      mesh->device.setStream(solver->dataStream);
			halo->o_gatherTmp.copyFrom(halo->gatherTmp);
			
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
#if 0
		if(nonHalo->Ngather)
			mesh->gatherScatterKernel(nonHalo->Ngather, nonHalo->o_gatherOffsets, nonHalo->o_gatherLocalIds, o_Aq);
#endif
	}
	else{
		// should not be hard coded
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
	
	//  occaTimerToc(mesh->device,"AxKernel");
}


void ellipticMatrixFreeAx(void **args, occa::memory o_q, occa::memory o_Aq, const char* options) {

	solver_t *solver = (solver_t *) args[0];
	dfloat  *lambda  = (dfloat *)  args[1];
	
	mesh_t *mesh = solver->mesh;
	dfloat *sendBuffer = solver->sendBuffer;
	dfloat *recvBuffer = solver->recvBuffer;
	
	// compute local element operations and store result in o_Aq
	if(strstr(options, "CONTINUOUS")){
		// BROKEN
		//    mesh->AxKernel(mesh->Nelements, mesh->o_ggeo, mesh->o_D, lambda, o_q, o_Aq);
		solver->AxKernel(mesh->Nelements, solver->o_gjGeo, solver->o_gjD2, solver->o_gjI, lambda, o_q, o_Aq,
		                 solver->o_invDegree, solver->o_grad);
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
	
	// convergence tolerance (currently absolute)
	const dfloat tol = 1e-8;
	
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
		ellipticParallelGatherScatter(mesh, solver->ogs, o_r, o_r, dfloatString, "add");
		
	// compute A*x
	ellipticOperator3D(solver, lambda, o_x, o_Ax, options);
	
	// subtract r = b - A*x
	ellipticScaledAdd(solver, -1.f, o_Ax, 1.f, o_r);
	
	occaTimerTic(mesh->device,"Preconditioner");
	if(strstr(options,"PCG")){
		// Precon^{-1} (b-A*x)
		ellipticPreconditioner3D(solver, o_r, o_zP, o_z, options); // r => rP => zP => z
		
		// p = z
		o_p.copyFrom(o_z); // PCG
	}
	else{
		// p = r
		o_p.copyFrom(o_r); // CG
	}
	occaTimerToc(mesh->device,"Preconditioner");
	
	// dot(r,r)
	dfloat rdotr0 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_r, options);
	dfloat rdotz0 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_z, options);
	dfloat rdotr1 = 0;
	dfloat rdotz1 = 0;
	int Niter = 0;
	dfloat alpha, beta, pAp;
	
	//  if(rank==0)
	//    printf("rdotr0 = %g, rdotz0 = %g\n", rdotr0, rdotz0);
	
	while(Niter<maxIterations){ // rdotr0>(tol*tol) &&
		// A*p
		ellipticOperator3D(solver, lambda, o_p, o_Ap, options);
		
		// dot(p,A*p)
		// these are only equivalent if p is continuous
#if 0
		if(strstr(options,"CONTINUOUS")){
			dfloat localpAp = 0;
			solver->o_pAp.copyTo(&localpAp);
			MPI_Allreduce(&localpAp, &pAp, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
			printf("pAp=%g\n", pAp);
		}else
#endif
		{
			pAp = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_p, o_Ap, options);
		}
		
		if(strstr(options,"PCG"))
			// alpha = dot(r,z)/dot(p,A*p)
			alpha = rdotz0/pAp;
		else
			// alpha = dot(r,r)/dot(p,A*p)
			alpha = rdotr0/pAp;
			
		// x <= x + alpha*p
		ellipticScaledAdd(solver,  alpha, o_p,  1.f, o_x);
		
		// r <= r - alpha*A*p
		ellipticScaledAdd(solver, -alpha, o_Ap, 1.f, o_r);
		
		// dot(r,r)
		rdotr1 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_r, options);
		
		//    if(rdotr1 < tol*tol) break;
		
		occaTimerTic(mesh->device,"Preconditioner");
		if(strstr(options,"PCG")){
		
			// z = Precon^{-1} r
			ellipticPreconditioner3D(solver, o_r, o_zP, o_z, options);
			
			// dot(r,z)
			rdotz1 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_z, options);
			
			// flexible pcg beta = (z.(-alpha*Ap))/zdotz0
			if(strstr(options,"FLEXIBLE")){
				dfloat zdotAp = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_z, o_Ap, options);
				beta = -alpha*zdotAp/rdotz0;
			}
			else{
				beta = rdotz1/rdotz0;
			}
			
			// p = z + beta*p
			ellipticScaledAdd(solver, 1.f, o_z, beta, o_p);
			
			// switch rdotz0 <= rdotz1
			rdotz0 = rdotz1;
		}
		else{
			beta = rdotr1/rdotr0;
			
			// p = r + beta*p
			ellipticScaledAdd(solver, 1.f, o_r, beta, o_p);
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
	
	//  if(rank==0)
	//    printf("elapsed = %g iter=%05d pAp = %g norm(r) = %g\n",
	//	   gElapsed, Niter, pAp, sqrt(rdotr0));
	
	occa::printTimer();
	
	//  printf("total number of nodes: %d\n", mesh->Np*mesh->Nelements);
	
	return Niter;
}