#include "ellipticTri2D.h"

void ellipticStartHaloExchange2D(mesh2D *mesh, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer);

void ellipticEndHaloExchange2D(mesh2D *mesh, occa::memory &o_q, dfloat *recvBuffer);

void ellipticParallelGatherScatterTri2D(mesh2D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op);

void ellipticOperator2D(solver_t *solver, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *options){

  mesh_t *mesh = solver->mesh;
  
  mesh->device.finish();
  occa::tic("AxKernel");

  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;
  
  if(strstr(options, "IPDG")){
    // halo exchange prior to A*q
    ellipticStartHaloExchange2D(mesh, o_q, sendBuffer, recvBuffer);
    
    ellipticEndHaloExchange2D(mesh, o_q, recvBuffer);
    
    iint allNelements = mesh->Nelements+mesh->totalHaloPairs; 
    mesh->gradientKernel(allNelements,
			 mesh->o_vgeo,
			 mesh->o_DrT,
			 mesh->o_DsT,
			 o_q,
			 solver->o_grad);
    
    // TW NOTE WAS 2 !
    dfloat tau = 2.f*(mesh->N+1)*(mesh->N+1); // 1/h factor built into kernel 
    mesh->ipdgKernel(mesh->Nelements,
		     mesh->o_vmapM,
		     mesh->o_vmapP,
		     lambda,
		     tau,
		     mesh->o_vgeo,
		     mesh->o_sgeo,
		     mesh->o_DrT,
		     mesh->o_DsT,
		     mesh->o_LIFTT,
		     mesh->o_MM,
		     solver->o_grad,
		     o_Aq);
  }
  else{
    
    
  }
  mesh->device.finish();
  occa::toc("AxKernel");
  
}

dfloat ellipticScaledAdd(solver_t *solver, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b){

  mesh_t *mesh = solver->mesh;
  
  iint Ntotal = mesh->Nelements*mesh->Np;

  mesh->device.finish();
  occa::tic("scaledAddKernel");
  
  // b[n] = alpha*a[n] + beta*b[n] n\in [0,Ntotal)
  mesh->scaledAddKernel(Ntotal, alpha, o_a, beta, o_b);

  mesh->device.finish();
  occa::toc("scaledAddKernel");
  
}

dfloat ellipticInnerProduct(solver_t *solver,
			    occa::memory &o_a,
			    occa::memory &o_b){


  mesh_t *mesh = solver->mesh;
  dfloat *tmp = solver->tmp;
  iint Nblock = solver->Nblock;
  iint Ntotal = mesh->Nelements*mesh->Np;

  occa::memory &o_tmp = solver->o_tmp;
  
  mesh->device.finish();
  occa::tic("weighted inner product2");

  mesh->innerProductKernel(Ntotal, o_a, o_b, o_tmp);
  
  mesh->device.finish();
  occa::toc("weighted inner product2");
  
  o_tmp.copyTo(tmp);

  dfloat ab = 0;
  for(iint n=0;n<Nblock;++n){
    ab += tmp[n];
  }
      
  dfloat globalab = 0;
  MPI_Allreduce(&ab, &globalab, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

  return globalab;
}

void ellipticPreconditioner2D(solver_t *solver,
			      occa::memory &o_r,
			      occa::memory &o_zP,
			      occa::memory &o_z,
			      const char *options){

  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  ogs_t    *ogs = solver->ogs; // C0 Gather ScatterTri info
  
  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;
  
  if(strstr(options,"PROJECT")){
    // S*G*r
    ellipticParallelGatherScatterTri2D(mesh, ogs, o_r, o_r, dfloatString, "add");

    // W*S*G*r (Project^t*r)
    mesh->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->o_projectL2, o_r, o_r);
  }
  
  if(strstr(options, "OAS")){
    
    ellipticStartHaloExchange2D(mesh, o_r, sendBuffer, recvBuffer);
    
    ellipticEndHaloExchange2D(mesh, o_r, recvBuffer);

    mesh->device.finish();
    occa::tic("preconKernel");

    precon->preconKernel(mesh->Nelements,
			 mesh->o_vmapP,
			 precon->o_oasForwardDgT,
			 precon->o_oasDiagInvOpDg,
			 precon->o_oasBackDgT,
			 o_r,
			 o_zP);
    
    mesh->device.finish();
    occa::toc("preconKernel");
    
    // OAS => additive Schwarz => sum up patch solutions
    ellipticParallelGatherScatterTri2D(mesh, precon->ogsDg, o_zP, o_zP, solver->type, "add");

    // extract block interiors on DEVICE
    mesh->device.finish();
    occa::tic("restrictKernel");

    precon->restrictKernel(mesh->Nelements, o_zP, o_z);

    mesh->device.finish();
    occa::toc("restrictKernel");   
    
    if(strstr(options, "COARSEGRID")){ // should split into two parts

      mesh->device.finish();
      occa::tic("coarseGrid");
      
      // Z1*Z1'*PL1*(Z1*z1) = (Z1*rL)  HMMM
      precon->coarsenKernel(mesh->Nelements, precon->o_coarseInvDegree, precon->o_V1, o_r, precon->o_r1);

      // do we need to gather (or similar) here ?
      precon->o_r1.copyTo(precon->r1); 
      
      xxtSolve(precon->z1, precon->xxt,precon->r1);

      precon->o_z1.copyFrom(precon->z1);
      
      precon->prolongateKernel(mesh->Nelements, precon->o_V1, precon->o_z1, precon->o_ztmp);

      dfloat one = 1.;
      ellipticScaledAdd(solver, one, precon->o_ztmp, one, o_z);

      occa::toc("coarseGrid");
    }

    if(strstr(options,"PROJECT")){
      mesh->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->o_projectL2, o_z, o_z);
      ellipticParallelGatherScatterTri2D(mesh, ogs, o_z, o_z, dfloatString, "add");
    }
    
  }
  else if(strstr(options, "JACOBI")){

    occa::tic("dotDivideKernel");   
    mesh->device.finish();
    
    // Jacobi preconditioner
    iint Ntotal = mesh->Np*mesh->Nelements;
    mesh->dotDivideKernel(Ntotal, o_r, precon->o_diagA, o_z);

    occa::toc("dotDivideKernel");   
    mesh->device.finish();
  }
  else{ // turn off preconditioner
    o_z.copyFrom(o_r);
  }
}


int ellipticSolveTri2D(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const char *options){

  mesh_t *mesh = solver->mesh;
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // convergence tolerance (currently absolute)
  const dfloat tol = 1e-6;

  // placeholder conjugate gradient:
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method
  
  // placeholder preconditioned conjugate gradient
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method

  occa::memory &o_p  = solver->o_p;
  occa::memory &o_z  = solver->o_z;
  occa::memory &o_zP = solver->o_zP;
  occa::memory &o_Ap = solver->o_Ap;
  occa::memory &o_Ax = solver->o_Ax;
  
  // compute A*x
  ellipticOperator2D(solver, lambda, o_x, solver->o_Ax, options);
  
  // subtract r = b - A*x
  ellipticScaledAdd(solver, -1.f, o_Ax, 1.f, o_r);

  // Precon^{-1} (b-A*x)
  ellipticPreconditioner2D(solver, o_r, o_zP, o_z, options); // r => rP => zP => z
  
  // p = z
  solver->o_p.copyFrom(solver->o_z);

  // dot(r,r)
  dfloat rdotr0 = ellipticInnerProduct(solver, o_r, o_r);
  dfloat rdotz0 = ellipticInnerProduct(solver, o_r, o_z);
  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;
  iint Niter = 0;
  dfloat alpha, beta;
  
  if(rank==0)
    printf("rdotr0 = %g, rdotz0 = %g\n", rdotr0, rdotz0);
  
  do{
    
    // A*p
    ellipticOperator2D(solver, lambda, o_p, o_Ap, options); 
    
    // dot(p,A*p)
    dfloat pAp =  ellipticInnerProduct(solver, o_p, o_Ap);
    
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
    rdotr1 = ellipticInnerProduct(solver, o_r, o_r);
    
    if(rdotr1 < tol*tol) break;

    // z = Precon^{-1} r
    ellipticPreconditioner2D(solver, o_r, o_zP, o_z, options);
    
    // dot(r,z)
    rdotz1 = ellipticInnerProduct(solver, o_r, o_z);
    
    beta = rdotz1/rdotz0;
    
    // p = z + beta*p
    ellipticScaledAdd(solver, 1.f, o_z, beta, o_p);
    
    // switch rdotz0,rdotr0 <= rdotz1,rdotr1
    rdotz0 = rdotz1;
    rdotr0 = rdotr1;
    
    if(rank==0)
      printf("iter=%05d pAp = %g norm(r) = %g\n", Niter, pAp, sqrt(rdotr0));

    ++Niter;
    
  }while(rdotr0>(tol*tol));

  occa::printTimer();

  printf("total number of nodes: %d\n", mesh->Np*mesh->Nelements);

  return Niter;
}
