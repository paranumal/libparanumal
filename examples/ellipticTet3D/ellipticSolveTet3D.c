#include "ellipticTet3D.h"

void ellipticStartHaloExchange3D(mesh3D *mesh, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer);

void ellipticEndHaloExchange3D(mesh3D *mesh, occa::memory &o_q, dfloat *recvBuffer);

void ellipticParallelGatherScatterTet3D(mesh3D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op);

void ellipticOperator3D(solver_t *solver, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *options){

  mesh_t *mesh = solver->mesh;

  mesh->device.finish();
  occa::tic("AxKernel");

  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;
  
  if(strstr(options, "IPDG")){
    // halo exchange prior to A*q
    ellipticStartHaloExchange3D(mesh, o_q, sendBuffer, recvBuffer);
    
    ellipticEndHaloExchange3D(mesh, o_q, recvBuffer);

    mesh->device.finish();
    occa::tic("gradientKernel");    
    
    iint allNelements = mesh->Nelements+mesh->totalHaloPairs; 
    mesh->gradientKernel(allNelements,
			 mesh->o_vgeo,
			 mesh->o_DrT,
			 mesh->o_DsT,
       mesh->o_DtT,
			 o_q,
			 solver->o_grad);

    mesh->device.finish();
    occa::toc("gradientKernel");
    occa::tic("ipdgKernel");
    
    // TW NOTE WAS 2 !
    dfloat tau = 2.f*(mesh->N+1)*(mesh->N+3)/3.; // 1/h factor built into kernel 
    mesh->ipdgKernel(mesh->Nelements,
		     mesh->o_vmapM,
		     mesh->o_vmapP,
		     lambda,
		     tau,
		     mesh->o_vgeo,
		     mesh->o_sgeo,
		     mesh->o_DrT,
		     mesh->o_DsT,
         mesh->o_DtT,
		     mesh->o_LIFTT,
		     mesh->o_MM,
		     solver->o_grad,
		     o_Aq);

    mesh->device.finish();
    occa::toc("ipdgKernel");
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

  if(strstr(options,"PROJECT")){
      mesh->device.finish();
      occa::tic("Project");

      // r <= S*G*r 
      ellipticParallelGatherScatterTet3D(mesh, solver->ogs, o_r, o_r, dfloatString, "add");
      // r <= W*S*G*r (Project^t*r)
      mesh->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->o_projectL2, o_r, o_r);
      mesh->device.finish();
      occa::toc("Project");
    }

  if(strstr(options, "OAS")){
    
    ellipticStartHaloExchange3D(mesh, o_r, sendBuffer, recvBuffer);
    
    ellipticEndHaloExchange3D(mesh, o_r, recvBuffer);

    mesh->device.finish();
    occa::tic("OASpreconKernel");

    precon->preconKernel(mesh->Nelements,
			 mesh->o_vmapP,
			 precon->o_oasForwardDgT,
			 precon->o_oasDiagInvOpDg,
			 precon->o_oasBackDgT,
			 o_r,
			 o_zP);
    
    // OAS => additive Schwarz => sum up patch solutions
    ellipticParallelGatherScatterTet3D(mesh, precon->ogsDg, o_zP, o_zP, solver->type, "add");

    mesh->device.finish();
    occa::toc("OASpreconKernel");

    // extract block interiors on DEVICE
    mesh->device.finish();
    occa::tic("restrictKernel");

    precon->restrictKernel(mesh->Nelements, o_zP, o_z);

    mesh->device.finish();
    occa::toc("restrictKernel");   
    
    if(strstr(options, "COARSEGRID")){ // should split into two parts

      mesh->device.finish();
      occa::tic("coarseGrid");

      mesh->device.finish();
      occa::tic("coarsenKernel");
      
      // Z1*Z1'*PL1*(Z1*z1) = (Z1*rL)  HMMM
      precon->coarsenKernel(mesh->Nelements, precon->o_coarseInvDegree, precon->o_V1, o_r, precon->o_r1);

      mesh->device.finish();
      occa::toc("coarsenKernel");
      

      // solve coarse problem using xxt
      if(strstr(options, "XXT")){
        precon->o_r1.copyTo(precon->r1); 
        occa::tic("xxtSolve");
        xxtSolve(precon->z1, precon->xxt,precon->r1);
        occa::toc("xxtSolve");
        precon->o_z1.copyFrom(precon->z1);
      }

      if(strstr(options,"GLOBALALMOND")){
        precon->o_r1.copyTo(precon->r1); 
      	occa::tic("parALMOND");
      	almondSolve(precon->z1, precon->parAlmond, precon->r1);
      	occa::toc("parALMOND");
        precon->o_z1.copyFrom(precon->z1);
      }

      if(strstr(options,"LOCALALMOND")){
        precon->o_r1.copyTo(precon->r1); 
        occa::tic("ALMOND");
        almondSolve(precon->z1, precon->almond, precon->r1);
        occa::toc("ALMOND");
        precon->o_z1.copyFrom(precon->z1);
      }
      
      // prolongate from P1 to PN
      occa::tic("prolongateKernel");
      precon->prolongateKernel(mesh->Nelements, precon->o_V1, precon->o_z1, precon->o_ztmp);
      mesh->device.finish();
      occa::toc("prolongateKernel");

      // increment z
      dfloat one = 1.;
      ellipticScaledAdd(solver, one, precon->o_ztmp, one, o_z);
      
      mesh->device.finish();
      occa::toc("coarseGrid");
    }

    if(strstr(options,"PROJECT")){
      mesh->device.finish();
      occa::tic("Project");
      mesh->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->o_projectL2, o_z, o_z);
      ellipticParallelGatherScatterTet3D(mesh, solver->ogs, o_z, o_z, dfloatString, "add");
      mesh->device.finish();
      occa::toc("Project");
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


int ellipticSolveTet3D(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const char *options){

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
  
  mesh->device.finish();
  occa::tic("PCG");

  // compute A*x
  ellipticOperator3D(solver, lambda, o_x, solver->o_Ax, options);
  
  // subtract r = b - A*x
  ellipticScaledAdd(solver, -1.f, o_Ax, 1.f, o_r);

  mesh->device.finish();
  occa::tic("Preconditioner");
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
  mesh->device.finish();
  occa::toc("Preconditioner");

  
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
    ellipticOperator3D(solver, lambda, o_p, o_Ap, options); 
   
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

    mesh->device.finish();
    occa::tic("Preconditioner");
    if(strstr(options,"PCG")){

      // z = Precon^{-1} r
      ellipticPreconditioner3D(solver, o_r, o_zP, o_z, options);

      // dot(r,z)
      rdotz1 = ellipticInnerProduct(solver, o_r, o_z);

      // flexible pcg beta = (z.(-alpha*Ap))/zdotz0
      if(strstr(options,"FLEXIBLE")){
        dfloat zdotAp = ellipticInnerProduct(solver, o_z, o_Ap);
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
    mesh->device.finish();
    occa::toc("Preconditioner");
    
    // switch rdotz0,rdotr0 <= rdotz1,rdotr1
    rdotr0 = rdotr1;
    
    if(rank==0)
      printf("iter=%05d pAp = %g norm(r) = %g\n", Niter, pAp, sqrt(rdotr0));

    ++Niter;
    
  }while(rdotr0>(tol*tol));

  mesh->device.finish();
  occa::toc("PCG");

  occa::printTimer();

  printf("total number of nodes: %d\n", mesh->Np*mesh->Nelements);

  return Niter;
}
