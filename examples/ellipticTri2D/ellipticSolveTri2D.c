#include "ellipticTri2D.h"

void ellipticStartHaloExchange2D(mesh2D *mesh, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer);

void ellipticEndHaloExchange2D(mesh2D *mesh, occa::memory &o_q, dfloat *recvBuffer);

void ellipticParallelGatherScatterTri2D(mesh2D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op);

void ellipticOperator2D(solver_t *solver, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *options){

  mesh_t *mesh = solver->mesh;

  occaTimerTic(mesh->device,"AxKernel");

  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;

  if(strstr(options, "CONTINUOUS")){
    // compute local element operations and store result in o_Aq
    solver->AxKernel(mesh->Nelements,
                   mesh->o_ggeo,
                   mesh->o_SrrT,
                   mesh->o_SrsT,
                   mesh->o_SsrT,
                   mesh->o_SssT,
                   mesh->o_MM,
                   lambda,
                   o_q,
                   o_Aq);

  } else{

    ellipticStartHaloExchange2D(mesh, o_q, sendBuffer, recvBuffer);

    ellipticEndHaloExchange2D(mesh, o_q, recvBuffer);

    occaTimerTic(mesh->device,"gradientKernel");

    iint allNelements = mesh->Nelements+mesh->totalHaloPairs;
    solver->gradientKernel(allNelements,
       mesh->o_vgeo,
       mesh->o_DrT,
       mesh->o_DsT,
       o_q,
       solver->o_grad);

    occaTimerToc(mesh->device,"gradientKernel");
    occaTimerTic(mesh->device,"ipdgKernel");

    solver->ipdgKernel(mesh->Nelements,
		       mesh->o_vmapM,
		       mesh->o_vmapP,
		       lambda,
		       solver->tau,
		       mesh->o_vgeo,
		       mesh->o_sgeo,
		       mesh->o_EToB,
		       mesh->o_DrT,
		       mesh->o_DsT,
		       mesh->o_LIFTT,
		       mesh->o_MM,
		       solver->o_grad,
		       o_Aq);

    occaTimerToc(mesh->device,"ipdgKernel");
  }

  if(strstr(options, "CONTINUOUS")||strstr(options, "PROJECT"))
    // parallel gather scatter
    ellipticParallelGatherScatterTri2D(mesh, solver->ogs, o_Aq, o_Aq, dfloatString, "add");

  occaTimerToc(mesh->device,"AxKernel");
}

void ellipticMatrixFreeAx(void **args, occa::memory o_q, occa::memory o_Aq, const char* options) {

  solver_t *solver = (solver_t *) args[0];
  dfloat  *lambda  = (dfloat *)  args[1];

  mesh_t *mesh = solver->mesh;
  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;

  if(strstr(options, "CONTINUOUS")){
    // compute local element operations and store result in o_Aq
    solver->AxKernel(mesh->Nelements,
                   mesh->o_ggeo,
                   mesh->o_SrrT,
                   mesh->o_SrsT,
                   mesh->o_SsrT,
                   mesh->o_SssT,
                   mesh->o_MM,
                   lambda,
                   o_q,
                   o_Aq);

  } else{

    ellipticStartHaloExchange2D(mesh, o_q, sendBuffer, recvBuffer);

    ellipticEndHaloExchange2D(mesh, o_q, recvBuffer);

    occaTimerTic(mesh->device,"gradientKernel");

    iint allNelements = mesh->Nelements+mesh->totalHaloPairs;
    solver->gradientKernel(allNelements,
       mesh->o_vgeo,
       mesh->o_DrT,
       mesh->o_DsT,
       o_q,
       solver->o_grad);

    occaTimerToc(mesh->device,"gradientKernel");
    occaTimerTic(mesh->device,"ipdgKernel");

    // TW NOTE WAS 2 !
    solver->ipdgKernel(mesh->Nelements,
		       mesh->o_vmapM,
		       mesh->o_vmapP,
		       lambda,
		       solver->tau,
		       mesh->o_vgeo,
		       mesh->o_sgeo,
		       mesh->o_EToB,
		       mesh->o_DrT,
		       mesh->o_DsT,
		       mesh->o_LIFTT,
		       mesh->o_MM,
		       solver->o_grad,
		       o_Aq);

    occaTimerToc(mesh->device,"ipdgKernel");
  }
}

dfloat ellipticScaledAdd(solver_t *solver, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b){

  mesh_t *mesh = solver->mesh;

  iint Ntotal = mesh->Nelements*mesh->Np;

  // b[n] = alpha*a[n] + beta*b[n] n\in [0,Ntotal)
  occaTimerTic(mesh->device,"scaledAddKernel");
  solver->scaledAddKernel(Ntotal, alpha, o_a, beta, o_b);
  occaTimerToc(mesh->device,"scaledAddKernel");
}

dfloat ellipticWeightedInnerProduct(solver_t *solver,
            occa::memory &o_w,
            occa::memory &o_a,
            occa::memory &o_b,
            const char *options){

  mesh_t *mesh = solver->mesh;
  dfloat *tmp = solver->tmp;
  iint Nblock = solver->Nblock;
  iint Ntotal = mesh->Nelements*mesh->Np;

  occa::memory &o_tmp = solver->o_tmp;

  occaTimerTic(mesh->device,"weighted inner product2");

  if(strstr(options,"CONTINUOUS")||strstr(options, "PROJECT"))
    solver->weightedInnerProduct2Kernel(Ntotal, o_w, o_a, o_b, o_tmp);
  else
    solver->innerProductKernel(Ntotal, o_a, o_b, o_tmp);

  occaTimerToc(mesh->device,"weighted inner product2");

  o_tmp.copyTo(tmp);

  dfloat wab = 0;
  for(iint n=0;n<Nblock;++n){
    wab += tmp[n];
  }

  dfloat globalwab = 0;
  MPI_Allreduce(&wab, &globalwab, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

  return globalwab;
}

dfloat ellipticInnerProduct(solver_t *solver,
			    occa::memory &o_a,
			    occa::memory &o_b){


  mesh_t *mesh = solver->mesh;
  dfloat *tmp = solver->tmp;
  iint Nblock = solver->Nblock;
  iint Ntotal = mesh->Nelements*mesh->Np;

  occa::memory &o_tmp = solver->o_tmp;

  occaTimerTic(mesh->device,"inner product");
  solver->innerProductKernel(Ntotal, o_a, o_b, o_tmp);
  occaTimerToc(mesh->device,"inner product");

  o_tmp.copyTo(tmp);

  dfloat ab = 0;
  for(iint n=0;n<Nblock;++n){
    ab += tmp[n];
  }

  dfloat globalab = 0;
  MPI_Allreduce(&ab, &globalab, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

  return globalab;
}

int ellipticSolveTri2D(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const char *options){

  mesh_t *mesh = solver->mesh;

  iint rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // convergence tolerance (currently absolute)
  const dfloat tol = 1e-6;

  // placeholder conjugate gradient:
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method

  // placeholder preconditioned conjugate gradient
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method

  occa::memory &o_p  = solver->o_p;
  occa::memory &o_z  = solver->o_z;
  occa::memory &o_Ap = solver->o_Ap;
  occa::memory &o_Ax = solver->o_Ax;

  occaTimerTic(mesh->device,"PCG");

  // gather-scatter
  if(strstr(options, "CONTINUOUS")||strstr(options, "PROJECT"))
    ellipticParallelGatherScatterTri2D(mesh, solver->ogs, o_r, o_r, dfloatString, "add");

  // compute A*x
  ellipticOperator2D(solver, lambda, o_x, solver->o_Ax, options);

  // subtract r = b - A*x
  ellipticScaledAdd(solver, -1.f, o_Ax, 1.f, o_r);

  dfloat rdotr0 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_r, options);
  dfloat rdotz0 = 0;
  iint Niter = 0;
  //sanity check
  if (rdotr0<=(tol*tol)) {
    printf("iter=0 norm(r) = %g\n", sqrt(rdotr0));
    return 0; 
  }


  occaTimerTic(mesh->device,"Preconditioner");
  if(strstr(options,"PCG")){
    // Precon^{-1} (b-A*x)
    ellipticPreconditioner2D(solver, lambda, o_r, o_z, options); // r => rP => zP => z

    // p = z
    o_p.copyFrom(o_z); // PCG
  }
  else{
    // p = r
    o_p.copyFrom(o_r); // CG
  }
  occaTimerToc(mesh->device,"Preconditioner");


  // dot(r,r)
  rdotz0 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_z, options);
  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;
  
  dfloat alpha, beta, pAp = 0;

  if(rank==0)
    printf("rdotr0 = %g, rdotz0 = %g\n", rdotr0, rdotz0);

  while(rdotr0>(tol*tol)){

    // A*p
    ellipticOperator2D(solver, lambda, o_p, o_Ap, options);

    // dot(p,A*p)
    pAp =  ellipticWeightedInnerProduct(solver, solver->o_invDegree,o_p, o_Ap, options);

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

    if(rdotr1 < tol*tol) {
      rdotr0 = rdotr1;
      Niter++;
      break;
    }

    occaTimerTic(mesh->device,"Preconditioner");
    if(strstr(options,"PCG")){

      // z = Precon^{-1} r
      ellipticPreconditioner2D(solver, lambda, o_r, o_z, options);

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

    // switch rdotz0,rdotr0 <= rdotz1,rdotr1
    rdotr0 = rdotr1;

    if((rank==0)&&(strstr(options,"VERBOSE")))
      printf("iter=%05d pAp = %g norm(r) = %g\n", Niter, pAp, sqrt(rdotr0));

    ++Niter;

  }

  if(rank==0)
    printf("iter=%05d pAp = %g norm(r) = %g\n", Niter, pAp, sqrt(rdotr0));

  occaTimerToc(mesh->device,"PCG");

  occa::printTimer();

  return Niter;

}
