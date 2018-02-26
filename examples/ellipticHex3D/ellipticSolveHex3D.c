#include "ellipticHex3D.h"


void ellipticStartHaloExchange3D(mesh3D *mesh, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer);

void ellipticEndHaloExchange3D(mesh3D *mesh, occa::memory &o_q, dfloat *recvBuffer);

void ellipticParallelGatherScatterHex3D(mesh3D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op);


void ellipticOperator3D(solver_t *solver, dfloat lambda,
      occa::memory &o_q, occa::memory &o_Aq, const char *options){

  mesh_t *mesh = solver->mesh;

  occaTimerTic(mesh->device,"AxKernel");
  
  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;

  // compute local element operations and store result in o_Aq
  if(strstr(options, "CONTINUOUS")){
    mesh->AxKernel(mesh->Nelements, mesh->o_ggeo, mesh->o_D, lambda, o_q, o_Aq);
  }
  else{

    ellipticStartHaloExchange3D(mesh, o_q, sendBuffer, recvBuffer);
    
    ellipticEndHaloExchange3D(mesh, o_q, recvBuffer);

    int allNelements = mesh->Nelements+mesh->totalHaloPairs;
    mesh->gradientKernel(allNelements, mesh->o_vgeo, mesh->o_D, o_q, solver->o_grad);

    dfloat tau = 2.f*(mesh->Nq)*(mesh->Nq+2)/3.;
    mesh->ipdgKernel(mesh->Nelements,
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

  if(strstr(options, "CONTINUOUS")||strstr(options, "PROJECT"))
    // parallel gather scatter
    ellipticParallelGatherScatterHex3D(mesh, solver->ogs, o_Aq, o_Aq, dfloatString, "add");
 
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
    mesh->AxKernel(mesh->Nelements, mesh->o_ggeo, mesh->o_D, lambda, o_q, o_Aq);
  }
  else{

    ellipticStartHaloExchange3D(mesh, o_q, sendBuffer, recvBuffer);
    
    ellipticEndHaloExchange3D(mesh, o_q, recvBuffer);

    int allNelements = mesh->Nelements+mesh->totalHaloPairs;
    mesh->gradientKernel(allNelements, mesh->o_vgeo, mesh->o_D, o_q, solver->o_grad);

    dfloat tau = 2.f*(mesh->Nq)*(mesh->Nq+2)/3.;
    mesh->ipdgKernel(mesh->Nelements,
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

  if(strstr(options, "OAS")){
    //L2 project weighting
    if(strstr(options,"CONTINUOUS")||strstr(options,"PROJECT")) {
      ellipticParallelGatherScatterHex3D(mesh,ogs,o_r,o_r,dfloatString,"add");
      mesh->dotMultiplyKernel(mesh->Np*mesh->Nelements,mesh->o_projectL2,o_r,o_r);
    }

    ellipticStartHaloExchange3D(mesh, o_r, sendBuffer, recvBuffer);
    
    ellipticEndHaloExchange3D(mesh, o_r, recvBuffer);
    
    occaTimerTic(mesh->device,"OASpreconKernel");
    if(strstr(options, "CONTINUOUS")) {
      // compute local precon on DEVICE
      precon->preconKernel(mesh->Nelements,
         precon->o_vmapPP,
         precon->o_faceNodesP,
         precon->o_oasForward,
         precon->o_oasDiagInvOp,
         precon->o_oasBack,
         o_r,
         o_zP);
      ellipticParallelGatherScatterHex3D(mesh, precon->ogsP, o_zP, o_zP, dfloatString, "add");
      mesh->dotMultiplyKernel(mesh->NqP*mesh->NqP*mesh->NqP*mesh->Nelements,precon->o_invDegreeP,o_zP,o_zP);
    } else{
      precon->preconKernel(mesh->Nelements,
         mesh->o_vmapP,
         precon->o_faceNodesP,
         precon->o_oasForwardDg,
         precon->o_oasDiagInvOpDg,
         precon->o_oasBackDg,
         o_r,
         o_zP);
      ellipticParallelGatherScatterHex3D(mesh, precon->ogsDg, o_zP, o_zP, dfloatString, "add");
      mesh->dotMultiplyKernel(mesh->NqP*mesh->NqP*mesh->NqP*mesh->Nelements,precon->o_invDegreeDGP,o_zP,o_zP);
    }
    occaTimerToc(mesh->device,"OASpreconKernel");
    
    // extract block interiors on DEVICE
    occaTimerTic(mesh->device,"restrictKernel");
    precon->restrictKernel(mesh->Nelements, o_zP, o_z);
    occaTimerToc(mesh->device,"restrictKernel");   


    if(strstr(options, "COARSEGRID")){ // should split into two parts
      occaTimerTic(mesh->device,"coarseGrid");

      // Z1*Z1'*PL1*(Z1*z1) = (Z1*rL)  HMMM
      occaTimerTic(mesh->device,"coarsenKernel");
      precon->coarsenKernel(mesh->Nelements, precon->o_V1, o_r, precon->o_r1);
      occaTimerToc(mesh->device,"coarsenKernel");

      if(strstr(options,"XXT")){
        precon->o_r1.copyTo(precon->r1); 
        xxtSolve(precon->z1, precon->xxt,precon->r1);
        precon->o_z1.copyFrom(precon->z1);
      }

      if(strstr(options,"ALMOND")){
        occaTimerTic(mesh->device,"Almond");
        parAlmondPrecon(precon->o_z1, precon->parAlmond, precon->o_r1);
        occaTimerToc(mesh->device,"Almond");
      }

      occaTimerTic(mesh->device,"prolongateKernel");
      precon->prolongateKernel(mesh->Nelements, precon->o_V1, precon->o_z1, precon->o_ztmp);
      occaTimerToc(mesh->device,"prolongateKernel");

      dfloat one = 1.;
      ellipticScaledAdd(solver, one, precon->o_ztmp, one, o_z);
      occaTimerToc(mesh->device,"coarseGrid");
    }

    //project weighting
    if(strstr(options,"CONTINUOUS")||strstr(options,"PROJECT")) {
      mesh->dotMultiplyKernel(mesh->Np*mesh->Nelements,mesh->o_projectL2,o_z,o_z);
      ellipticParallelGatherScatterHex3D(mesh, ogs, o_z, o_z, dfloatString, "add");
    }
  } else if (strstr(options, "FULLALMOND")) {

    occaTimerTic(mesh->device,"parALMOND");
    parAlmondPrecon(o_z, precon->parAlmond, o_r);
    occaTimerToc(mesh->device,"parALMOND");
  
  } else if(strstr(options, "JACOBI")){

    int Ntotal = mesh->Np*mesh->Nelements;
    // Jacobi preconditioner
    occaTimerTic(mesh->device,"dotDivideKernel");   
    mesh->dotDivideKernel(Ntotal, o_r, precon->o_diagA, o_z);
    occaTimerToc(mesh->device,"dotDivideKernel");   
  }
  else // turn off preconditioner
    o_z.copyFrom(o_r);
  
}

int ellipticSolveHex3D(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const char *options){

  mesh_t *mesh = solver->mesh;
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // convergence tolerance (currently absolute)
  const dfloat tol = 1e-6;

  occa::memory &o_p  = solver->o_p;
  occa::memory &o_z  = solver->o_z;
  occa::memory &o_zP = solver->o_zP;
  occa::memory &o_Ap = solver->o_Ap;
  occa::memory &o_Ax = solver->o_Ax;

  occa::streamTag startTag = mesh->device.tagStream();
  
  occaTimerTic(mesh->device,"PCG");

  // gather-scatter
  if(strstr(options,"CONTINUOUS")||strstr(options, "PROJECT"))
    ellipticParallelGatherScatterHex3D(mesh, solver->ogs, o_r, o_r, dfloatString, "add");

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
    
  if(rank==0)
    printf("rdotr0 = %g, rdotz0 = %g\n", rdotr0, rdotz0);

  do{
    // A*p
    ellipticOperator3D(solver, lambda, o_p, o_Ap, options); 

    // dot(p,A*p)
    pAp = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_p, o_Ap, options);

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
    
    if(rdotr1 < tol*tol) break;

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

    if(rank==0)
      printf("iter=%05d pAp = %g norm(r) = %g\n", Niter, pAp, sqrt(rdotr0));
    
    ++Niter;
    
  }while(rdotr0>(tol*tol));

  occaTimerToc(mesh->device,"PCG");

  occa::streamTag stopTag = mesh->device.tagStream();

  double elapsed = mesh->device.timeBetween(startTag, stopTag);
  double gElapsed;
  MPI_Allreduce(&elapsed, &gElapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
  if(rank==0)
    printf("elapsed = %g iter=%05d pAp = %g norm(r) = %g\n",
	   gElapsed, Niter, pAp, sqrt(rdotr0));

  occa::printTimer();

  printf("total number of nodes: %d\n", mesh->Np*mesh->Nelements);

  return Niter;
}
