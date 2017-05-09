#include "ellipticQuad2D.h"


void ellipticStartHaloExchange2D(mesh2D *mesh, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer);

void ellipticEndHaloExchange2D(mesh2D *mesh, occa::memory &o_q, dfloat *recvBuffer);

void ellipticParallelGatherScatterQuad2D(mesh2D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op);

void ellipticOperator2D(solver_t *solver, dfloat lambda,
      occa::memory &o_q, occa::memory &o_Aq, const char *options){

  mesh_t *mesh = solver->mesh;

  mesh->device.finish();
  occa::tic("AxKernel");

  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;

  if(strstr(options, "CONTINUOUS")){
    // compute local element operations and store result in o_Aq
    mesh->AxKernel(mesh->Nelements, mesh->o_ggeo, mesh->o_D, lambda, o_q, o_Aq);
    
  } else{

    ellipticStartHaloExchange2D(mesh, o_q, sendBuffer, recvBuffer);
    
    ellipticEndHaloExchange2D(mesh, o_q, recvBuffer);

    // need start/end elements then can split into two parts
    iint allNelements = mesh->Nelements+mesh->totalHaloPairs; 
    mesh->gradientKernel(allNelements, mesh->o_vgeo, mesh->o_D, o_q, solver->o_grad);

    // TW NOTE WAS 2 !
    dfloat tau = 2.f*mesh->Nq*mesh->Nq; // 1/h factor built into kernel 
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
    ellipticParallelGatherScatterQuad2D(mesh, solver->ogs, o_Aq, o_Aq, dfloatString, "add");

  mesh->device.finish();
  occa::toc("AxKernel");
}

void ellipticScaledAdd(solver_t *solver, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b){

  mesh_t *mesh = solver->mesh;

  iint Ntotal = mesh->Nelements*mesh->Np;

  mesh->device.finish();
  occa::tic("scaledAddKernel");
  
  // b[n] = alpha*a[n] + beta*b[n] n\in [0,Ntotal)
  mesh->scaledAddKernel(Ntotal, alpha, o_a, beta, o_b);

  mesh->device.finish();
  occa::toc("scaledAddKernel");
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

  mesh->device.finish();
  occa::tic("weighted inner product2");

  if(strstr(options,"CONTINUOUS")||strstr(options, "PROJECT"))
    mesh->weightedInnerProduct2Kernel(Ntotal, o_w, o_a, o_b, o_tmp);
  else
    mesh->innerProductKernel(Ntotal, o_a, o_b, o_tmp);
  
  mesh->device.finish();
  occa::toc("weighted inner product2");
  
  o_tmp.copyTo(tmp);

  dfloat wab = 0;
  for(iint n=0;n<Nblock;++n){
    wab += tmp[n];
  }
      
  dfloat globalwab = 0;
  MPI_Allreduce(&wab, &globalwab, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

  return globalwab;
}

dfloat ellipticLocalInnerProduct(solver_t *solver,
         occa::memory &o_a,
         occa::memory &o_b){

  mesh_t *mesh = solver->mesh;
  dfloat *tmp = solver->tmp;
  iint Nblock = solver->Nblock;
  iint Ntotal = mesh->Nelements*mesh->Np;

  occa::memory &o_tmp = solver->o_tmp;
  
  mesh->device.finish();
  occa::tic("inner product2");

  mesh->innerProductKernel(Ntotal, o_a, o_b, o_tmp);
  
  mesh->device.finish();
  occa::toc("inner product2");
  
  o_tmp.copyTo(tmp);
  
  dfloat ab = 0;
  for(iint n=0;n<Nblock;++n)
    ab += tmp[n];

  return ab;
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


  if(strstr(options, "OAS")){
    
    ellipticStartHaloExchange2D(mesh, o_r, sendBuffer, recvBuffer);
    
    ellipticEndHaloExchange2D(mesh, o_r, recvBuffer);

    //    diagnostic(mesh->Np*mesh->Nelements, o_r, "o_r");

    // compute local precon on DEVICE
    mesh->device.finish();
    occa::tic("OASpreconKernel");
    if(strstr(options, "CONTINUOUS")||strstr(options, "PROJECT")) {

      precon->preconKernel(mesh->Nelements,
         precon->o_vmapPP,
         precon->o_faceNodesP,
         precon->o_oasForward,
         precon->o_oasDiagInvOp,
         precon->o_oasBack,
         o_r,
         o_zP);
      ellipticParallelGatherScatterQuad2D(mesh, precon->ogsP, o_zP, o_zP, dfloatString, "add");
    }
    else{
      
      precon->preconKernel(mesh->Nelements,
         mesh->o_vmapP,
         precon->o_faceNodesP,
         precon->o_oasForwardDg,
         precon->o_oasDiagInvOpDg,
         precon->o_oasBackDg,
         o_r,
         o_zP);
      ellipticParallelGatherScatterQuad2D(mesh, precon->ogsDg, o_zP, o_zP, dfloatString, "add");
    } 
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

      // halo exchange to make sure each vertex patch has available halo
      ellipticStartHaloExchange2D(mesh, o_r, sendBuffer, recvBuffer);
      
      ellipticEndHaloExchange2D(mesh, o_r, recvBuffer);
      
      mesh->device.finish();
      occa::tic("coarsenKernel");
      // Z1*Z1'*PL1*(Z1*z1) = (Z1*rL)  HMMM
      precon->coarsenKernel(mesh->Nelements, precon->o_coarseInvDegree, precon->o_V1, o_r, precon->o_r1);
      mesh->device.finish();
      occa::toc("coarsenKernel");

      if(strstr(options,"XXT")){
        precon->o_r1.copyTo(precon->r1); 
        xxtSolve(precon->z1, precon->xxt,precon->r1);
        precon->o_z1.copyFrom(precon->z1);
      }

      if(strstr(options,"GLOBALALMOND")){
        // should eliminate these copies
        mesh->device.finish();
        occa::tic("parAlmond");
        precon->o_r1.copyTo(precon->r1); 
        almondSolve(precon->z1, precon->parAlmond, precon->r1);
        precon->o_z1.copyFrom(precon->z1);
        mesh->device.finish();
        occa::toc("parAlmond");
      }      

      if(strstr(options,"LOCALALMOND")){
        // should eliminate these copies
        mesh->device.finish();
        occa::tic("Almond");
        precon->o_r1.copyTo(precon->r1); 
        almondSolve(precon->z1, precon->almond, precon->r1);
        precon->o_z1.copyFrom(precon->z1);
        mesh->device.finish();
        occa::toc("Almond");
      }

      mesh->device.finish();
      occa::tic("prolongateKernel");
      precon->prolongateKernel(mesh->Nelements, precon->o_V1, precon->o_z1, precon->o_ztmp);
      mesh->device.finish();
      occa::toc("prolongateKernel");

      // do we have to DG gatherscatter here 
      dfloat one = 1.;
      ellipticScaledAdd(solver, one, precon->o_ztmp, one, o_z);

      mesh->device.finish();
      occa::toc("coarseGrid");
    }
    
  } else if (strstr(options, "FULLALMOND")) {

    o_r.copyTo(precon->r1); 
    occa::tic("parALMOND");
    almondSolve(precon->z1, precon->parAlmond, precon->r1);
    occa::toc("parALMOND");
    o_z.copyFrom(precon->z1);

  } else if(strstr(options, "JACOBI")){

    mesh->device.finish();
    occa::tic("dotDivideKernel");   
    // Jacobi preconditioner
    iint Ntotal = mesh->Np*mesh->Nelements;
    mesh->dotDivideKernel(Ntotal, o_r, precon->o_diagA, o_z);
    mesh->device.finish();
    occa::toc("dotDivideKernel");   
  }
  else {// turn off preconditioner
    o_z.copyFrom(o_r); 
  }
}


int ellipticSolveQuad2D(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const char *options){

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
  occa::memory &o_zP = solver->o_zP;
  occa::memory &o_Ap = solver->o_Ap;
  occa::memory &o_Ax = solver->o_Ax;

  mesh->device.finish();
  occa::tic("PCG");

  // compute A*x
  ellipticOperator2D(solver, lambda, o_x, o_Ax, options);
  
  // subtract r = b - A*x
  ellipticScaledAdd(solver, -1.f, o_Ax, 1.f, o_r);

  // gather-scatter 
  if(strstr(options, "CONTINUOUS")||strstr(options, "PROJECT"))
    ellipticParallelGatherScatterQuad2D(mesh, solver->ogs, o_r, o_r, dfloatString, "add");
  
  mesh->device.finish();
  occa::tic("Preconditioner");
  if(strstr(options,"PCG")){

    // Precon^{-1} (b-A*x)
    ellipticPreconditioner2D(solver, o_r, o_zP, o_z, options); // r => rP => zP => z
    
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
  dfloat rdotr0 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_r, options);
  dfloat rdotz0 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_z, options);
  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;
  iint Niter = 0;
  dfloat alpha, beta;

  if(rank==0)
    printf("rdotr0 = %g, rdotz0 = %g\n", rdotr0, rdotz0);
  
  do{
    //diagnostic(Ntotal, o_p, "o_p");
    
    // A*p
    ellipticOperator2D(solver, lambda, o_p, o_Ap, options); 

    //diagnostic(Ntotal, o_p, "o_Ap");
    
    // dot(p,A*p)
    dfloat pAp = ellipticWeightedInnerProduct(solver, solver->o_invDegree,o_p, o_Ap, options);
    
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

    mesh->device.finish();
    occa::tic("Preconditioner");
    if(strstr(options,"PCG")){

      // z = Precon^{-1} r
      ellipticPreconditioner2D(solver, o_r, o_zP, o_z, options);

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
    mesh->device.finish();
    occa::toc("Preconditioner");

    // switch rdotr0 <= rdotr1
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