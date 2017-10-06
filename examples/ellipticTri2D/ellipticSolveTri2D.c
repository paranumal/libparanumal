#include "ellipticTri2D.h"

void ellipticOperator2D(solver_t *solver, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *options){

  mesh_t *mesh = solver->mesh;

  occaTimerTic(mesh->device,"AxKernel");

  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;
  dfloat *gradSendBuffer = solver->gradSendBuffer;
  dfloat *gradRecvBuffer = solver->gradRecvBuffer;

  if(strstr(options, "CONTINUOUS")){
    ogs_t *nonHalo = solver->nonHalo;
    ogs_t *halo = solver->halo;

    if(halo->Ngather) {
      solver->partialAxKernel(solver->NglobalGatherElements, solver->o_globalGatherElementList,
                              mesh->o_ggeo, mesh->o_SrrT, mesh->o_SrsT, mesh->o_SsrT, mesh->o_SssT,
                              mesh->o_MM, lambda, o_q, o_Aq);
      mesh->gatherKernel(halo->Ngather, halo->o_gatherOffsets, halo->o_gatherLocalIds, o_Aq, halo->o_gatherTmp);
      halo->o_gatherTmp.copyTo(halo->gatherTmp);
    }
    if(nonHalo->Ngather)
      solver->partialAxKernel(solver->NlocalGatherElements, solver->o_localGatherElementList,
                              mesh->o_ggeo, mesh->o_SrrT, mesh->o_SrsT, mesh->o_SsrT, mesh->o_SssT,
                              mesh->o_MM, lambda, o_q, o_Aq);

    // C0 halo gather-scatter (on data stream)
    if(halo->Ngather) {
      occa::streamTag tag;

      // MPI based gather scatter using libgs
      gsParallelGatherScatter(halo->gatherGsh, halo->gatherTmp, dfloatString, "add");

      // copy totally gather halo data back from HOST to DEVICE
      mesh->device.setStream(solver->dataStream);
      halo->o_gatherTmp.asyncCopyFrom(halo->gatherTmp);

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
    if(nonHalo->Ngather)
      mesh->gatherScatterKernel(nonHalo->Ngather, nonHalo->o_gatherOffsets, nonHalo->o_gatherLocalIds, o_Aq);


  } else{

    iint offset = 0;
    dfloat alpha = 0., alphaG =0.;
    iint Nblock = solver->Nblock;
    dfloat *tmp = solver->tmp;
    occa::memory &o_tmp = solver->o_tmp;

    dfloat BRPenalty = 1.0;

    ellipticStartHaloExchange2D(solver, o_q, mesh->Np, sendBuffer, recvBuffer);
    ellipticInterimHaloExchange2D(solver, o_q, mesh->Np, sendBuffer, recvBuffer);    

#if USE_BERN
    solver->BRGradientVolumeKernel(mesh->Nelements,
                                          mesh->o_vgeo,
                                          mesh->o_D1ids,
                                          mesh->o_D2ids,
                                          mesh->o_D3ids,
                                          mesh->o_Dvals,
                                          o_q,
                                          solver->o_grad);
#else
    solver->BRGradientVolumeKernel(mesh->Nelements,
                                          mesh->o_vgeo,
                                          mesh->o_DrT,
                                          mesh->o_DsT,
                                          o_q,
                                          solver->o_grad);
#endif

    ellipticEndHaloExchange2D(solver, o_q, mesh->Np, recvBuffer);
    
#if USE_BERN
    solver->BRGradientSurfaceKernel(mesh->Nelements,
                               mesh->o_vmapM,
                               mesh->o_vmapP,
                               mesh->o_sgeo,
                               solver->o_EToB,
                               mesh->o_L0vals,
                               mesh->o_ELids,
                               mesh->o_ELvals,
                               o_q,
                               solver->o_grad);    
#else
    solver->BRGradientSurfaceKernel(mesh->Nelements,
                               mesh->o_vmapM,
                               mesh->o_vmapP,
                               mesh->o_sgeo,
                               solver->o_EToB,
                               mesh->o_LIFTT,
                               o_q,
                               solver->o_grad);
#endif

    //Start the rank 1 augmentation if all BCs are Neumann
    //TODO this could probably be moved inside the Ax kernel for better performance
    if(solver->allNeumann)  {
      mesh->sumKernel(mesh->Nelements*mesh->Np, o_q, o_tmp);
    
      o_tmp.copyTo(tmp);

      for(iint n=0;n<Nblock;++n)
        alpha += tmp[n];
      
      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
      alphaG *= solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
    }
    
    ellipticStartHaloExchange2D(solver, solver->o_grad, 2*mesh->Np, gradSendBuffer, gradRecvBuffer);
    ellipticInterimHaloExchange2D(solver, solver->o_grad, 2*mesh->Np, gradSendBuffer, gradRecvBuffer);    

#if USE_BERN
    solver->BRDivergenceVolumeKernel(mesh->Nelements,
                                          mesh->o_vgeo,
                                          mesh->o_D1ids,
                                          mesh->o_D2ids,
                                          mesh->o_D3ids,
                                          mesh->o_Dvals,
                                          solver->o_grad,
                                          o_Aq);
#else
    solver->BRDivergenceVolumeKernel(mesh->Nelements,
                                          mesh->o_vgeo,
                                          mesh->o_DrT,
                                          mesh->o_DsT,
                                          solver->o_grad,
                                          o_Aq);
#endif

    ellipticEndHaloExchange2D(solver, solver->o_grad, 2*mesh->Np, gradRecvBuffer);
    
#if USE_BERN
    solver->BRDivergenceSurfaceKernel(mesh->Nelements,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                lambda,
                                BRPenalty,
                                mesh->o_vgeo,
                                mesh->o_sgeo,
                                solver->o_EToB,
                                mesh->o_L0vals,
                                mesh->o_ELids,
                                mesh->o_ELvals,
                                mesh->o_BBMM,
                                o_q,
                                solver->o_grad,
                                o_Aq);
#else
    solver->BRDivergenceSurfaceKernel(mesh->Nelements,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                lambda,
                                BRPenalty,
                                mesh->o_vgeo,
                                mesh->o_sgeo,
                                solver->o_EToB,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                o_q,
                                solver->o_grad,
                                o_Aq);
#endif

    if(solver->allNeumann) 
      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);
  }

  if(strstr(options, "CONTINUOUS"))
    // parallel gather scatter
    ellipticParallelGatherScatterTri2D(mesh, solver->ogs, o_Aq, o_Aq, dfloatString, "add");

  occaTimerToc(mesh->device,"AxKernel");
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

  if(strstr(options,"CONTINUOUS"))
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

int ellipticSolveTri2D(solver_t *solver, dfloat lambda, dfloat tol, 
                        occa::memory &o_r, occa::memory &o_x, const char *options){

  mesh_t *mesh = solver->mesh;

  iint rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // placeholder conjugate gradient:
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method

  // placeholder preconditioned conjugate gradient
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method

  occa::memory &o_p  = solver->o_p;
  occa::memory &o_z  = solver->o_z;
  occa::memory &o_Ap = solver->o_Ap;
  occa::memory &o_Ax = solver->o_Ax;

  occaTimerTic(mesh->device,"PCG");

  //start timer
  mesh->device.finish();
  MPI_Barrier(MPI_COMM_WORLD);
  double tic = MPI_Wtime();

  // gather-scatter
  if(strstr(options, "CONTINUOUS"))
    ellipticParallelGatherScatterTri2D(mesh, solver->ogs, o_r, o_r, dfloatString, "add");

  // compute A*x
  ellipticOperator2D(solver, lambda, o_x, solver->o_Ax, options);

  // subtract r = b - A*x
  ellipticScaledAdd(solver, -1.f, o_Ax, 1.f, o_r);

  dfloat rdotr0 = ellipticWeightedInnerProduct(solver, solver->o_invDegree, o_r, o_r, options);
  
  // 
  dfloat n2b     = rdotr0>(1e-12*1e-12) ? rdotr0:1.0; 
  dfloat ABS_TOL = 1e-10*1e-10; // absolute tolerance 10^-10
  dfloat REL_TOL = tol*tol*n2b; // 
  //
  //dfloat TOL     = ABS_TOL>REL_TOL ? ABS_TOL:REL_TOL; 


  dfloat *Ax = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  dfloat *x = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  for (int i=0;i<mesh->Nelements*mesh->Np;i++) {
    x[i] = 1.;
    o_x.copyFrom(x);
    ellipticOperator2D(solver, lambda, o_x, o_Ax, options);
    o_Ax.copyTo(Ax);
    for (int j =0;j<mesh->Nelements*mesh->Np;j++) {
      printf("%4.2f \t", Ax[j]);
    }
    printf("\n");
    x[i] = 0.;
  }

  dfloat TOL     = tol*tol; 
  
  dfloat rdotz0 = 0;
  iint Niter = 0;
  //sanity check
  if (rdotr0<=(TOL)) {
   // printf("iter=0 norm(r) = %g\n", sqrt(rdotr0));
    occaTimerToc(mesh->device,"PCG");
    return 0;
  }


  occaTimerTic(mesh->device,"Preconditioner");
  if(strstr(options,"PCG")){
    // Precon^{-1} (b-A*x)
    ellipticPreconditioner2D(solver, lambda, o_r, o_z, options);

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

  if((rank==0)&&strstr(options,"VERBOSE"))
    printf("rdotr0 = %g, rdotz0 = %g\n", rdotr0, rdotz0);

  while(rdotr0>(TOL)){

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

    if(rdotr1 < TOL) {
      rdotr0 = rdotr1;
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
     
    //if((rank==0)&&(strstr(options,"VERBOSE")))
    //  printf("iter=%05d pAp = %g norm(r) = %g\n", Niter, pAp, sqrt(rdotr0)/sqrt(n2b));

    ++Niter;

  }


   //printf("iter=%05d pAp = %g norm(r) = %g relnorm(r) = %g\n", Niter, pAp, sqrt(rdotr0), sqrt(rdotr0)/sqrt(n2b));
  //if((rank==0)&&strstr(options,"VERBOSE"))
  //  printf("iter=%05d pAp = %g norm(r) = %g\n", Niter, pAp, sqrt(rdotr0));


  if(strstr(options,"VERBOSE")){
    mesh->device.finish();
    double toc = MPI_Wtime();
    double localElapsed = toc-tic;

    iint size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    iint   localDofs = mesh->Np*mesh->Nelements;
    iint   localElements = mesh->Nelements;
    double globalElapsed;
    iint   globalDofs;
    iint   globalElements;

    MPI_Reduce(&localElapsed, &globalElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_IINT,   MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce(&localElements,&globalElements,1, MPI_IINT,   MPI_SUM, 0, MPI_COMM_WORLD );

    if (rank==0)
      printf("%02d %02d %d %d %d %17.15lg %3.5g \t [ RANKS N NELEMENTS DOFS ITERATIONS ELAPSEDTIME PRECONMEMORY] \n",
           size, mesh->N, globalElements, globalDofs, Niter, globalElapsed, solver->precon->preconBytes/(1E9));
  }

  occaTimerToc(mesh->device,"PCG");

 //occa::printTimer();

  return Niter;

}
