#include "ellipticQuad2D.h"

void ellipticOperator2D(solver_t *solver, dfloat lambda,
      occa::memory &o_q, occa::memory &o_Aq, const char *options){

  mesh_t *mesh = solver->mesh;

  occaTimerTic(mesh->device,"AxKernel");

  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;
  dfloat *gradSendBuffer = solver->gradSendBuffer;
  dfloat *gradRecvBuffer = solver->gradRecvBuffer;

  dfloat alpha = 0., alphaG = 0.;
  dlong Nblock = solver->Nblock;
  dfloat *tmp = solver->tmp;
  occa::memory &o_tmp = solver->o_tmp;

  if(strstr(options, "CONTINUOUS")){
    ogs_t *ogs = solver->mesh->ogs;

    if(solver->allNeumann)
      //solver->innerProductKernel(mesh->Nelements*mesh->Np, solver->o_invDegree,o_q, o_tmp);
      mesh->sumKernel(mesh->Nelements*mesh->Np, o_q, o_tmp);

    if(ogs->NhaloGather) {
      solver->partialAxKernel(solver->NglobalGatherElements, solver->o_globalGatherElementList,
                              mesh->o_ggeo, mesh->o_D, lambda, o_q, o_Aq);
      mesh->device.finish();
      mesh->device.setStream(solver->dataStream);
      mesh->gatherKernel(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherLocalIds, o_Aq, ogs->o_haloGatherTmp);
      ogs->o_haloGatherTmp.asyncCopyTo(ogs->haloGatherTmp);
      mesh->device.setStream(solver->defaultStream);
    }

    if(solver->NlocalGatherElements){
      solver->partialAxKernel(solver->NlocalGatherElements, solver->o_localGatherElementList,
                              mesh->o_ggeo, mesh->o_D, lambda, o_q, o_Aq);
    }

    if(solver->allNeumann) {
      o_tmp.copyTo(tmp);

      for(dlong n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
      alphaG *= solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
    }

    if(ogs->NnonHaloGather) 
      mesh->gatherScatterKernel(ogs->NnonHaloGather, ogs->o_nonHaloGatherOffsets, ogs->o_nonHaloGatherLocalIds, o_Aq);

    // C0 halo gather-scatter (on data stream)
    if(ogs->NhaloGather) {
      mesh->device.setStream(solver->dataStream);
      mesh->device.finish();

      // MPI based gather scatter using libgs
      gsParallelGatherScatter(ogs->haloGsh, ogs->haloGatherTmp, dfloatString, "add");

      // copy totally gather halo data back from HOST to DEVICE
      ogs->o_haloGatherTmp.asyncCopyFrom(ogs->haloGatherTmp);
    
      // do scatter back to local nodes
      mesh->scatterKernel(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherLocalIds, ogs->o_haloGatherTmp, o_Aq);
      mesh->device.setStream(solver->defaultStream);
    }

    if(solver->allNeumann) {
      mesh->addScalarKernel((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np, alphaG, o_Aq);
    }

    mesh->device.finish();    
    mesh->device.setStream(solver->dataStream);
    mesh->device.finish();    
    mesh->device.setStream(solver->defaultStream);

    //post-mask
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_Aq);

  } else if(strstr(options, "IPDG")) {
    dlong offset = 0;
    dfloat alpha = 0., alphaG =0.;
    dlong Nblock = solver->Nblock;
    dfloat *tmp = solver->tmp;
    occa::memory &o_tmp = solver->o_tmp;

    ellipticStartHaloExchange2D(solver, o_q, mesh->Np, sendBuffer, recvBuffer);

    solver->partialGradientKernel(mesh->Nelements,
                                  offset,
                                  mesh->o_vgeo,
                                  mesh->o_D,
                                  o_q,
                                  solver->o_grad);

    ellipticInterimHaloExchange2D(solver, o_q, mesh->Np, sendBuffer, recvBuffer);

    //Start the rank 1 augmentation if all BCs are Neumann
    //TODO this could probably be moved inside the Ax kernel for better performance
    if(solver->allNeumann)
      mesh->sumKernel(mesh->Nelements*mesh->Np, o_q, o_tmp);

    if(mesh->NinternalElements) {
      solver->partialIpdgKernel(mesh->NinternalElements,
                                mesh->o_internalElementIds,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                lambda,
                                solver->tau,
                                mesh->o_vgeo,
                                mesh->o_sgeo,
                                solver->o_EToB,
                                mesh->o_D,
                                solver->o_grad,
                                o_Aq);
    }

    if(solver->allNeumann) {
      o_tmp.copyTo(tmp);

      for(dlong n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
      alphaG *= solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
    }

    ellipticEndHaloExchange2D(solver, o_q, mesh->Np, recvBuffer);

    if(mesh->totalHaloPairs){
      offset = mesh->Nelements;
      solver->partialGradientKernel(mesh->totalHaloPairs,
                                    offset,
                                    mesh->o_vgeo,
                                    mesh->o_D,
                                    o_q,
                                    solver->o_grad);
    }

    if(mesh->NnotInternalElements) {
      solver->partialIpdgKernel(mesh->NnotInternalElements,
                                mesh->o_notInternalElementIds,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                lambda,
                                solver->tau,
                                mesh->o_vgeo,
                                mesh->o_sgeo,
                                solver->o_EToB,
                                mesh->o_D,
                                solver->o_grad,
                                o_Aq);
    }

    if(solver->allNeumann)
      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);
  }

  occaTimerToc(mesh->device,"AxKernel");
}

void ellipticScaledAdd(solver_t *solver, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b){

  mesh_t *mesh = solver->mesh;

  dlong Ntotal = mesh->Nelements*mesh->Np;

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
  dlong Nblock = solver->Nblock;
  dlong Ntotal = mesh->Nelements*mesh->Np;

  occa::memory &o_tmp = solver->o_tmp;

  occaTimerTic(mesh->device,"weighted inner product2");

  if(strstr(options,"CONTINUOUS"))
    solver->weightedInnerProduct2Kernel(Ntotal, o_w, o_a, o_b, o_tmp);
  else
    solver->innerProductKernel(Ntotal, o_a, o_b, o_tmp);

  occaTimerToc(mesh->device,"weighted inner product2");

  o_tmp.copyTo(tmp);

  dfloat wab = 0;
  for(dlong n=0;n<Nblock;++n){
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
  dlong Nblock = solver->Nblock;
  dlong Ntotal = mesh->Nelements*mesh->Np;

  occa::memory &o_tmp = solver->o_tmp;

  occaTimerTic(mesh->device,"inner product2");
  solver->innerProductKernel(Ntotal, o_a, o_b, o_tmp);
  occaTimerToc(mesh->device,"inner product2");

  o_tmp.copyTo(tmp);

  dfloat ab = 0;
  for(dlong n=0;n<Nblock;++n)
    ab += tmp[n];

  return ab;
}


int ellipticSolveQuad2D(solver_t *solver, dfloat lambda, dfloat tol,
                        occa::memory &o_r, occa::memory &o_x, const char *options){

  mesh_t *mesh = solver->mesh;

  int Niter = 0;
  int maxIter = 5000; 

  double start = 0.0, end =0.0;

  if(strstr(options,"VERBOSE")){
    mesh->device.finish();
    start = MPI_Wtime(); 
  }

  occaTimerTic(mesh->device,"Linear Solve");
  Niter = pcg (solver, options, lambda, o_r, o_x, tol, maxIter);
  occaTimerToc(mesh->device,"Linear Solve");

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(strstr(options,"VERBOSE")){
    mesh->device.finish();
    end = MPI_Wtime();
    double localElapsed = end-start;

    occa::printTimer();

    if(rank==0) printf("Solver converged in %d iters \n", Niter );

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    hlong   localDofs = (hlong) mesh->Np*mesh->Nelements;
    hlong   localElements = (hlong) mesh->Nelements;
    double globalElapsed;
    hlong   globalDofs;
    hlong   globalElements;

    MPI_Reduce(&localElapsed, &globalElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_HLONG,   MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce(&localElements,&globalElements,1, MPI_HLONG,   MPI_SUM, 0, MPI_COMM_WORLD );

    if (rank==0){
      printf("%02d %02d "hlongFormat" "hlongFormat" %d %17.15lg %3.5g \t [ RANKS N NELEMENTS DOFS ITERATIONS ELAPSEDTIME PRECONMEMORY] \n",
             size, mesh->N, globalElements, globalDofs, Niter, globalElapsed, solver->precon->preconBytes/(1E9));
    }
  }
  return Niter;
}