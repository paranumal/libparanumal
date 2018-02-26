#include "ellipticTri2D.h"

void ellipticOperator2D(solver_t *solver, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *options){

  mesh_t *mesh = solver->mesh;

  occaTimerTic(mesh->device,"AxKernel");

  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;
  dfloat *gradSendBuffer = solver->gradSendBuffer;
  dfloat *gradRecvBuffer = solver->gradRecvBuffer;

  dfloat alpha = 0., alphaG =0.;
  int Nblock = solver->Nblock;
  dfloat *tmp = solver->tmp;
  occa::memory &o_tmp = solver->o_tmp;

  if(strstr(options, "CONTINUOUS")){
    ogs_t *ogs = solver->mesh->ogs;

    if(solver->allNeumann)
      //solver->innerProductKernel(mesh->Nelements*mesh->Np, solver->o_invDegree,o_q, o_tmp);
      mesh->sumKernel(mesh->Nelements*mesh->Np, o_q, o_tmp);

    if(ogs->NhaloGather) {
      if (strstr(options, "SPARSE")){
        solver->partialAxKernel(solver->NglobalGatherElements, solver->o_globalGatherElementList,
            mesh->o_ggeo, mesh->o_sparseStackedNZ, mesh->o_sparseSrrT, mesh->o_sparseSrsT, mesh->o_sparseSssT,
            mesh->o_MM, lambda, o_q, o_Aq);
      } else{
        solver->partialAxKernel(solver->NglobalGatherElements, solver->o_globalGatherElementList,
            mesh->o_ggeo, mesh->o_SrrT, mesh->o_SrsT, mesh->o_SsrT, mesh->o_SssT,
            mesh->o_MM, lambda, o_q, o_Aq);
      }
      mesh->device.finish();
      mesh->device.setStream(solver->dataStream);
      mesh->gatherKernel(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherLocalIds, o_Aq, ogs->o_haloGatherTmp);
      ogs->o_haloGatherTmp.asyncCopyTo(ogs->haloGatherTmp);
      mesh->device.setStream(solver->defaultStream);
    }
    if(solver->NlocalGatherElements){
      if (strstr(options, "SPARSE")){
        solver->partialAxKernel(solver->NlocalGatherElements, solver->o_localGatherElementList,
            mesh->o_ggeo, mesh->o_sparseStackedNZ, mesh->o_sparseSrrT, mesh->o_sparseSrsT, mesh->o_sparseSssT,
            mesh->o_MM, lambda, o_q, o_Aq);
      } else {
        solver->partialAxKernel(solver->NlocalGatherElements, solver->o_localGatherElementList,
            mesh->o_ggeo, mesh->o_SrrT, mesh->o_SrsT, mesh->o_SsrT, mesh->o_SssT,
            mesh->o_MM, lambda, o_q, o_Aq);
      }
    }
    if(solver->allNeumann) {
      o_tmp.copyTo(tmp);

      for(int n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
      alphaG *= solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
    }

    // finalize gather using local and global contributions
    if (strstr(options,"SPARSE")) solver->dotMultiplyKernel(mesh->Np*mesh->Nelements, o_Aq, mesh->o_mapSgn, o_Aq);
    if(ogs->NnonHaloGather) mesh->gatherScatterKernel(ogs->NnonHaloGather, ogs->o_nonHaloGatherOffsets, ogs->o_nonHaloGatherLocalIds, o_Aq);
    if (strstr(options,"SPARSE")) solver->dotMultiplyKernel(mesh->Np*mesh->Nelements, o_Aq, mesh->o_mapSgn, o_Aq);

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
      //dfloat one = 1.f;
      //solver->scaledAddKernel(mesh->Nelements*mesh->Np, alphaG, solver->o_invDegree, one, o_Aq);
    }

    mesh->device.finish();    
    mesh->device.setStream(solver->dataStream);
    mesh->device.finish();    
    mesh->device.setStream(solver->defaultStream);

    //post-mask
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_Aq);

  } else if(strstr(options, "IPDG")) {
    int offset = 0;
    dfloat alpha = 0., alphaG =0.;
    int Nblock = solver->Nblock;
    dfloat *tmp = solver->tmp;
    occa::memory &o_tmp = solver->o_tmp;

    ellipticStartHaloExchange2D(solver, o_q, mesh->Np, sendBuffer, recvBuffer);

    if(strstr(options, "NODAL")) {
      solver->partialGradientKernel(mesh->Nelements,
          offset,
          mesh->o_vgeo,
          mesh->o_DrT,
          mesh->o_DsT,
          o_q,
          solver->o_grad);
    } else if(strstr(options, "BERN")) {
      solver->partialGradientKernel(mesh->Nelements,
          offset,
          mesh->o_vgeo,
          mesh->o_D1ids,
          mesh->o_D2ids,
          mesh->o_D3ids,
          mesh->o_Dvals,
          o_q,
          solver->o_grad);
    }

    ellipticInterimHaloExchange2D(solver, o_q, mesh->Np, sendBuffer, recvBuffer);

    //Start the rank 1 augmentation if all BCs are Neumann
    //TODO this could probably be moved inside the Ax kernel for better performance
    if(solver->allNeumann)
      mesh->sumKernel(mesh->Nelements*mesh->Np, o_q, o_tmp);

    if(mesh->NinternalElements) {
      if(strstr(options, "NODAL")) {
        solver->partialIpdgKernel(mesh->NinternalElements,
            mesh->o_internalElementIds,
            mesh->o_vmapM,
            mesh->o_vmapP,
            lambda,
            solver->tau,
            mesh->o_vgeo,
            mesh->o_sgeo,
            solver->o_EToB,
            mesh->o_DrT,
            mesh->o_DsT,
            mesh->o_LIFTT,
            mesh->o_MM,
            solver->o_grad,
            o_Aq);
      } else if(strstr(options, "BERN")) {
        solver->partialIpdgKernel(mesh->NinternalElements,
            mesh->o_internalElementIds,
            mesh->o_vmapM,
            mesh->o_vmapP,
            lambda,
            solver->tau,
            mesh->o_vgeo,
            mesh->o_sgeo,
            solver->o_EToB,
            mesh->o_D1ids,
            mesh->o_D2ids,
            mesh->o_D3ids,
            mesh->o_Dvals,
            mesh->o_L0vals,
            mesh->o_ELids,
            mesh->o_ELvals,
            mesh->o_BBMM,
            solver->o_grad,
            o_Aq);
      }
    }

    if(solver->allNeumann) {
      o_tmp.copyTo(tmp);

      for(int n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
      alphaG *= solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
    }

    ellipticEndHaloExchange2D(solver, o_q, mesh->Np, recvBuffer);

    if(mesh->totalHaloPairs){
      offset = mesh->Nelements;
      if(strstr(options, "NODAL")) {
        solver->partialGradientKernel(mesh->totalHaloPairs,
            offset,
            mesh->o_vgeo,
            mesh->o_DrT,
            mesh->o_DsT,
            o_q,
            solver->o_grad);
      } else if(strstr(options, "BERN")) {
        solver->partialGradientKernel(mesh->totalHaloPairs,
            offset,
            mesh->o_vgeo,
            mesh->o_D1ids,
            mesh->o_D2ids,
            mesh->o_D3ids,
            mesh->o_Dvals,
            o_q,
            solver->o_grad);
      }
    }

    if(mesh->NnotInternalElements) {
      if(strstr(options, "NODAL")) {
        solver->partialIpdgKernel(mesh->NnotInternalElements,
            mesh->o_notInternalElementIds,
            mesh->o_vmapM,
            mesh->o_vmapP,
            lambda,
            solver->tau,
            mesh->o_vgeo,
            mesh->o_sgeo,
            solver->o_EToB,
            mesh->o_DrT,
            mesh->o_DsT,
            mesh->o_LIFTT,
            mesh->o_MM,
            solver->o_grad,
            o_Aq);
      } else if(strstr(options, "BERN")) {
        solver->partialIpdgKernel(mesh->NnotInternalElements,
            mesh->o_notInternalElementIds,
            mesh->o_vmapM,
            mesh->o_vmapP,
            lambda,
            solver->tau,
            mesh->o_vgeo,
            mesh->o_sgeo,
            solver->o_EToB,
            mesh->o_D1ids,
            mesh->o_D2ids,
            mesh->o_D3ids,
            mesh->o_Dvals,
            mesh->o_L0vals,
            mesh->o_ELids,
            mesh->o_ELvals,
            mesh->o_BBMM,
            solver->o_grad,
            o_Aq);
      }
    }

    if(solver->allNeumann)
      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);

  } else if (strstr(options, "BRDG")){

    int offset = 0;
    dfloat alpha = 0., alphaG =0.;
    int Nblock = solver->Nblock;
    dfloat *tmp = solver->tmp;
    occa::memory &o_tmp = solver->o_tmp;

    ellipticStartHaloExchange2D(solver, o_q, mesh->Np, sendBuffer, recvBuffer);
    ellipticInterimHaloExchange2D(solver, o_q, mesh->Np, sendBuffer, recvBuffer);

    if(strstr(options, "NODAL")) {
      solver->BRGradientVolumeKernel(mesh->Nelements,
          mesh->o_vgeo,
          mesh->o_DrT,
          mesh->o_DsT,
          o_q,
          solver->o_grad);
    } else if(strstr(options, "BERN")) {
      solver->BRGradientVolumeKernel(mesh->Nelements,
          mesh->o_vgeo,
          mesh->o_D1ids,
          mesh->o_D2ids,
          mesh->o_D3ids,
          mesh->o_Dvals,
          o_q,
          solver->o_grad);
    }

    ellipticEndHaloExchange2D(solver, o_q, mesh->Np, recvBuffer);

    if(strstr(options, "NODAL")) {
      solver->BRGradientSurfaceKernel(mesh->Nelements,
          mesh->o_vmapM,
          mesh->o_vmapP,
          mesh->o_sgeo,
          solver->o_EToB,
          mesh->o_LIFTT,
          o_q,
          solver->o_grad);
    } else if(strstr(options, "BERN")) {
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
    }

    //Start the rank 1 augmentation if all BCs are Neumann
    //TODO this could probably be moved inside the Ax kernel for better performance
    if(solver->allNeumann)  {
      mesh->sumKernel(mesh->Nelements*mesh->Np, o_q, o_tmp);

      o_tmp.copyTo(tmp);

      for(int n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
      alphaG *= solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
    }

    ellipticStartHaloExchange2D(solver, solver->o_grad, 2*mesh->Np, gradSendBuffer, gradRecvBuffer);
    ellipticInterimHaloExchange2D(solver, solver->o_grad, 2*mesh->Np, gradSendBuffer, gradRecvBuffer);

    if(strstr(options, "NODAL")) {
      solver->BRDivergenceVolumeKernel(mesh->Nelements,
          mesh->o_vgeo,
          mesh->o_DrT,
          mesh->o_DsT,
          solver->o_grad,
          o_Aq);
    } else if(strstr(options, "BERN")) {
      solver->BRDivergenceVolumeKernel(mesh->Nelements,
          mesh->o_vgeo,
          mesh->o_D1ids,
          mesh->o_D2ids,
          mesh->o_D3ids,
          mesh->o_Dvals,
          solver->o_grad,
          o_Aq);
    }

    ellipticEndHaloExchange2D(solver, solver->o_grad, 2*mesh->Np, gradRecvBuffer);

    if(strstr(options, "NODAL")) {
      solver->BRDivergenceSurfaceKernel(mesh->Nelements,
          mesh->o_vmapM,
          mesh->o_vmapP,
          lambda,
          solver->tau,
          mesh->o_vgeo,
          mesh->o_sgeo,
          solver->o_EToB,
          mesh->o_LIFTT,
          mesh->o_MM,
          o_q,
          solver->o_grad,
          o_Aq);
    } else if(strstr(options, "BERN")) {
      solver->BRDivergenceSurfaceKernel(mesh->Nelements,
          mesh->o_vmapM,
          mesh->o_vmapP,
          lambda,
          solver->tau,
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
    }

    if(solver->allNeumann)
      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);
  }

  occaTimerToc(mesh->device,"AxKernel");
}

dfloat ellipticScaledAdd(solver_t *solver, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b){

  mesh_t *mesh = solver->mesh;

  int Ntotal = mesh->Nelements*mesh->Np;

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
  int Nblock = solver->Nblock;
  int Ntotal = mesh->Nelements*mesh->Np;

  occa::memory &o_tmp = solver->o_tmp;

  occaTimerTic(mesh->device,"weighted inner product2");

  if(strstr(options,"CONTINUOUS"))
    solver->weightedInnerProduct2Kernel(Ntotal, o_w, o_a, o_b, o_tmp);
  else
    solver->innerProductKernel(Ntotal, o_a, o_b, o_tmp);

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

dfloat ellipticInnerProduct(solver_t *solver,
    occa::memory &o_a,
    occa::memory &o_b){


  mesh_t *mesh = solver->mesh;
  dfloat *tmp = solver->tmp;
  int Nblock = solver->Nblock;
  int Ntotal = mesh->Nelements*mesh->Np;

  occa::memory &o_tmp = solver->o_tmp;

  occaTimerTic(mesh->device,"inner product");
  solver->innerProductKernel(Ntotal, o_a, o_b, o_tmp);
  occaTimerToc(mesh->device,"inner product");

  o_tmp.copyTo(tmp);

  dfloat ab = 0;
  for(int n=0;n<Nblock;++n){
    ab += tmp[n];
  }

  dfloat globalab = 0;
  MPI_Allreduce(&ab, &globalab, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

  return globalab;
}

int ellipticSolveTri2D(solver_t *solver, dfloat lambda, dfloat tol,
    occa::memory &o_r, occa::memory &o_x, const char *options){

  mesh2D *mesh = solver->mesh;

/*
  dfloat *t = (dfloat *) calloc(mesh->Np*mesh->Nelements,sizeof(dfloat));
  dfloat *Pt = (dfloat *) calloc(mesh->Np*mesh->Nelements,sizeof(dfloat));
  occa::memory o_t = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat),t);
  occa::memory o_Pt = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat),t);
  
  char fname[BUFSIZ];
  sprintf(fname, "Precon.m");
  FILE *fp = fopen(fname, "w");

  fprintf(fp, "P = [");
  agmgLevel **levels = solver->precon->parAlmond->levels;
  for (int n=0;n<mesh->Np*mesh->Nelements;n++) {
    for (int m=0;m<mesh->Np*mesh->Nelements;m++) t[m] =0;
    o_Pt.copyFrom(t);

    //if (mesh->mask[n] == 0) {
      //for (int m=0;m<mesh->Np*mesh->Nelements;m++) Pt[m] =0;
    //} else {
      t[n] = 1;

      //gsParallelGatherScatter(solver->hostGsh, t, dfloatString, "add");

      o_t.copyFrom(t);
      //ellipticPreconditioner2D(solver, lambda, o_t, o_Pt, options);
      //levels[0]->device_smooth(levels[0]->smoothArgs, o_t, o_Pt, true);
      //ellipticOperator2D(solver, lambda, o_t, o_Pt, options);

      //levels[1]->device_coarsen(levels[1]->coarsenArgs, o_t, levels[1]->o_Srhs);
      levels[0]->device_gather (levels[0]->gatherArgs,  o_t, levels[0]->o_rhs);
      levels[0]->device_scatter(levels[0]->scatterArgs,  levels[0]->o_rhs, o_Pt);
      //levels[1]->device_prolongate(levels[1]->prolongateArgs, levels[1]->o_Srhs, o_Pt);

      o_Pt.copyTo(Pt);
    //}

    for (int m=0;m<mesh->Np*mesh->Nelements;m++) 
      fprintf(fp, "%.10e,", Pt[m]);

    fprintf(fp, "\n");
  }
  fprintf(fp, "];\n");
  fclose(fp);
*/

  int Niter;
  int maxIter = 5000; 

  double start, end;

  if(strstr(options,"VERBOSE")){
    mesh->device.finish();
    start = MPI_Wtime(); 
  }

  occaTimerTic(mesh->device,"Linear Solve");
  if(strstr(options, "GMRES")) {
    Niter = pgmresm  (solver, options, lambda, o_r, o_x, tol, maxIter);
  } else if(strstr(options, "BiCGStab")) {
    Niter = pbicgstab(solver, options, lambda, o_r, o_x, tol, maxIter);
  } else if(strstr(options, "CG")) {
    Niter = pcg      (solver, options, lambda, o_r, o_x, tol, maxIter);
  }
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

    int   localDofs = mesh->Np*mesh->Nelements;
    int   localElements = mesh->Nelements;
    double globalElapsed;
    int   globalDofs;
    int   globalElements;

    MPI_Reduce(&localElapsed, &globalElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_INT,   MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce(&localElements,&globalElements,1, MPI_INT,   MPI_SUM, 0, MPI_COMM_WORLD );

    if (rank==0){
      printf("%02d %02d %d %d %d %17.15lg %3.5g \t [ RANKS N NELEMENTS DOFS ITERATIONS ELAPSEDTIME PRECONMEMORY] \n",
             size, mesh->N, globalElements, globalDofs, Niter, globalElapsed, solver->precon->preconBytes/(1E9));
    }
  }
  return Niter;

}
