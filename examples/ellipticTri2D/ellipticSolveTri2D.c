#include "ellipticTri2D.h"
//this is bad coding practice, added for quick testing
dfloat timeAx;
int NAx;
//end
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
    ogs_t *nonHalo = solver->nonHalo;
    ogs_t *halo = solver->halo;

    if(solver->allNeumann)
      //solver->innerProductKernel(mesh->Nelements*mesh->Np, solver->o_invDegree,o_q, o_tmp);
      mesh->sumKernel(mesh->Nelements*mesh->Np, o_q, o_tmp);

    if(halo->Ngather) {
      occa::streamTag startAxtime = mesh->device.tagStream();

      if (strstr(options, "SPARSE")){
        solver->partialAxKernel(solver->NglobalGatherElements, solver->o_globalGatherElementList,
            mesh->o_ggeo, mesh->o_SrrT, mesh->o_SrsT, mesh->o_IndTchar,mesh->o_SssT,
            mesh->o_MM, lambda, o_q, o_Aq);
      }
      else{
        solver->partialAxKernel(solver->NglobalGatherElements, solver->o_globalGatherElementList,
            mesh->o_ggeo, mesh->o_SrrT, mesh->o_SrsT, mesh->o_SsrT, mesh->o_SssT,
            mesh->o_MM, lambda, o_q, o_Aq);
      }
      occa::streamTag stopAxtime = mesh->device.tagStream();
      NAx++;
      timeAx += mesh->device.timeBetween(startAxtime, stopAxtime);
      mesh->gatherKernel(halo->Ngather, halo->o_gatherOffsets, halo->o_gatherLocalIds, o_Aq, halo->o_gatherTmp);
      halo->o_gatherTmp.copyTo(halo->gatherTmp);
    }
    if(nonHalo->Ngather){
      occa::streamTag startAxtime = mesh->device.tagStream();
      if (strstr(options, "SPARSE")){
        solver->partialAxKernel(solver->NlocalGatherElements, solver->o_localGatherElementList,
            mesh->o_ggeo, mesh->o_SrrT, mesh->o_SrsT, mesh->o_IndTchar,mesh->o_SssT,
            mesh->o_MM, lambda, o_q, o_Aq);
      }
      else {
        solver->partialAxKernel(solver->NlocalGatherElements, solver->o_localGatherElementList,
            mesh->o_ggeo, mesh->o_SrrT, mesh->o_SrsT, mesh->o_SsrT, mesh->o_SssT,
            mesh->o_MM, lambda, o_q, o_Aq);
      }
      occa::streamTag stopAxtime = mesh->device.tagStream();
      NAx++;    
      timeAx += mesh->device.timeBetween(startAxtime, stopAxtime);

    }
    if(solver->allNeumann) {
      o_tmp.copyTo(tmp);

      for(int n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
      alphaG *= solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
    }

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

    if(solver->allNeumann) {
      mesh->addScalarKernel((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np, alphaG, o_Aq);
      //dfloat one = 1.f;
      //solver->scaledAddKernel(mesh->Nelements*mesh->Np, alphaG, solver->o_invDegree, one, o_Aq);
    }

    // finalize gather using local and global contributions
    mesh->device.setStream(solver->defaultStream);
    if(nonHalo->Ngather)
      mesh->gatherScatterKernel(nonHalo->Ngather, nonHalo->o_gatherOffsets, nonHalo->o_gatherLocalIds, o_Aq);


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
  // o_Aq.copyFrom(o_q);

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
    occa::memory &o_r, occa::memory &o_x, const char *options, int NblockV, int NnodesV){

  mesh2D *mesh = solver->mesh;
  printf("N=%d \n", mesh->N);
  //KS for testing
  /*  int Nbytes = mesh->Np*(6+mesh->Np*5);
      int flops = mesh->Np* (mesh->Np *10+15); */
  //TW testing
  //  int Nbytes = mesh->Np*(6+mesh->Np);//assume the operators are cached
  int Nbytes = sizeof(dfloat)*mesh->Np*(6+0*mesh->Np);
  //replace with mesh->Np*(6+4*mesh->Np); if they are NOT`
  //TW v1:
  //  int flops = mesh->Np*(mesh->Np*8+7);
  //TW v2
  int flops = mesh->Np*(mesh->Np*6+5);



  Nbytes /= 2; 
  printf("copyint %d elements \n", Nbytes*mesh->Nelements);  
  occa::memory o_foo = mesh->device.malloc(Nbytes*mesh->Nelements);
  occa::memory o_bah = mesh->device.malloc(Nbytes*mesh->Nelements);

  mesh->device.finish();

  occa::streamTag startCopy = mesh->device.tagStream();
  int Ntrials = 10;        
  for(int it=0;it<Ntrials;++it){
    o_bah.copyTo(o_foo);
  }
  occa::streamTag endCopy = mesh->device.tagStream();
  double  copyElapsed = mesh->device.timeBetween(startCopy, endCopy);
printf("copy elapsed %f ", copyElapsed);

  timeAx = 0.0f;


  // gather-scatter
  if(strstr(options, "CONTINUOUS"))
    ellipticParallelGatherScatterTri2D(mesh, solver->ogs, o_r, o_r, dfloatString, "add");

  int Niter;
  int maxIter = 500; 

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
  printf("time per Ax kernel %16.17f (executed %d times) \n", timeAx/NAx, NAx);
  double kernelElapsed = timeAx/NAx;


  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(strstr(options,"VERBOSE")){
    mesh->device.finish();
    end = MPI_Wtime();
    double localElapsed = end-start;

    occa::printTimer();

    printf("Solver converged in %d iters \n", Niter );

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int   localDofs = mesh->Np*mesh->Nelements;
    int   localElements = mesh->Nelements;
    double globalElapsed;
    int   globalDofs;
    int   globalElements;
    double globalCopyElapsed;

    MPI_Reduce(&localElapsed, &globalElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_int,   MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce(&localElements,&globalElements,1, MPI_int,   MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce(&copyElapsed,&globalCopyElapsed,1, MPI_DOUBLE,   MPI_MAX, 0, MPI_COMM_WORLD );



    if (rank==0){

      double copyBandwidth   = mesh->Nelements*((Nbytes*Ntrials*2.)/(1e9*globalCopyElapsed));
      double kernelBandwidth = mesh->Nelements*((Nbytes*2)/(1e9*kernelElapsed));
      double kernelGFLOPS = mesh->Nelements*flops/(1e9*kernelElapsed);  
      //globalElements*flops*iterations/(1024*1024*1024.*globalElapsed);      
      double roofline = mesh->Nelements*flops*Ntrials/(1e9*globalCopyElapsed);

      // double roofline = mesh->Nelements*flops*Ntrials/(1e9*copyElapsed);
      //( BWfromcopy512(N)*W)/(8*D);
      // double roofline = copyBandwidth*flops/(Nbytes*2*8);
      printf("copy BW %16.17g achieved BW %16.17g\n", copyBandwidth, kernelBandwidth);
      printf("ROOFLINE %16.17g \n", roofline);
      printf("GFLOPS %16.17f \n", kernelGFLOPS);
      printf("PARAMETERS %d %d %16.17f \n ", NblockV, NnodesV, kernelGFLOPS );
      printf("%02d %02d %d %d %d %17.15lg %3.5g \t [ RANKS N NELEMENTS DOFS ITERATIONS ELAPSEDTIME PRECONMEMORY] \n",
          size, mesh->N, globalElements, globalDofs, Niter, globalElapsed, solver->precon->preconBytes/(1E9));
    }
  }
  return Niter;

}
