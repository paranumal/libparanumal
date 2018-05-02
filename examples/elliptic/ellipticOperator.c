#include "elliptic.h"

void ellipticOperator(solver_t *solver, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *options){

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
          mesh->o_ggeo, mesh->o_SrrT, mesh->o_SrsT, mesh->o_SsrT, mesh->o_SssT,
          mesh->o_MM, lambda, o_q, o_Aq);
      mesh->device.finish();
      mesh->device.setStream(solver->dataStream);
      mesh->gatherKernel(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherLocalIds, o_Aq, ogs->o_haloGatherTmp);
      ogs->o_haloGatherTmp.asyncCopyTo(ogs->haloGatherTmp);
      mesh->device.setStream(solver->defaultStream);
    }
    if(solver->NlocalGatherElements){
        solver->partialAxKernel(solver->NlocalGatherElements, solver->o_localGatherElementList,
            mesh->o_ggeo, mesh->o_SrrT, mesh->o_SrsT, mesh->o_SsrT, mesh->o_SssT,
            mesh->o_MM, lambda, o_q, o_Aq);
    }
    if(solver->allNeumann) {
      o_tmp.copyTo(tmp);

      for(dlong n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
      alphaG *= solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
    }

    // finalize gather using local and global contributions
    if(ogs->NnonHaloGather) mesh->gatherScatterKernel(ogs->NnonHaloGather, ogs->o_nonHaloGatherOffsets, ogs->o_nonHaloGatherLocalIds, o_Aq);

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
    dlong offset = 0;
    dfloat alpha = 0., alphaG =0.;
    dlong Nblock = solver->Nblock;
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

      for(dlong n=0;n<Nblock;++n)
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
  } 

  occaTimerToc(mesh->device,"AxKernel");
}
