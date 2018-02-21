#include "ellipticTet3D.h"

void ellipticStartHaloExchange3D(mesh3D *mesh, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer);
void ellipticEndHaloExchange3D(mesh3D *mesh, occa::memory &o_q, dfloat *recvBuffer);
void ellipticParallelGatherScatterTet3D(mesh3D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op);
dfloat ellipticScaledAdd(solver_t *solver, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b);
void ellipticOperator3D(solver_t *solver, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *options);


void ellipticPreconditioner3D(solver_t *solver,
            dfloat lambda,
            occa::memory &o_r,
            occa::memory &o_z,
            const char *options){

  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  ogs_t    *ogs = solver->ogs; // C0 Gather ScatterTri info

  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;

  if (strstr(options, "FULLALMOND")||strstr(options, "MULTIGRID")) {

    occaTimerTic(mesh->device,"parALMOND");
    parAlmondPrecon(precon->parAlmond, o_z, o_r);
    occaTimerToc(mesh->device,"parALMOND");


  } else if(strstr(options, "MASSMATRIX")){

    dfloat invLambda = 1./lambda;

    if (strstr(options,"IPDG")) {
      occaTimerTic(mesh->device,"blockJacobiKernel");
      precon->blockJacobiKernel(mesh->Nelements, invLambda, mesh->o_vgeo, precon->o_invMM, o_r, o_z);
      occaTimerToc(mesh->device,"blockJacobiKernel");
    } else if (strstr(options,"CONTINUOUS")) {
      ogs_t *ogs = solver->mesh->ogs;

      solver->dotMultiplyKernel(mesh->Nelements*mesh->Np, ogs->o_invDegree, o_r, solver->o_rtmp);
      //pre-mask
      if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_r);

      if(ogs->NhaloGather) {
        precon->partialblockJacobiKernel(solver->NglobalGatherElements, 
                                solver->o_globalGatherElementList,
                                invLambda, mesh->o_vgeo, precon->o_invMM, solver->o_rtmp, o_z);
        mesh->gatherKernel(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherLocalIds, o_z, ogs->o_haloGatherTmp);
        ogs->o_haloGatherTmp.copyTo(ogs->haloGatherTmp);
      }
      if(solver->NlocalGatherElements){
        precon->partialblockJacobiKernel(solver->NlocalGatherElements, 
                                solver->o_localGatherElementList,
                                invLambda, mesh->o_vgeo, precon->o_invMM, solver->o_rtmp, o_z);
      }

      // C0 halo gather-scatter (on data stream)
      if(ogs->NhaloGather) {
        occa::streamTag tag;

        // MPI based gather scatter using libgs
        gsParallelGatherScatter(ogs->haloGsh, ogs->haloGatherTmp, dfloatString, "add");

        // copy totally gather halo data back from HOST to DEVICE
        mesh->device.setStream(solver->dataStream);
        ogs->o_haloGatherTmp.asyncCopyFrom(ogs->haloGatherTmp);

        // wait for async copy
        tag = mesh->device.tagStream();
        mesh->device.waitFor(tag);

        // do scatter back to local nodes
        mesh->scatterKernel(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherLocalIds, ogs->o_haloGatherTmp, o_z);

        // make sure the scatter has finished on the data stream
        tag = mesh->device.tagStream();
        mesh->device.waitFor(tag);
      }

      // finalize gather using local and global contributions
      mesh->device.setStream(solver->defaultStream);
      if(ogs->NnonHaloGather) mesh->gatherScatterKernel(ogs->NnonHaloGather, ogs->o_nonHaloGatherOffsets, ogs->o_nonHaloGatherLocalIds, o_z);

      solver->dotMultiplyKernel(mesh->Nelements*mesh->Np, ogs->o_invDegree, o_z, o_z);

      //post-mask
      if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_z);
    }

  } else if(strstr(options, "JACOBI")){

    iint Ntotal = mesh->Np*mesh->Nelements;
    // Jacobi preconditioner
    occaTimerTic(mesh->device,"dotDivideKernel");
    solver->dotMultiplyKernel(Ntotal, o_r, precon->o_invDiagA, o_z);
    occaTimerToc(mesh->device,"dotDivideKernel");
  }
  else{ // turn off preconditioner
    o_z.copyFrom(o_r);
  }
}
