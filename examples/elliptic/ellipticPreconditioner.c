#include "ellipticTri2D.h"


void ellipticPreconditioner2D(solver_t *solver,
            dfloat lambda,
            occa::memory &o_r,
            occa::memory &o_z,
            const char *options){

  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  
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
      //if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_r);

      if(ogs->NhaloGather) {
        precon->partialblockJacobiKernel(solver->NglobalGatherElements, 
                                solver->o_globalGatherElementList,
                                invLambda, mesh->o_vgeo, precon->o_invMM, solver->o_rtmp, o_z);
        mesh->device.finish();
        mesh->device.setStream(solver->dataStream);
        mesh->gatherKernel(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherLocalIds, o_z, ogs->o_haloGatherTmp);
        ogs->o_haloGatherTmp.asyncCopyTo(ogs->haloGatherTmp);
        mesh->device.setStream(solver->defaultStream);
      }
      if(solver->NlocalGatherElements){
        precon->partialblockJacobiKernel(solver->NlocalGatherElements, 
                                solver->o_localGatherElementList,
                                invLambda, mesh->o_vgeo, precon->o_invMM, solver->o_rtmp, o_z);
      }

      // finalize gather using local and global contributions
      if(ogs->NnonHaloGather) mesh->gatherScatterKernel(ogs->NnonHaloGather, ogs->o_nonHaloGatherOffsets, ogs->o_nonHaloGatherLocalIds, o_z);

      // C0 halo gather-scatter (on data stream)
      if(ogs->NhaloGather) {
        mesh->device.setStream(solver->dataStream);
        mesh->device.finish();

        // MPI based gather scatter using libgs
        gsParallelGatherScatter(ogs->haloGsh, ogs->haloGatherTmp, dfloatString, "add");

        // copy totally gather halo data back from HOST to DEVICE
        ogs->o_haloGatherTmp.asyncCopyFrom(ogs->haloGatherTmp);

        // do scatter back to local nodes
        mesh->scatterKernel(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherLocalIds, ogs->o_haloGatherTmp, o_z);
        mesh->device.setStream(solver->defaultStream);
      }

      mesh->device.finish();    
      mesh->device.setStream(solver->dataStream);
      mesh->device.finish();    
      mesh->device.setStream(solver->defaultStream);

      solver->dotMultiplyKernel(mesh->Nelements*mesh->Np, ogs->o_invDegree, o_z, o_z);

      //post-mask
      if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_z);
    }

  } else if (strstr(options, "SEMFEM")) {

    o_z.copyFrom(o_r);
    solver->dotMultiplyKernel(mesh->Nelements*mesh->Np, solver->o_invDegree, o_z, o_z);
    precon->SEMFEMInterpKernel(mesh->Nelements,mesh->o_SEMFEMAnterp,o_z,precon->o_rFEM);
    meshParallelGather(mesh, precon->FEMogs, precon->o_rFEM, precon->o_GrFEM);
    occaTimerTic(mesh->device,"parALMOND");
    parAlmondPrecon(precon->parAlmond, precon->o_GzFEM, precon->o_GrFEM);
    occaTimerToc(mesh->device,"parALMOND");
    meshParallelScatter(mesh, precon->FEMogs, precon->o_GzFEM, precon->o_zFEM);
    precon->SEMFEMAnterpKernel(mesh->Nelements,mesh->o_SEMFEMAnterp,precon->o_zFEM,o_z);
    solver->dotMultiplyKernel(mesh->Nelements*mesh->Np, solver->o_invDegree, o_z, o_z);

    ellipticParallelGatherScatterTri2D(mesh, mesh->ogs, o_z, dfloatString, "add");

  } else if(strstr(options, "JACOBI")){

    dlong Ntotal = mesh->Np*mesh->Nelements;
    // Jacobi preconditioner
    occaTimerTic(mesh->device,"dotDivideKernel");
    solver->dotMultiplyKernel(Ntotal, o_r, precon->o_invDiagA, o_z);
    occaTimerToc(mesh->device,"dotDivideKernel");
  
  } else{ // turn off preconditioner
    o_z.copyFrom(o_r);
  }
}
