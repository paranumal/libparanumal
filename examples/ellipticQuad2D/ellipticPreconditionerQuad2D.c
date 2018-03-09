#include "ellipticQuad2D.h"

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

  } else if(strstr(options, "OAS")){
    
    //patch solve
    //ellipticPatchSmootherQuad2D(solver,o_r,o_z,options);
    smoothQuad2D(precon->OASsmoothArgs, o_r, o_z,true);

    occaTimerTic(mesh->device,"coarseGrid");

    // Z1*Z1'*PL1*(Z1*z1) = (Z1*rL)  HMMM
    occaTimerTic(mesh->device,"coarsenKernel");
    precon->coarsenKernel(mesh->Nelements, precon->o_V1, o_r, precon->o_r1);
    occaTimerToc(mesh->device,"coarsenKernel");

    occaTimerTic(mesh->device,"ALMOND");
    parAlmondPrecon(precon->parAlmond, precon->o_z1, precon->o_r1);
    occaTimerToc(mesh->device,"ALMOND");

    // prolongate from P1 to PN, adding patch and coarse solves together
    occaTimerTic(mesh->device,"prolongateKernel");
    precon->prolongateKernel(mesh->Nelements, precon->o_V1, precon->o_z1, solver->o_z);
    occaTimerToc(mesh->device,"prolongateKernel");

    occaTimerToc(mesh->device,"coarseGrid");

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

    ellipticParallelGatherScatterQuad2D(mesh, mesh->ogs, o_z, dfloatString, "add");

  } else if(strstr(options, "JACOBI")){

    dlong Ntotal = mesh->Np*mesh->Nelements;
    // Jacobi preconditioner
    occaTimerTic(mesh->device,"dotDivideKernel");
    solver->dotMultiplyKernel(Ntotal, o_r, precon->o_invDiagA, o_z);
    occaTimerToc(mesh->device,"dotDivideKernel");
  
  }  else{ // turn off preconditioner
    o_z.copyFrom(o_r); 
  }
}
