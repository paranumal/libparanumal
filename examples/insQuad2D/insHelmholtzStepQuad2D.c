#include "insQuad2D.h"

void insHelmholtzStepQuad2D(ins_t *ins, int tstep, char   * options){
  
  mesh2D *mesh = ins->mesh; 
  solver_t *solver = ins->vSolver; 
  
  dfloat t = tstep*ins->dt + ins->dt;

  dlong offset = mesh->Nelements+mesh->totalHaloPairs;
  int subcycling = (strstr(options,"SUBCYCLING")) ? 1:0; //TODO: Move this to kernel #define?

  occaTimerTic(mesh->device,"HelmholtzRhsForcing"); 
  // compute all forcing i.e. f^(n+1) - grad(Pr)
  ins->helmholtzRhsForcingKernel(mesh->Nelements,
                                 subcycling,
                                 mesh->o_vgeo,
                                 ins->idt,
                                 ins->inu,
                                 ins->a0,
                                 ins->a1,
                                 ins->a2,
                                 ins->b0,
                                 ins->b1,
                                 ins->b2,
                                 ins->c0,
                                 ins->c1,
                                 ins->c2,
                                 ins->index,
                                 offset,
                                 ins->o_U,
                                 ins->o_V,
                                 ins->o_NU,
                                 ins->o_NV,
                                 ins->o_Px,
                                 ins->o_Py,
                                 ins->o_rhsU,
                                 ins->o_rhsV);
  occaTimerToc(mesh->device,"HelmholtzRhsForcing"); 

  if (strstr(ins->vSolverOptions,"CONTINUOUS")) {
    ins->helmholtzRhsBCKernel(mesh->Nelements,
                              mesh->o_ggeo,
                              mesh->o_sgeo,
                              mesh->o_D,
                              mesh->o_vmapM,
                              ins->lambda,
                              t,
                              mesh->o_x,
                              mesh->o_y,
                              ins->o_VmapB,
                              ins->o_rhsU,
                              ins->o_rhsV);

    // gather-scatter
    ellipticParallelGatherScatterQuad2D(mesh, mesh->ogs, ins->o_rhsU, dfloatString, "add");  
    ellipticParallelGatherScatterQuad2D(mesh, mesh->ogs, ins->o_rhsV, dfloatString, "add");  
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_rhsU);
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_rhsV);

  } else if (strstr(ins->vSolverOptions,"IPDG")) {
    occaTimerTic(mesh->device,"HelmholtzRhsIpdg");   
    ins->helmholtzRhsIpdgBCKernel(mesh->Nelements,
                                  mesh->o_vmapM,
                                  ins->tau,
                                  t,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  mesh->o_EToB,
                                  mesh->o_D,
                                  ins->o_rhsU,
                                  ins->o_rhsV);
    occaTimerToc(mesh->device,"HelmholtzRhsIpdg");   
  }

  //use intermediate buffer for solve storage TODO: fix this later. Should be able to pull out proper buffer in elliptic solve
  int Ntotal = offset*mesh->Np;
  ins->o_UH.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
  ins->o_VH.copyFrom(ins->o_V,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));


  occaTimerToc(mesh->device,"Ux-Solve");
  ins->NiterU = ellipticSolveQuad2D( solver, ins->lambda, ins->velTOL, ins->o_rhsU, ins->o_UH, ins->vSolverOptions);
  occaTimerToc(mesh->device,"Ux-Solve"); 

  occaTimerToc(mesh->device,"Uy-Solve");
  ins->NiterV = ellipticSolveQuad2D(solver, ins->lambda, ins->velTOL, ins->o_rhsV, ins->o_VH, ins->vSolverOptions);
  occaTimerToc(mesh->device,"Uy-Solve");

  if (strstr(ins->vSolverOptions,"CONTINUOUS")) {
    ins->helmholtzAddBCKernel(mesh->Nelements,
                            t,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_vmapM,
                            ins->o_VmapB,
                            ins->o_UH,
                            ins->o_VH);
  }

  //copy into next stage's storage
  int index1 = (ins->index+1)%3; //hard coded for 3 stages
  ins->o_UH.copyTo(ins->o_U,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);
  ins->o_VH.copyTo(ins->o_V,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);    
}
