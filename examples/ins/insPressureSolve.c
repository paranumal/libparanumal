#include "ins.h"

// solve Ap = rhsP
void insPressureSolve(ins_t *ins, dfloat time){

  mesh_t *mesh = ins->mesh;
  elliptic_t *solver = ins->pSolver;
  
  if (ins->pOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    ins->poissonRhsBCKernel(mesh->Nelements,
                            mesh->o_ggeo,
                            mesh->o_sgeo,
                            mesh->o_Dmatrices,
                            mesh->o_Smatrices,
                            mesh->o_MM,
                            mesh->o_vmapM,
                            mesh->o_sMT,
                            time,
                            ins->dt,
                            ins->c0, ins->c1, ins->c2,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            ins->o_PmapB,
                            ins->o_rhsP);
  } else if (ins->pOptions.compareArgs("DISCRETIZATION","IPDG")) {
    occaTimerTic(mesh->device,"PoissonRhsIpdg"); 
    ins->poissonRhsIpdgBCKernel(mesh->Nelements,
                                  mesh->o_vmapM,
                                  solver->tau,
                                  time,
                                  ins->dt,
                                  ins->c0, ins->c1, ins->c2,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  mesh->o_EToB,
                                  mesh->o_Dmatrices,
                                  mesh->o_LIFTT,
                                  mesh->o_MM,
                                  ins->o_rhsP);
    occaTimerToc(mesh->device,"PoissonRhsIpdg");
  }

  //keep current PI as the initial guess?

  // gather-scatter
  if(ins->pOptions.compareArgs("DISCRETIZATION","CONTINUOUS")){
    ellipticParallelGatherScatter(mesh, mesh->ogs, ins->o_rhsP, dfloatString, "add");  
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_rhsP);
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_PI);
  }

  occaTimerTic(mesh->device,"Pr Solve");
  ins->NiterP = ellipticSolve(solver, 0.0, ins->presTOL, ins->o_rhsP, ins->o_PI); 
  occaTimerToc(mesh->device,"Pr Solve"); 

  if (ins->pOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    ins->poissonAddBCKernel(mesh->Nelements,
                            time,
                            ins->dt,
                            ins->c0, ins->c1, ins->c2,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vmapM,
                            ins->o_PmapB,
                            ins->o_PI);
  }
}
