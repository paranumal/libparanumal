#include "insHex3D.h"

void insHelmholtzStepHex3D(ins_t *ins, int tstep, char   * options){
  
  mesh3D *mesh = ins->mesh; 
  solver_t *usolver = ins->uSolver; 
  solver_t *vsolver = ins->vSolver; 
  solver_t *wsolver = ins->wSolver; 
  
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
                                 ins->o_W,
                                 ins->o_NU,
                                 ins->o_NV,
                                 ins->o_NW,
                                 ins->o_Px,
                                 ins->o_Py,
                                 ins->o_Pz,
                                 ins->o_rhsU,
                                 ins->o_rhsV,
                                 ins->o_rhsW);
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
                              mesh->o_z,
                              ins->o_VmapB,
                              ins->o_rhsU,
                              ins->o_rhsV,
                              ins->o_rhsW);

    // gather-scatter
    ellipticParallelGatherScatterHex3D(mesh, mesh->ogs, ins->o_rhsU, dfloatString, "add");  
    ellipticParallelGatherScatterHex3D(mesh, mesh->ogs, ins->o_rhsV, dfloatString, "add");  
    ellipticParallelGatherScatterHex3D(mesh, mesh->ogs, ins->o_rhsW, dfloatString, "add");  
    if (usolver->Nmasked) mesh->maskKernel(usolver->Nmasked, usolver->o_maskIds, ins->o_rhsU);
    if (vsolver->Nmasked) mesh->maskKernel(vsolver->Nmasked, vsolver->o_maskIds, ins->o_rhsV);
    if (wsolver->Nmasked) mesh->maskKernel(wsolver->Nmasked, wsolver->o_maskIds, ins->o_rhsW);

  } else if (strstr(ins->vSolverOptions,"IPDG")) {
    occaTimerTic(mesh->device,"HelmholtzRhsIpdg");   
    ins->helmholtzRhsIpdgBCKernel(mesh->Nelements,
                                  mesh->o_vmapM,
                                  ins->tau,
                                  t,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  mesh->o_EToB,
                                  mesh->o_D,
                                  ins->o_rhsU,
                                  ins->o_rhsV,
                                  ins->o_rhsW);
    occaTimerToc(mesh->device,"HelmholtzRhsIpdg");   
  }

  //use intermediate buffer for solve storage TODO: fix this later. Should be able to pull out proper buffer in elliptic solve
  dlong Ntotal = offset*mesh->Np;
  ins->o_UH.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
  ins->o_VH.copyFrom(ins->o_V,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
  ins->o_WH.copyFrom(ins->o_W,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));

  if (strstr (ins->vSolverOptions,"CONTINUOUS")) {
    if (usolver->Nmasked) mesh->maskKernel(usolver->Nmasked, usolver->o_maskIds, ins->o_UH);
    if (vsolver->Nmasked) mesh->maskKernel(vsolver->Nmasked, vsolver->o_maskIds, ins->o_VH);
    if (wsolver->Nmasked) mesh->maskKernel(wsolver->Nmasked, wsolver->o_maskIds, ins->o_WH);
  }

  occaTimerTic(mesh->device,"Ux-Solve");
  ins->NiterU = ellipticSolveHex3D(usolver, ins->lambda, ins->velTOL, ins->o_rhsU, ins->o_UH, ins->vSolverOptions);
  occaTimerToc(mesh->device,"Ux-Solve"); 

  occaTimerTic(mesh->device,"Uy-Solve");
  ins->NiterV = ellipticSolveHex3D(vsolver, ins->lambda, ins->velTOL, ins->o_rhsV, ins->o_VH, ins->vSolverOptions);
  occaTimerToc(mesh->device,"Uy-Solve");

  occaTimerTic(mesh->device,"Uz-Solve");
  ins->NiterW = ellipticSolveHex3D(wsolver, ins->lambda, ins->velTOL, ins->o_rhsW, ins->o_WH, ins->vSolverOptions);
  occaTimerToc(mesh->device,"Uz-Solve");

  if (strstr(ins->vSolverOptions,"CONTINUOUS")) {
    ins->helmholtzAddBCKernel(mesh->Nelements,
                            t,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vmapM,
                            ins->o_VmapB,
                            ins->o_UH,
                            ins->o_VH,
                            ins->o_WH);
  }

  //copy into next stage's storage
  int index1 = (ins->index+1)%3; //hard coded for 3 stages
  ins->o_UH.copyTo(ins->o_U,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);
  ins->o_VH.copyTo(ins->o_V,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);    
  ins->o_WH.copyTo(ins->o_W,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);    
}
