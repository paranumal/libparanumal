#include "ins.h"

// complete a time step using LSERK4
void insHelmholtzStep(ins_t *ins, dfloat time){
  
  mesh_t *mesh = ins->mesh; 
  solver_t *usolver = ins->uSolver; 
  solver_t *vsolver = ins->vSolver; 
  solver_t *wsolver = ins->wSolver; 
  
  //dfloat t = tstep*ins->dt + ins->dt;

  dlong offset = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
  int subcycling = (strstr(options,"SUBCYCLING")) ? 1:0;

  occaTimerTic(mesh->device,"HelmholtzRhsForcing"); 
  // compute all forcing i.e. f^(n+1) - grad(Pr)
  ins->helmholtzRhsForcingKernel(mesh->Nelements,
                                 subcycling,
                                 mesh->o_vgeo,
                                 mesh->o_MM,
                                 ins->idt,
                                 ins->inu,
                                 ins->a0, ins->a1, ins->a2,
                                 ins->b0, ins->b1, ins->b2,
                                 ins->c0, ins->c1, ins->c2,
                                 offset,
                                 ins->o_U,
                                 ins->o_NU,
                                 ins->o_gradP,
                                 ins->o_rhsU,
                                 ins->o_rhsV,
                                 ins->o_rhsW);
  occaTimerToc(mesh->device,"HelmholtzRhsForcing"); 

  if (strstr(ins->vSolverOptions,"CONTINUOUS")) {
    ins->helmholtzRhsBCKernel(mesh->Nelements,
                              mesh->o_ggeo,
                              mesh->o_sgeo,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              mesh->o_vmapM,
                              mesh->o_sMT,
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
    ellipticParallelGatherScatter(mesh, mesh->ogs, ins->o_rhsU, dfloatString, "add");  
    ellipticParallelGatherScatter(mesh, mesh->ogs, ins->o_rhsV, dfloatString, "add");  
    if (ins->dim==3)
      ellipticParallelGatherScatter(mesh, mesh->ogs, ins->o_rhsW, dfloatString, "add");  
    if (usolver->Nmasked) mesh->maskKernel(usolver->Nmasked, usolver->o_maskIds, ins->o_rhsU);
    if (vsolver->Nmasked) mesh->maskKernel(vsolver->Nmasked, vsolver->o_maskIds, ins->o_rhsV);
    if (ins->dim==3)
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
                                  mesh->o_Dmatrices,
                                  mesh->o_LIFTT,
                                  mesh->o_MM,
                                  ins->o_rhsU,
                                  ins->o_rhsV,
                                  ins->o_rhsW);
    occaTimerToc(mesh->device,"HelmholtzRhsIpdg");   
  }

  //copy current velocity fields as initial guess
  dlong Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
  ins->o_UH.copyFrom(ins->o_U,Ntotal*sizeof(dfloat));
  ins->o_VH.copyFrom(ins->o_V,Ntotal*sizeof(dfloat));
  if (ins->dim==3)
    ins->o_WH.copyFrom(ins->o_W,Ntotal*sizeof(dfloat));

  if (strstr (ins->vSolverOptions,"CONTINUOUS")) {
    if (usolver->Nmasked) mesh->maskKernel(usolver->Nmasked, usolver->o_maskIds, ins->o_UH);
    if (vsolver->Nmasked) mesh->maskKernel(vsolver->Nmasked, vsolver->o_maskIds, ins->o_VH);
    if (ins->dim==3)
      if (wsolver->Nmasked) mesh->maskKernel(wsolver->Nmasked, wsolver->o_maskIds, ins->o_WH);

  }
  
  occaTimerTic(mesh->device,"Ux-Solve");
  ins->NiterU = ellipticSolve(usolver, ins->lambda, ins->velTOL, ins->o_rhsU, ins->o_UH);
  occaTimerToc(mesh->device,"Ux-Solve"); 

  occaTimerTic(mesh->device,"Uy-Solve");
  ins->NiterV = ellipticSolve(vsolver, ins->lambda, ins->velTOL, ins->o_rhsV, ins->o_VH);
  occaTimerToc(mesh->device,"Uy-Solve");

  if (ins->dim==3) {
    occaTimerTic(mesh->device,"Uz-Solve");
    ins->NiterW = ellipticSolve(wsolver, ins->lambda, ins->velTOL, ins->o_rhsW, ins->o_WH);
    occaTimerToc(mesh->device,"Uz-Solve");
  }

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

  //copy into intermediate stage storage
  ins->o_UH.copyTo(ins->o_Uhat,Ntotal*sizeof(dfloat),0*Ntotal*sizeof(dfloat),0);
  ins->o_VH.copyTo(ins->o_Uhat,Ntotal*sizeof(dfloat),1*Ntotal*sizeof(dfloat),0);    
  if (ins->dim==3)
    ins->o_WH.copyTo(ins->o_Uhat,Ntotal*sizeof(dfloat),2*Ntotal*sizeof(dfloat),0);    
}
