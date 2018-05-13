#include "ins.h"

void insPressureRhs(ins_t *ins, dfloat time, int stage){

  mesh_t *mesh = ins->mesh;

  // rhsP = Div Uhat
  insDivergence(ins, time, ins->o_rkU, ins->o_rhsP);
  
  // rhsP = -MM*Div Uhat/pa_ss dt
  //dfloat g0 = 1.0/ins->prkA[stage->ins->Nstages+stage];
  occaTimerTic(mesh->device,"PoissonRhsForcing");
  ins->pressureRhsKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_MM,
                              ins->idt,  
                              ins->o_rhsP);
  occaTimerToc(mesh->device,"PoissonRhsForcing");

#if 0
  //add penalty from jumps in previous pressure
  ins->poissonPenaltyKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_vgeo,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                ins->tau,
                                mesh->o_x,
                                mesh->o_y,
                                t,
                                ins->dt,
                                ins->c0,
                                ins->c1,
                                ins->c2,
                                ins->index,
                                (mesh->Nelements+mesh->totalHaloPairs),
                                ins->o_P,
                                ins->o_rhsP);
  #endif
}
