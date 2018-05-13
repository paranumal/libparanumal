#include "ins.h"

void insVelocityUpdate(ins_t *ins, dfloat time, int stage, 
                        occa::memory o_rkGP,
                        occa::memory o_rkU){

  mesh_t *mesh = ins->mesh;

  // U^s = Uhat - dt * GP^s + dt*\sum^s-1 pa_si GP^i
  occaTimerTic(mesh->device,"VelocityUpdate");
  ins->velocityUpdateKernel(mesh->Nelements,
                              stage,
                              ins->dt,
                              ins->fieldOffset,
                              ins->o_prkA,
                              o_rkGP,
                              ins->o_GP,
                              o_rkU);
  occaTimerToc(mesh->device,"VelocityUpdate");
}
