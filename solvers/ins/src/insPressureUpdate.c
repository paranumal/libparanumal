#include "ins.h"

void insPressureUpdate(ins_t *ins, dfloat time, int stage, occa::memory o_rkP){

  mesh_t *mesh = ins->mesh;
  
  // P^s = PI + \sum^s-1 prk_si P^i 
  occaTimerTic(mesh->device,"PressureUpdate");
  ins->pressureUpdateKernel(mesh->Nelements,
                            stage,
		                        ins->o_prkB,
                            ins->fieldOffset,
                            ins->o_PI,
                            ins->o_P,
                            o_rkP);
  occaTimerToc(mesh->device,"PressureUpdate");
}
