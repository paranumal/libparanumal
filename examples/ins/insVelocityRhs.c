#include "ins.h"

void insVelocityRhs(ins_t *ins, dfloat time, int stage){
  
  mesh_t *mesh = ins->mesh; 

  // compute all forcing i.e. 
  // rhsU^s = MM*(U^n - \sum^s-1 ea_si N(U^i) + \sum^s-1 ia_si LU^i - \sum^s-1 pa_si GP^i)/ia_ss nu dt
  ins->velocityRhsKernel(mesh->Nelements,
                         stage,
                         mesh->o_vgeo,
                         mesh->o_MM,
                         ins->idt,
                         ins->inu,
                         ins->o_erkA,
                         ins->o_irkA,
                         ins->o_prkA,
                         ins->fieldOffset,
                         ins->o_U,
                         ins->o_NU,
                         ins->o_LU,
                         ins->o_GP,
                         ins->o_rhsU,
                         ins->o_rhsV,
                         ins->o_rhsW);
}
