#include "ins.h"

void insVelocityRhs(ins_t *ins, dfloat time, int stage, occa::memory o_rhsU, occa::memory o_rhsV, occa::memory o_rhsW){
  
  mesh_t *mesh = ins->mesh; 

  if (ins->options.compareArgs("TIME INTEGRATOR", "ARK")) { 
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
                           o_rhsU,
                           o_rhsV,
                           o_rhsW);
  } else if (ins->options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    // rhsU^s = MM*(\sum^s b_i U^n-i - \sum^s-1 a_i N(U^n-i) + \sum^s-1 c_i GP^n-i)/nu dt
    ins->velocityRhsKernel(mesh->Nelements,
                           mesh->o_vgeo,
                           mesh->o_MM,
                           ins->idt,
                           ins->inu,
                           ins->o_extbdfA,
                           ins->o_extbdfB,
                           ins->o_extbdfC,
                           ins->fieldOffset,
                           ins->o_U,
                           ins->o_NU,
                           ins->o_GP,
                           o_rhsU,
                           o_rhsV,
                           o_rhsW);
  }
}
