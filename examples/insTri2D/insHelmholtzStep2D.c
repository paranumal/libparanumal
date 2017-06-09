#include "ins2D.h"

// complete a time step using LSERK4
void insHelmholtzStep2D(ins_t *ins, iint tstep,  iint haloBytes,
				               dfloat * sendBuffer, dfloat * recvBuffer, 
                       char   * options){

	mesh2D *mesh = ins->mesh; 
  solver_t *solver = ins->vSolver; 
	dfloat t = tstep*ins->dt;

 	

  
  iint stokesSteps = 10;
  dfloat sc = (tstep>=stokesSteps) ? 1: 0; // switch off advection for first steps

  // compute all forcing i.e. f^(n+1) - grad(Pr)
  ins->helmholtzRhsForcingKernel(mesh->Nelements,
                                 mesh->o_vgeo,
                                 mesh->o_MM,
                                 ins->dt, 
                                 ins->a0,
                                 ins->a1,
                                 ins->a2,
                                 sc*ins->b0,
                                 sc*ins->b1,
                                 sc*ins->b2,
                                 ins->o_U,
                                 ins->o_V,
                                 ins->o_NU,
                                 ins->o_NV,
                                 ins->o_Px,
                                 ins->o_Py,
                                 ins->o_rhsU,
                                 ins->o_rhsV);
  
  ins->helmholtzRhsIpdgBCKernel(mesh->Nelements,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                ins->tau,
                                t,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_vgeo,
                                mesh->o_sgeo,
                                mesh->o_EToB,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                ins->o_rhsU,
                                ins->o_rhsV);

  printf("Solving for Ux \n");
  ellipticSolveTri2D( solver, ins->lambda, ins->o_rhsU, ins->o_U, ins->vSolverOptions);

  printf("Solving for Uy \n");
  ellipticSolveTri2D(solver, ins->lambda, ins->o_rhsV, ins->o_V, ins->vSolverOptions);
}
