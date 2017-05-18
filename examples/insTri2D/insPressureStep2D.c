#include "ins2D.h"

// complete a time step using LSERK4
void insPressureStep2D(solver_t *ins, iint tstep, iint haloBytes,
				       dfloat * sendBuffer, dfloat * recvBuffer, 
				        char   * options){

mesh2D *mesh = ins->mesh; 

dfloat t = tstep*ins->dt;

  // Compute Pressure RHS
  
 // ins->pressureRhsVolumeKernel(	mesh->Nelements,
	// 							ins->g0,
	// 							ins->dt,
	// 							mesh->o_vgeo,
	// 							mesh->o_DrT,
	// 							mesh->o_DsT,
	// 							ins->o_UI,
	// 							ins->o_rhsPr);


  // Solve for Pressure









   
}
