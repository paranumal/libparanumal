#include "ins2D.h"

// complete a time step using LSERK4
void insHelmholtzStepSS2D(ins_t *ins, iint tstep,  iint haloBytes,
			dfloat * sendBuffer, dfloat * recvBuffer, 
			char   * options){
  
  mesh2D *mesh = ins->mesh; 
  solver_t *solver = ins->vSolver; 
 // dfloat t = tstep*ins->dt;
  dfloat t = tstep*ins->dt + ins->dt;
  
  iint offset = mesh->Nelements+mesh->totalHaloPairs;
  iint rhsPackingMode = (strstr(options, "VECTORHELMHOLTZ")) ? 1:0;
    
  if(strstr(options,"SUBCYCLING")){
     // compute all forcing i.e. f^(n+1) - grad(Pr)
    ins->helmholtzRhsForcingKernel(mesh->Nelements,
                                   ins->dt,
                                   ins->g0,
                                   mesh->o_vgeo,
                                   mesh->o_MM,
                                   ins->o_Ut,
                                   ins->o_Vt,
                                   ins->o_rhsU,
                                   ins->o_rhsV);
  }
  else{
    ins->helmholtzRhsForcingKernel(mesh->Nelements,
                                   ins->dt,
                                   ins->g0,
                                   mesh->o_vgeo,
                                   mesh->o_MM,
                                   ins->o_Ut,
                                   ins->o_Vt,
                                   ins->o_rhsU,
                                   ins->o_rhsV);
  }
  
  ins->helmholtzRhsIpdgBCKernel(mesh->Nelements,
				                        rhsPackingMode,
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

  
    iint Niter;
    printf("Solving for Ux: Niter= ");
    Niter = ellipticSolveTri2D( solver, ins->lambda, ins->o_rhsU, ins->o_Ut, ins->vSolverOptions);
    printf("%d \n",Niter);

    printf("Solving for Uy: Niter= ");
    Niter = ellipticSolveTri2D(solver, ins->lambda, ins->o_rhsV, ins->o_Vt, ins->vSolverOptions);
    printf("%d \n",Niter);
    //copy into next stage's storage
    ins->index = (ins->index+1)%3; //hard coded for 3 stages
    const iint Ntotal = (mesh->Nelements + mesh->totalHaloPairs)*mesh->Np; 
    ins->o_Ut.copyTo(ins->o_U,Ntotal*sizeof(dfloat),ins->index*Ntotal*sizeof(dfloat),0);
    ins->o_Vt.copyTo(ins->o_V,Ntotal*sizeof(dfloat),ins->index*Ntotal*sizeof(dfloat),0);  

}
