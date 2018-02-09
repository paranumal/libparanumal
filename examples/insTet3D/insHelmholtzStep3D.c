#include "ins3D.h"

// complete a time step using LSERK4
void insHelmholtzStep3D(ins_t *ins, iint tstep,  iint haloBytes,
			dfloat * sendBuffer, dfloat * recvBuffer, 
			char   * options){
  
  mesh3D *mesh = ins->mesh; 
  solver_t *solver = ins->vSolver; 
  
  dfloat t = tstep*ins->dt + ins->dt;
  
  iint offset = mesh->Nelements+mesh->totalHaloPairs;
  int subcycling = (strstr(options,"SUBCYCLING")) ? 1:0;
  
  
   // compute all forcing i.e. f^(n+1) - grad(Pr)
  ins->helmholtzRhsForcingKernel(mesh->Nelements,
			                           subcycling,
                                 mesh->o_vgeo,
                                 mesh->o_MM,
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
  
  ins->helmholtzRhsIpdgBCKernel(mesh->Nelements,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                ins->tau,
                                t,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                mesh->o_vgeo,
                                mesh->o_sgeo,
                                mesh->o_EToB,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_DtT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                ins->o_rhsU,
                                ins->o_rhsV,
                                ins->o_rhsW);

  //use intermediate buffer for solve storage TODO: fix this later. Should be able to pull out proper buffer in elliptic solve
  iint Ntotal = offset*mesh->Np;
  ins->o_UH.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
  ins->o_VH.copyFrom(ins->o_V,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
  ins->o_WH.copyFrom(ins->o_W,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));

  //printf("Solving for Ux \n");
  ins->NiterU = ellipticSolveTet3D(solver, ins->lambda, ins->velTOL, ins->o_rhsU, ins->o_UH, ins->vSolverOptions);
  
  //printf("Solving for Uy \n");
  ins->NiterV = ellipticSolveTet3D(solver, ins->lambda, ins->velTOL, ins->o_rhsV, ins->o_VH, ins->vSolverOptions);
  
  //printf("Solving for Uz \n");
  ins->NiterW = ellipticSolveTet3D(solver, ins->lambda, ins->velTOL, ins->o_rhsW, ins->o_WH, ins->vSolverOptions);
  
  //copy into next stage's storage
  int index1 = (ins->index+1)%3; //hard coded for 3 stages
  ins->o_UH.copyTo(ins->o_U,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);
  ins->o_VH.copyTo(ins->o_V,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);  
  ins->o_WH.copyTo(ins->o_W,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);  
}
