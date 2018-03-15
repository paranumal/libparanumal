#include "ins2D.h"

// complete a time step using LSERK4
void insPoissonStepSS2D(ins_t *ins, int tstep, int haloBytes,
				               dfloat * sendBuffer, dfloat * recvBuffer,
				                char   * options){

  mesh2D *mesh = ins->mesh;
  solver_t *solver = ins->pSolver;

  // dfloat t0 = tstep*ins->dt;
  dfloat t = tstep*ins->dt + ins->dt;

  //hard coded for 3 stages.
  int index0 = (ins->index+0)%3; // Still in current time
  int Ntotal = (mesh->Nelements+mesh->totalHaloPairs);
  int offset = index0*Ntotal;


// Compute Curl(Curl(U)) and store in rhsU 
  ins->poissonRhsCurlKernel(mesh->Nelements,
                      offset,
                      mesh->o_vgeo,
                      mesh->o_DrT,
                      mesh->o_DsT,
                      ins->o_U,
                      ins->o_V,
                      ins->o_rhsU,
                      ins->o_rhsV);


  // Compute derived Neumann data : result is lifted and multiplied with invJ
  ins->poissonRhsNeumannKernel(mesh->Nelements,
                                ins->index,
                                Ntotal,
                                ins->dt,
                                ins->a0,
                                ins->a1,
                                ins->a2,
                                mesh->o_sgeo,
                                mesh->o_LIFTT,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                t,
                                mesh->o_x,
                                mesh->o_y,
                                ins->o_NU,
                                ins->o_NV,
                                ins->o_rhsU, // holds curl of Vorticity X
                                ins->o_rhsV, // holds curl of Vorticity Y
                                ins->o_WN,   //
                                ins->o_rhsP); 


if(strstr(options,"BROKEN")){
  const int voffset = 0; 
  // compute div(ut) and store on rhsU
  ins->divergenceVolumeKernel(mesh->Nelements,
                             mesh->o_vgeo,
                             mesh->o_DrT,
                             mesh->o_DsT,
                             voffset,
                             ins->o_Ut,
                             ins->o_Vt,
                             ins->o_rhsU);
}






  
 // compute all forcing i.e. f^(n+1) - grad(Pr)
  ins->poissonRhsForcingKernel(mesh->Nelements,
                              ins->dt,  
                              ins->g0,
                              mesh->o_vgeo,
                              mesh->o_MM,
                              ins->o_rhsU, // div(Ut)
                              ins->o_rhsP);
 //
  const int pressure_solve = 1; // Solve for Pressure
  ins->poissonRhsIpdgBCKernel(mesh->Nelements,
                                pressure_solve,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                ins->tau,
                                t,
                                ins->dt,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_vgeo,
                                mesh->o_sgeo,
                                mesh->o_EToB,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                ins->o_rhsP);



  ins->o_Pt.copyFrom(ins->o_P,Ntotal*mesh->Np*sizeof(dfloat),0,ins->index*Ntotal*mesh->Np*sizeof(dfloat));
  int Niter;
  printf("Solving for Pr: Niter= ");
  Niter= ellipticSolveTri2D(solver, 0.0, ins->presTOL, ins->o_rhsP, ins->o_Pt,  ins->pSolverOptions); 
  printf("%d\n", Niter); 

  int index1 = (ins->index+1)%3; //hard coded for 3 stages
  ins->o_Pt.copyTo(ins->o_P,Ntotal*mesh->Np*sizeof(dfloat),index1*Ntotal*mesh->Np*sizeof(dfloat),0);


if(strstr(options,"BROKEN")){
  const int poffset = index1*Ntotal*mesh->Np; 
  // Compute Volume Contribution of gradient of pressure gradient
  ins->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_DrT,
                            mesh->o_DsT,
                            poffset,  
                            ins->o_P,
                            ins->o_rhsU,
                            ins->o_rhsV);
}
  

ins->poissonUpdateKernel(mesh->Nelements,
                            ins->dt,
                            ins->g0,
                            ins->o_rhsU,
                            ins->o_rhsV,
                            ins->o_Ut,
                            ins->o_Vt);

}
