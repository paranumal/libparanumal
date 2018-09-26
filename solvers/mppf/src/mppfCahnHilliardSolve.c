#include "mppf.h"

// solve lambda*U + A*U = rhsU
void mppfCahnHilliardSolve(mppf_t *mppf, dfloat time){
  
  mesh_t *mesh = mppf->mesh; 
  elliptic_t *phiSolver = mppf->phiSolver; 
  elliptic_t *psiSolver = mppf->psiSolver; 
  
  if (mppf->phiOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
     mppf->phaseFieldRhsBCKernel(mesh->Nelements,
                              mesh->o_ggeo,
                              mesh->o_sgeo,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              mesh->o_vmapM,
                              mesh->o_EToB,
                              mesh->o_sMT,
                              mppf->lambdaPsi,
                              time,
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              mppf->o_VmapB,
                              mppf->o_rhsPhi);
    
    // gather-scatter
    // ellipticParallelGatherScatter(mesh, mesh->ogs, mppf->o_rhsPhi, dfloatString, "add");
    ogsGatherScatter(mppf->o_rhsPhi, ogsDfloat, ogsAdd, mesh->ogs);

    if (psiSolver->Nmasked) mesh->maskKernel(psiSolver->Nmasked, psiSolver->o_maskIds, mppf->o_rhsPhi);
  } else if (mppf->phiOptions.compareArgs("DISCRETIZATION","IPDG")) {
    // Currently we do not need that deuto homegenous bcs need to be used for more complex bcs
    occaTimerTic(mesh->device,"CahnHilliardRhsIpdgBC");    
    mppf->phaseFieldRhsIpdgBCKernel(mesh->Nelements,
                                  mesh->o_vmapM,
                                  psiSolver->tau,
                                  time,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  mesh->o_EToB,
                                  mesh->o_Dmatrices,
                                  mesh->o_LIFTT,
                                  mesh->o_MM,
                                  mppf->o_rhsPhi);
    occaTimerToc(mesh->device,"CahnHilliardRhsIpdgBC");   
  }

  //copy current field as initial guess? (could use Uhat or beter guess)
  dlong Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
  mppf->o_rkPhi.copyFrom(mppf->o_Phi,Ntotal*sizeof(dfloat),0, 0*mppf->fieldOffset*sizeof(dfloat));
  
  occaTimerTic(mesh->device,"Psi-Solve");
  mppf->NiterPsi = ellipticSolve(psiSolver, mppf->lambdaPsi, mppf->phiTOL, mppf->o_rhsPhi, mppf->o_Psi);
  occaTimerToc(mesh->device,"Psi-Solve"); 


if (mppf->phiOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    mppf->phaseFieldAddBCKernel(mesh->Nelements,
                            time,
                            mesh->o_sgeo,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vmapM,
                            mppf->o_VmapB,
                            mppf->o_Psi);
  }




  // PHASE-FIELD SECOND SOLVE for PHI

  // Compute Rhs for Phi Solve i.e. rhs =  -M*J*Psi
  mppf->phaseFieldRhsSolve2Kernel(mesh->Nelements, mesh->o_vgeo, mesh->o_MM, mppf->o_Psi, mppf->o_rhsPhi);

 if (mppf->phiOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    // NADA



  } else if (mppf->phiOptions.compareArgs("DISCRETIZATION","IPDG")) {
    // Currently we do not need that deuto homegenous bcs need to be used for more complex bcs
    occaTimerTic(mesh->device,"CahnHilliardRhsIpdgBC");    
    mppf->phaseFieldRhsIpdgBCKernel(mesh->Nelements,
                                  mesh->o_vmapM,
                                  phiSolver->tau,
                                  time,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  mesh->o_EToB,
                                  mesh->o_Dmatrices,
                                  mesh->o_LIFTT,
                                  mesh->o_MM,
                                  mppf->o_rhsPhi);
    occaTimerToc(mesh->device,"CahnHilliardRhsIpdgBC");   
  }

  occaTimerTic(mesh->device,"Phi-Solve");
  mppf->NiterPhi = ellipticSolve(phiSolver, mppf->lambdaPhi, mppf->phiTOL, mppf->o_rhsPhi, mppf->o_rkPhi);
  occaTimerToc(mesh->device,"Phi-Solve");

 

  if (mppf->vOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    // NADA




    
  }

  

}
