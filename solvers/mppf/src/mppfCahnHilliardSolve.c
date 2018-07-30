#include "mppf.h"

// solve lambda*U + A*U = rhsU
void mppfCahnHilliardSolve(mppf_t *mppf, dfloat time, occa::memory o_PhiHat){
  
  mesh_t *mesh = mppf->mesh; 
  elliptic_t *phiSolver = mppf->phiSolver; 
  elliptic_t *psiSolver = mppf->psiSolver; 
  
  if (mppf->phiOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    // mppf->velocityRhsBCKernel(mesh->Nelements,
    //                           mesh->o_ggeo,
    //                           mesh->o_sgeo,
    //                           mesh->o_Dmatrices,
    //                           mesh->o_Smatrices,
    //                           mesh->o_MM,
    //                           mesh->o_vmapM,
    //                           mesh->o_sMT,
    //                           mppf->lambda,
    //                           time,
    //                           mesh->o_x,
    //                           mesh->o_y,
    //                           mesh->o_z,
    //                           mppf->o_VmapB,
    //                           o_rhsU,
    //                           o_rhsV,
    //                           o_rhsW);
    
    // // gather-scatter
    // ellipticParallelGatherScatter(mesh, mesh->ogs, o_rhsU, dfloatString, "add");  
    // ellipticParallelGatherScatter(mesh, mesh->ogs, o_rhsV, dfloatString, "add");  
    // if (mppf->dim==3)
    //   ellipticParallelGatherScatter(mesh, mesh->ogs, o_rhsW, dfloatString, "add");  
    // if (usolver->Nmasked) mesh->maskKernel(usolver->Nmasked, usolver->o_maskIds, o_rhsU);
    // if (vsolver->Nmasked) mesh->maskKernel(vsolver->Nmasked, vsolver->o_maskIds, o_rhsV);
    // if (mppf->dim==3)
    //   if (wsolver->Nmasked) mesh->maskKernel(wsolver->Nmasked, wsolver->o_maskIds, o_rhsW);

  } else if (mppf->phiOptions.compareArgs("DISCRETIZATION","IPDG")) {

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

  //copy current velocity fields as initial guess? (could use Uhat or beter guess)
  dlong Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
  mppf->o_PhiH.copyFrom(mppf->o_Phi,Ntotal*sizeof(dfloat),0, 0*mppf->fieldOffset*sizeof(dfloat));
  

  // if (mppf->vOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
  //   if (usolver->Nmasked) mesh->maskKernel(usolver->Nmasked, usolver->o_maskIds, mppf->o_UH);
  //   if (vsolver->Nmasked) mesh->maskKernel(vsolver->Nmasked, vsolver->o_maskIds, mppf->o_VH);
  //   if (mppf->dim==3)
  //     if (wsolver->Nmasked) mesh->maskKernel(wsolver->Nmasked, wsolver->o_maskIds, mppf->o_WH);

  // }
  
  occaTimerTic(mesh->device,"Psi-Solve");
  mppf->NiterPsi = ellipticSolve(psiSolver, mppf->lambdaPsi, mppf->phiTOL, mppf->o_rhsPhi, mppf->o_Psi);
  occaTimerToc(mesh->device,"Psi-Solve"); 
  
  // Compute Rhs for Phi Solve i.e. rhs =  -M*J*Psi
  mppf->phaseFieldRhsSolve2Kernel(mesh->Nelements, mesh->o_vgeo, mesh->o_MM, mppf->o_Psi, mppf->o_rhsPhi);

  occaTimerTic(mesh->device,"Phi-Solve");
  mppf->NiterPhi = ellipticSolve(phiSolver, mppf->lambdaPhi, mppf->phiTOL, mppf->o_rhsPhi, mppf->o_PhiH);
  occaTimerToc(mesh->device,"Phi-Solve");

 

  if (mppf->vOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    // mppf->velocityAddBCKernel(mesh->Nelements,
    //                         time,
    //                         mesh->o_sgeo,
    //                         mesh->o_x,
    //                         mesh->o_y,
    //                         mesh->o_z,
    //                         mesh->o_vmapM,
    //                         mppf->o_VmapB,
    //                         mppf->o_UH,
    //                         mppf->o_VH,
    //                         mppf->o_WH);
  }

  //copy into intermediate stage storage
   mppf->o_PhiH.copyTo(o_PhiHat,Ntotal*sizeof(dfloat),0*mppf->fieldOffset*sizeof(dfloat),0);
}
