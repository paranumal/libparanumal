#include "mppf.h"

void mppfCahnHilliardRhs(mppf_t *mppf, dfloat time){
  
  mesh_t *mesh = mppf->mesh; 

  // Velocity needed to be exchanged, Phi is already done on update function
  if(mesh->totalHaloPairs>0){
    mppf->velocityHaloExtractKernel(mesh->Nelements,
                                    mesh->totalHaloPairs,
                                    mesh->o_haloElementList,
                                    mppf->fieldOffset,
                                    mppf->o_U,
                                    mppf->o_vHaloBuffer);

    // copy extracted halo to HOST 
    mppf->o_vHaloBuffer.copyTo(mppf->vSendBuffer);           

    // start halo exchange
    meshHaloExchangeStart(mesh,
                         mesh->Np*(mppf->NVfields)*sizeof(dfloat),
                         mppf->vSendBuffer,
                         mppf->vRecvBuffer);
  }

  

mppf->phaseFieldAdvectionVolumeKernel(mesh->Nelements,
                                       mesh->o_vgeo,
                                       mesh->o_cubvgeo,
                                       mesh->o_cubDWmatrices,
                                       mesh->o_cubInterpT,
                                       mesh->o_cubProjectT,
                                       mppf->fieldOffset,
                                       mppf->o_U,
                                       mppf->o_cU,
                                       mppf->o_Phi,
                                       mppf->o_NPhi); 


  // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0){
    meshHaloExchangeFinish(mesh);

    mppf->o_vHaloBuffer.copyFrom(mppf->vRecvBuffer); 

    mppf->velocityHaloScatterKernel(mesh->Nelements,
                                    mesh->totalHaloPairs,
                                    mppf->fieldOffset,
                                    mppf->o_U,
                                    mppf->o_vHaloBuffer);
  }

  mppf->phaseFieldAdvectionSurfaceKernel(mesh->Nelements,
                                        mesh->o_vgeo,
                                        mesh->o_sgeo,
                                        mesh->o_cubsgeo,
                                        mesh->o_intInterpT,
                                        mesh->o_intLIFTT,
                                        mesh->o_cubInterpT,
                                        mesh->o_cubProjectT,
                                        mesh->o_vmapM,
                                        mesh->o_vmapP,
                                        mesh->o_EToB,
                                        time,
                                        mesh->o_intx,
                                        mesh->o_inty,
                                        mesh->o_intz,
                                        mppf->fieldOffset,
                                        mppf->o_U,
                                        mppf->o_Phi,
                                        mppf->o_NPhi);


#if 0
// dfloat tn1 = mppf->time + mppf->dt;
  // Compute laplace (h (phi^*)  - S/eta^2 phi^*) (replace with actial IPDG version)
  mppf->phaseFieldDivGradKernel(mesh->Nelements,
                                mesh->o_vgeo,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                mesh->o_Dmatrices,
                                time,
                                mppf->o_extbdfA,
                                mppf->fieldOffset,
                                mppf->inveta2,
                                mppf->chSeta2,
                                mppf->o_Phi,
                                mppf->o_HPhi,
                                mppf->o_lapPhi);
#else
// Upadte Hq energy term
const dlong Ntotal = mesh->Nelements+mesh->totalHaloPairs;
mppf->phaseFieldUpdateMixingEnergyKernel(Ntotal,
                                         time,
                                         mppf->inveta2,
                                         mppf->chSeta2,
                                         mppf->o_Phi,
                                         mppf->o_HPhi);




// Extrapolate and compute Gradient
mppf->phaseFieldAxGradKernel(mesh->Nelements,
                            mppf->fieldOffset,
                            mppf->o_extbdfA,
                            mesh->o_vgeo,
                            mesh->o_Dmatrices,
                            mppf->o_HPhi,
                            mppf->o_GHPhi);




// Velocity needed to be exchanged, Phi is already done on update function
  if(mesh->totalHaloPairs>0){
    const int gNfields = 4; 
    mppf->gradPhaseFieldHaloExtractKernel(mesh->Nelements,
                                          mesh->totalHaloPairs,
                                          mesh->o_haloElementList,
                                          mppf->o_GHPhi,
                                          mppf->o_gPhiHaloBuffer);

    // copy extracted halo to HOST 
    mppf->o_gPhiHaloBuffer.copyTo(mppf->gPhiSendBuffer);           

    // start halo exchange
    meshHaloExchangeStart(mesh,
                         mesh->Np*gNfields*sizeof(dfloat),
                         mppf->gPhiSendBuffer,
                         mppf->gPhiRecvBuffer);


    meshHaloExchangeFinish(mesh);

    mppf->o_gPhiHaloBuffer.copyFrom(mppf->gPhiRecvBuffer); 

    mppf->gradPhaseFieldHaloScatterKernel(mesh->Nelements,
                                          mesh->totalHaloPairs,
                                          mppf->o_GHPhi,
                                          mppf->o_gPhiHaloBuffer);
  }




const dfloat zero = 0.0; 

// laplace(HPhi) = - M^-1 * Ax(HPhi)
mppf->phaseFieldAxIpdgKernel(mesh->Nelements,
                        mesh->o_vmapM,
                        mesh->o_vmapP,
                        zero, // lambda,
                        mppf->psiSolver->tau,// zero, // mppf->tau,
                        mesh->o_vgeo,
                        mesh->o_sgeo,
                        mesh->o_EToB,
                        mesh->o_Dmatrices,
                        mesh->o_LIFTT,
                        mesh->o_MM,
                        mppf->o_GHPhi,
                        mppf->o_lapPhi);
#endif

  
  
  mppf->phaseFieldRhsSolve1Kernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_MM,
                            time,
                            mppf->idt,
                            mppf->o_extbdfA,
                            mppf->o_extbdfB,
                            mppf->fieldOffset,
                            mppf->o_Phi,
                            mppf->o_NPhi,
                            mppf->o_lapPhi,
                            mppf->o_rhsPhi);

}
