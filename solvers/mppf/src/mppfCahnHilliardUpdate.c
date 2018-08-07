#include "mppf.h"

void mppfCahnHilliardUpdate(mppf_t *mppf, dfloat time){
  
 mesh_t *mesh = mppf->mesh; 
 
  // Update History
  for (int s=mppf->Nstages;s>1;s--) {
    mppf->o_Phi.copyFrom(mppf->o_Phi, mppf->Ntotal*sizeof(dfloat), 
                                  (s-1)*mppf->Ntotal*sizeof(dfloat), 
                                  (s-2)*mppf->Ntotal*sizeof(dfloat));

    mppf->o_NPhi.copyFrom(mppf->o_NPhi, mppf->Ntotal*sizeof(dfloat), 
                                  (s-1)*mppf->Ntotal*sizeof(dfloat), 
                                  (s-2)*mppf->Ntotal*sizeof(dfloat));

    mppf->o_HPhi.copyFrom(mppf->o_HPhi, mppf->Ntotal*sizeof(dfloat), 
                                  (s-1)*mppf->Ntotal*sizeof(dfloat), 
                                  (s-2)*mppf->Ntotal*sizeof(dfloat));
  }

  // Update Phi
  mppf->o_Phi.copyFrom(mppf->o_rkPhi, mppf->Ntotal*sizeof(dfloat));

  // Compute grad(Phi)

  if(mesh->totalHaloPairs>0){
    mppf->phaseFieldHaloExtractKernel(mesh->Nelements,
                                      mesh->totalHaloPairs,
                                      mesh->o_haloElementList,
                                      mppf->o_Phi,
                                      mppf->o_phiHaloBuffer);

    // copy extracted halo to HOST 
    mppf->o_phiHaloBuffer.copyTo(mppf->phiSendBuffer);           
  
    // start halo exchange
    meshHaloExchangeStart(mesh,
                         mesh->Np*sizeof(dfloat),
                         mppf->phiSendBuffer,
                         mppf->phiRecvBuffer);
  }


  mppf->phaseFieldGradientVolumeKernel(mesh->Nelements,
                                      mesh->o_vgeo,
                                      mesh->o_Dmatrices,
                                      mppf->fieldOffset,
                                      mppf->o_Phi,
                                      mppf->o_GPhi);

   if(mesh->totalHaloPairs>0){

    meshHaloExchangeFinish(mesh);

    mppf->o_phiHaloBuffer.copyFrom(mppf->phiRecvBuffer);

    mppf->phaseFieldHaloScatterKernel(mesh->Nelements,
                                      mesh->totalHaloPairs,
                                      mppf->o_Phi,
                                      mppf->o_phiHaloBuffer);
  }


 mppf->phaseFieldGradientSurfaceKernel(mesh->Nelements,
                                      mesh->o_sgeo,
                                      mesh->o_LIFTT,
                                      mesh->o_vmapM,
                                      mesh->o_vmapP,
                                      mesh->o_EToB,
                                      mesh->o_x,
                                      mesh->o_y,
                                      mesh->o_z,
                                      time,
                                      mppf->fieldOffset,
                                      mppf->o_Phi,
                                      mppf->o_GPhi);


  

  // Smooth density and viscosity on device 
  mppf->setMaterialPropertyKernel(mesh->Nelements, mppf->o_Phi, mppf->o_Rho, mppf->o_Mu);


  
  // Compute Gradient of Mu

  mppf->phaseFieldGradientVolumeKernel(mesh->Nelements,
                                      mesh->o_vgeo,
                                      mesh->o_Dmatrices,
                                      mppf->fieldOffset,
                                      mppf->o_Mu,
                                      mppf->o_GMu);






}