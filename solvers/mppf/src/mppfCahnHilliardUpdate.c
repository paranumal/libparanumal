#include "mppf.h"

void mppfCahnHilliardUpdate(mppf_t *mppf, dfloat time){
  
 mesh_t *mesh = mppf->mesh; 
 
  // Update History// Very bad storage!!!! will try to decrease
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

// Feed exact Phi and Psi
#if 0
  // Set pahse field function on device
  mppf->setPhaseFieldKernel(mesh->Nelements,
                          time,
                          mppf->eta,
                          mesh->o_x,
                          mesh->o_y,
                          mesh->o_z,
                          mppf->o_rkPhi);

  for(int e=0; e<mesh->Nelements;e++){
    for(int n=0; n<mesh->Np; n++){
      const int id = e*mesh->Np + n;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];
      // mppf->Psi[id] =   -2*M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(time) + mppf->chA*cos(M_PI*x)*cos(M_PI*y)*sin(time);
      mppf->Psi[id] =   -2.0*M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(time);
    }
  }
mppf->o_Psi.copyFrom(mppf->Psi);
#endif

  

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


// Give exact GPhi
 #if 0

for(int e=0; e<mesh->Nelements;e++){
    for(int n=0; n<mesh->Np; n++){
      const int id = e*mesh->Np + n;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];
      //
      dfloat phix   = -M_PI*cos(M_PI*y)*sin(M_PI*x)*sin(time);
      dfloat phiy   = -M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(time);
      //
      mppf->rkU[id + 0*mppf->fieldOffset] = phix;
      mppf->rkU[id + 1*mppf->fieldOffset] = phiy;
    }
  }
  mppf->o_GPhi.copyFrom(mppf->rkU);

 #endif


  

  // Smooth density and viscosity on device 
  mppf->setMaterialPropertyKernel(mesh->Nelements, mppf->o_Phi, mppf->o_Rho, mppf->o_Mu);


}