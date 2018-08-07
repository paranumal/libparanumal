#include "mppf.h"

void mppfCahnHilliardRhs(mppf_t *mppf, dfloat time){
  
  mesh_t *mesh = mppf->mesh; 

  // Velocity needed to be axchange, Phi is already done on update function

  //Exctract Halo On Device for Phase Field function it is already done after Solve!!!!!
  // // Assumes Velocity is already halo exchanged  
  // if(mesh->totalHaloPairs>0){
  //   mppf->phaseFieldHaloExtractKernel(mesh->Nelements,
  //                                     mesh->totalHaloPairs,
  //                                     mesh->o_haloElementList,
  //                                     mppf->o_Phi,
  //                                     mppf->o_phiHaloBuffer);

  //   // copy extracted halo to HOST 
  //   mppf->o_phiHaloBuffer.copyTo(mppf->phiSendBuffer);           
  
  //   // start halo exchange
  //   meshHaloExchangeStart(mesh,
  //                        mesh->Np*sizeof(dfloat),
  //                        mppf->phiSendBuffer,
  //                        mppf->phiRecvBuffer);
  // }
  
  // 
  // 1-> compute NPhi =  u*grad(Phi) = div(u*Phi) on Cubature Nodes and update potential function HPhi
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
                                       mppf->o_NPhi, // Nonlinear convective Cahn-Hilliard term
                                       mppf->o_HPhi); // Potential function to be extrapolated: double-well currently


  // if(mesh->totalHaloPairs>0){

  //   meshHaloExchangeFinish(mesh);

  //   mppf->o_phiHaloBuffer.copyFrom(mppf->phiRecvBuffer);

  //   mppf->phaseFieldHaloScatterKernel(mesh->Nelements,
  //                                     mesh->totalHaloPairs,
  //                                     mppf->o_Phi,
  //                                     mppf->o_phiHaloBuffer);
  // }

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
