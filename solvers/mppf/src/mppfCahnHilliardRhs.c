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




  


// Upadte Hq energy term
const dlong Ntotal = mesh->Nelements+mesh->totalHaloPairs;
dfloat tt = time-mppf->dt;
mppf->phaseFieldUpdateMixingEnergyKernel(Ntotal,
                                         tt,
                                         mppf->inveta2,
                                         mppf->chSeta2,
                                         mesh->o_cubInterpT,
                                         mesh->o_cubProjectT,
                                         mesh->o_x,
                                         mesh->o_y,
                                         mesh->o_z,
                                         mppf->o_Phi,
                                         mppf->o_HPhi);


mppf->phaseFieldExtrapolateKernel(Ntotal, // Extrapolate the halo elements, no need to exchange!!!
                                  mppf->Nstages,
                                  mppf->fieldOffset,
                                  mppf->o_extbdfA,
                                  mppf->o_HPhi,
                                  mppf->o_HPhiE); // Extrapolated mixing energy


if(mppf->phiOptions.compareArgs("DISCRETIZATION", "IPDG")){

  mppf->phaseFieldAxGradVolumeKernel(mesh->Nelements,
                                  mesh->o_vgeo,
                                  mesh->o_Dmatrices,
                                  mppf->fieldOffset,
                                  mppf->o_HPhiE,
                                  mppf->o_GPhi);



   mppf->phaseFieldAxGradSurfaceKernel(mesh->Nelements,
                                      mppf->inveta2,
                                      mppf->chSeta2,
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
                                      mppf->o_HPhiE,
                                      mppf->o_GU,
                                      mppf->o_GPhi);

   // Just Testing postprocessing
#if 0
  dlong NtotalT = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
  for(int s=0; s<mppf->dim; s++){
    mppf->o_rkP.copyFrom(mppf->o_GPhi,NtotalT*sizeof(dfloat),0,s*mppf->fieldOffset*sizeof(dfloat));

    ogsGatherScatter(mppf->o_rkP, ogsDfloat, ogsAdd, mesh->ogs);  
    mppf->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, mppf->o_rkP, mppf->o_rkP);
    
    mppf->o_rkP.copyTo(mppf->o_GPhi,NtotalT*sizeof(dfloat),s*mppf->fieldOffset*sizeof(dfloat),0);
  }
#endif

const dfloat zero = 0.0; 
    // laplace(HPhi) = - M^-1 * Ax(HPhi)
  mppf->phaseFieldAxIpdgKernel(mesh->Nelements,
                        mesh->o_vmapM,
                        mesh->o_vmapP,
                        time,
                        mppf->inveta2,
                        mppf->chSeta2,
                        zero, // lambda,
                        zero, //mppf->psiSolver->tau, // zero ????
                        mppf->fieldOffset,
                        mesh->o_vgeo,
                        mesh->o_sgeo,
                        mesh->o_x,
                        mesh->o_y,
                        mesh->o_z,
                        mesh->o_EToB,
                        mesh->o_Dmatrices,
                        mesh->o_LIFTT,
                        mesh->o_MM,
                        mppf->o_HPhiE,
                        mppf->o_GPhi,
                        mppf->o_lapPhi);
}


if (mppf->phiOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {

mppf->phaseFieldAxKernel(mesh->Nelements,
                        mesh->o_ggeo,
                        time,          // no need just for exact solution
                        mppf->inveta2, // no need just for exact solution
                        mppf->chSeta2, // no need just for exact solution
                        mesh->o_Dmatrices,
                        mesh->o_Smatrices,
                        mppf->o_HPhiE,
                        mppf->o_lapPhi);
}




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