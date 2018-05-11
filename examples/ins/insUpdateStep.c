#include "ins.h"

void insUpdateStep(ins_t *ins, dfloat time){

  mesh_t *mesh = ins->mesh;
  //dfloat t = tstep*ins->dt + ins->dt;

  dlong offset = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;

  if (ins->pOptions.compareArgs("DISCRETIZATION","IPDG")) {
    if(mesh->totalHaloPairs>0){

      ins->pressureHaloExtractKernel(mesh->Nelements,
                                 mesh->totalHaloPairs,
                                 mesh->o_haloElementList,
                                 ins->o_PI,
                                 ins->o_pHaloBuffer);

      // copy extracted halo to HOST
      ins->o_pHaloBuffer.copyTo(ins->pSendBuffer);

      // start halo exchange
      meshHaloExchangeStart(mesh,
                           mesh->Np*sizeof(dfloat),
                           ins->pSendBuffer,
                           ins->pRecvBuffer);
    }
  }
  
  occaTimerTic(mesh->device,"GradientVolume");
  // Compute Volume Contribution of gradient of pressure gradient
  ins->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_Dmatrices,
                            offset,
                            ins->o_PI,
                            ins->o_gradPI);
  occaTimerToc(mesh->device,"GradientVolume");

  if (ins->pOptions.compareArgs("DISCRETIZATION","IPDG")) {
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);

      ins->o_pHaloBuffer.copyFrom(ins->pRecvBuffer);

      ins->pressureHaloScatterKernel(mesh->Nelements,
                                      mesh->totalHaloPairs,
                                      ins->o_PI,
                                      ins->o_pHaloBuffer);
    }
    
    const int solverid =1 ;

    occaTimerTic(mesh->device,"GradientSurface");
    // Compute Surface Contribution of gradient of pressure increment
    ins->gradientSurfaceKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_LIFTT,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                time,
                                offset,
                                solverid, // pressure increment BCs
                                ins->o_PI,
                                ins->o_gradPI);
    
    occaTimerToc(mesh->device,"GradientSurface");
  }

  
  // U <= U - dt/g0 * grad(pressure increment)
  occaTimerTic(mesh->device,"UpdateUpdate");
  ins->updateUpdateKernel(mesh->Nelements,
                              ins->dt,
                              ins->ig0,
			                        ins->c0,
                              ins->c1,
                              ins->c2,
                              offset,
                              ins->o_Uhat,
                              ins->o_PI,
                              ins->o_gradPI,
                              ins->o_U,
                              ins->o_NU,
                              ins->o_P,
                              ins->o_gradP);
  occaTimerToc(mesh->device,"UpdateUpdate");
}
