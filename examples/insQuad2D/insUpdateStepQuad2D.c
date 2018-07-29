/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "insQuad2D.h"

void insUpdateStepQuad2D(ins_t *ins, int tstep, char   * options){

  mesh2D *mesh = ins->mesh;
  dfloat t = tstep*ins->dt + ins->dt;

  dlong offset = (mesh->Nelements+mesh->totalHaloPairs);

  if (strstr(ins->pSolverOptions,"IPDG")) {
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
  dlong ioffset = 0;
  ins->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_D,
                            ioffset,  //no offset
                            ins->o_PI,
                            ins->o_PIx,
                            ins->o_PIy);
  occaTimerToc(mesh->device,"GradientVolume");

  if (strstr(ins->pSolverOptions,"IPDG")) {
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);

      ins->o_pHaloBuffer.copyFrom(ins->pRecvBuffer);

      ins->pressureHaloScatterKernel(mesh->Nelements,
                                      mesh->totalHaloPairs,
                                      mesh->o_haloElementList,
                                      ins->o_PI,
                                      ins->o_pHaloBuffer);
    }
    
    const int solverid =1 ;

    occaTimerTic(mesh->device,"GradientSurface");
    // // Compute Surface Contribution of gradient of pressure increment
    ins->gradientSurfaceKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                mesh->o_x,
                                mesh->o_y,
                                t,
                                ins->dt,
                                ins->c0,
                                ins->c1,
                                ins->c2,
                                ins->index,
                                offset,
                                solverid, // pressure increment BCs
                                ins->o_P,
                                ins->o_PI,
                                ins->o_PIx,
                                ins->o_PIy);
    
    occaTimerToc(mesh->device,"GradientSurface");
  }

  // U <= U - dt/g0 * d(pressure increment)/dx
  // V <= V - dt/g0 * d(pressure increment)/dy

  // occaTimerTic(mesh->device,"UpdateUpdate");
  ins->updateUpdateKernel(mesh->Nelements,
                              ins->dt,
                              ins->ig0,
                              ins->a0,
                              ins->a1,
                              ins->a2,
                              ins->c0,
                              ins->c1,
                              ins->c2,
                              ins->o_PI,
                              ins->o_PIx,
                              ins->o_PIy,
                              ins->index,
                              offset,
                              ins->o_U,
                              ins->o_V,
                              ins->o_P);
  occaTimerToc(mesh->device,"UpdateUpdate");

  ins->index = (ins->index+1)%3; //hard coded for 3 stages
}
