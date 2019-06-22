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

#include "ins.h"

void insPressureRhs(ins_t *ins, dfloat time, int stage){

  mesh_t *mesh = ins->mesh;

  // int TOMBO = 0;
  // TOMBO =ins->options.compareArgs("TIME INTEGRATOR", "TOMBO") ? 1:0; 

if(ins->TOMBO){
  
  // printf("TOMBO\n");
  // First Compute the F (o_rkU) = sum_j (alfa_j U_j)- sum_j ( beta_j (NU_j + NC_j) )   
  ins->pressureRhsKernel(mesh->Nelements,
                         mesh->o_vgeo,
                         mesh->o_MM,
                         ins->idt,
                         ins->nu,
                         ins->o_extbdfA,
                         ins->o_extbdfB,
                         ins->fieldOffset,
                         ins->o_U,
                         ins->o_NU,
                         ins->o_NC,
                         ins->o_FU,
                         ins->o_rkU); 


  
  // weak divergence term (not that no boundary contribution) -> (F, grad(phi))_sigma - g_0/dt*(n*U^n+1)_del_sigma
  insDivergence(ins, time, ins->o_rkU, ins->o_rhsP);

#if 0
    // ins->o_rkU.copyTo(ins->U);
   ogsGatherScatter(ins->o_rhsP, ogsDfloat, ogsAdd, mesh->ogs);
    ins->o_rhsP.copyTo(ins->P);

    char fname[BUFSIZ];
    string outName;
    sprintf(fname, "div_%04d.vtu",ins->frame++);
    insPlotVTU(ins, fname);
#endif

  
  // Add  -(grad P, grad phi) to rhsP
  const dfloat lambda = 0.0; 
  // simple AX kernel, will be modified later AK.....
  ins->pressureAxKernel(mesh->Nelements,
                        mesh->o_ggeo, 
                        mesh->o_Dmatrices, 
                        mesh->o_Smatrices, 
                        mesh->o_MM, 
                        lambda, 
                        ins->o_P, 
                        ins->o_rhsP); 
  
  // dfloat *test = (dfloat *)calloc(ins->Ntotal, sizeof(dfloat));
  // occa::memory o_test = mesh->device.malloc(ins->Ntotal*sizeof(dfloat), test);

  // ins->o_rhsP.copyFrom(o_test); 
  //  // simple AX kernel, will be modified later AK.....
  // ins->AxKernel(mesh->Nelements,
  //               mesh->o_ggeo, 
  //               mesh->o_Dmatrices, 
  //               mesh->o_Smatrices, 
  //               mesh->o_MM, 
  //               lambda, 
  //               ins->o_P, 
  //               ins->o_rhsP); 


#if 0
    ogsGatherScatter(ins->o_rhsP, ogsDfloat, ogsAdd, mesh->ogs);
    ins->o_rhsP.copyTo(ins->P);

    char fname[BUFSIZ];
    string outName;
    sprintf(fname, "divRhsP_%04d.vtu",ins->frame++);
    insPlotVTU(ins, fname);
#endif


  // Add laplace (Pn)

}else{
  // rhsP = Div Uhat
  insDivergence(ins, time, ins->o_rkU, ins->o_rhsP);
  
  // rhsP = -MM*Div Uhat/pa_ss dt
  //dfloat g0 = 1.0/ins->prkA[stage->ins->Nstages+stage];
  occaTimerTic(mesh->device,"PoissonRhsForcing");
  ins->pressureRhsKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_MM,
                              ins->idt, 
                              ins->g0, 
                              ins->o_rhsP);
  occaTimerToc(mesh->device,"PoissonRhsForcing");

#if 0
  //add penalty from jumps in previous pressure
  ins->poissonPenaltyKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_vgeo,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                ins->tau,
                                mesh->o_x,
                                mesh->o_y,
                                t,
                                ins->dt,
                                ins->c0,
                                ins->c1,
                                ins->c2,
                                ins->index,
                                (mesh->Nelements+mesh->totalHaloPairs),
                                ins->o_P,
                                ins->o_rhsP);
#endif  
  }
}
