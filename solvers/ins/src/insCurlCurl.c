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

void insCurlCurl(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_NC){

  mesh_t *mesh = ins->mesh;
  // No that UH, VH and WH is used as the temporary arrays to store curl(u)
  // Note that multiplied with Mass Matrix i.e. JW to prepare smoothing !!!!!
  ins->curlKernel(mesh->Nelements,
                  mesh->o_vgeo,
                  mesh->o_Dmatrices,
                  mesh->o_MM, 
                  ins->fieldOffset,
                  o_U,
                  ins->o_UH, // Wx
                  ins->o_VH, // Wy
                  ins->o_WH);// Wz

  // Gather-Scatter
  if(ins->dim==3){    
   ogsGatherScatter(ins->o_UH, ogsDfloat, ogsAdd, mesh->ogs);
   ogsGatherScatter(ins->o_VH, ogsDfloat, ogsAdd, mesh->ogs);
  } 
   ogsGatherScatter(ins->o_WH, ogsDfloat, ogsAdd, mesh->ogs);


  int nfield = 0; 
  if (ins->dim==3){
    nfield = 3; 
    ins->o_UH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),0*ins->fieldOffset*sizeof(dfloat),0);
    ins->o_VH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),1*ins->fieldOffset*sizeof(dfloat),0);    
    ins->o_WH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),2*ins->fieldOffset*sizeof(dfloat),0);    
  }else{
    nfield =1; 
    // if quad or tri copy Wz to first place
    ins->o_WH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),0*ins->fieldOffset*sizeof(dfloat),0); 
  }
  
 
  // Multiply With Inverse Mass Matrix
   ins->invMassMatrixKernel(mesh->Nelements,
                      ins->fieldOffset,
                      nfield,
                      mesh->o_vgeo,
                      ins->o_InvM, // mesh->o_MM, // should be invMM for tri/tet
                      ins->o_rkU);


  if(ins->dim==2){
  // Second curl on smoothed curl(u)
  ins->curlBKernel(mesh->Nelements,
                  mesh->o_vgeo,
                  mesh->o_Dmatrices,
                  mesh->o_MM, 
                  ins->fieldOffset,
                  ins->o_rkU,
                  ins->o_UH, // Wx
                  ins->o_VH, // Wy
                  ins->o_WH);// Wz
  }else{
    ins->curlKernel(mesh->Nelements,
                  mesh->o_vgeo,
                  mesh->o_Dmatrices,
                  mesh->o_MM, 
                  ins->fieldOffset,
                  ins->o_rkU,
                  ins->o_UH, // Wx
                  ins->o_VH, // Wy
                  ins->o_WH);// Wz
  }
  

  

  //  
   // Gather-Scatter
  ogsGatherScatter(ins->o_UH, ogsDfloat, ogsAdd, mesh->ogs);
  ogsGatherScatter(ins->o_VH, ogsDfloat, ogsAdd, mesh->ogs);
  if(ins->dim==3)       
   ogsGatherScatter(ins->o_WH, ogsDfloat, ogsAdd, mesh->ogs);


  //
  nfield = 0; 
  if (ins->dim==3){
    nfield = 3; 
    ins->o_UH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),0*ins->fieldOffset*sizeof(dfloat),0);
    ins->o_VH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),1*ins->fieldOffset*sizeof(dfloat),0);    
    ins->o_WH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),2*ins->fieldOffset*sizeof(dfloat),0);    
  }else{
    nfield =2; 
    ins->o_UH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),0*ins->fieldOffset*sizeof(dfloat),0);
    ins->o_VH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),1*ins->fieldOffset*sizeof(dfloat),0); 
  }
  //  Multiply With Inverse Mass Matrix
   ins->invMassMatrixKernel(mesh->Nelements,
                            ins->fieldOffset,
                            nfield,
                            mesh->o_vgeo,
                            ins->o_InvM, // mesh->o_MM, // should be invMM for tri/tet
                            ins->o_rkU);


  ins->o_rkU.copyTo(o_NC,ins->NVfields*ins->Ntotal*sizeof(dfloat),0*ins->fieldOffset*sizeof(dfloat),0);
}
