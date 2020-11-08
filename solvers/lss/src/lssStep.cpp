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

#include "lss.hpp"

//evaluate ODE rhs = f(q,t)
void lss_t::rhsf(occa::memory& o_Q, occa::memory& o_RHS, const dfloat T){
  // Solve Level-Set Advection  
  if(advection){
    Advection(o_Q, o_RHS, T);
  }
  // Solve Redistance Problem
  if(redistance){
   Redistance(o_Q, o_RHS, T); 
  }

}




//evaluate ODE rhs = f(q,t)
void lss_t::rhsf_subcell(occa::memory& o_Q, occa::memory &o_sQ, occa::memory& o_RHS,occa::memory& o_sRHS, const dfloat T){
  
   Redistance_Subcell(o_Q, o_sQ, o_RHS, o_sRHS, T); 
}


void lss_t::Advection(occa::memory& o_Q, occa::memory& o_RHS, const dfloat T){
 
// extract q halo on DEVICE
  traceHalo->ExchangeStart(o_Q, 1, ogs_dfloat);

  setFlowFieldKernel(mesh.Nelements, 
                    offset, 
                    T, 
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    o_U,
                    o_Q);

  if(cubature){
     advectionVolumeKernel(mesh.Nelements,
                         offset, 
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubDWmatrices,
                         mesh.o_cubInterpT,
                         mesh.o_cubProjectT,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         T,
                         o_U,
                         o_Q,
                         o_RHS);
  }else{
  advectionVolumeKernel(mesh.Nelements,
               offset, 
               mesh.o_vgeo,
               mesh.o_Dmatrices,
               T,
               mesh.o_x,
               mesh.o_y,
               mesh.o_z,
               o_U,
               o_Q,
               o_RHS);
 }

  traceHalo->ExchangeFinish(o_Q, 1, ogs_dfloat);
  
  if(cubature){
    advectionSurfaceKernel(mesh.Nelements,
                            offset, 
                            mesh.o_vgeo,
                            mesh.o_cubsgeo,
                            mesh.o_vmapM,
                            mesh.o_vmapP,
                            mesh.o_EToB,
                            mesh.o_intInterpT,
                            mesh.o_intLIFTT,
                            mesh.o_intx,
                            mesh.o_inty,
                            mesh.o_intz,
                            T,
                            o_U,
                            o_Q,
                            o_RHS);

  }else{
  advectionSurfaceKernel(mesh.Nelements,
                offset, 
                mesh.o_sgeo,
                mesh.o_LIFTT,
                mesh.o_vmapM,
                mesh.o_vmapP,
                mesh.o_EToB,
                T,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                o_U,
                o_Q,
                o_RHS);
 }
}



void lss_t::Redistance(occa::memory& o_Q, occa::memory& o_RHS, const dfloat T){

    traceHalo->ExchangeStart(o_Q, 1, ogs_dfloat);

    redistanceVolumeKernel(mesh.Nelements,
               offset, 
               mesh.o_vgeo,
               mesh.o_DWmatrices, // o_Dmatrices
               T,
               mesh.o_x,
               mesh.o_y,
               mesh.o_z,
               o_Q,
               o_gradq);


  traceHalo->ExchangeFinish(o_Q, 1, ogs_dfloat);
  
  redistanceSurfaceKernel(mesh.Nelements,
                offset, 
                mesh.o_sgeo,
                mesh.o_LIFTT,
                mesh.o_vmapM,
                mesh.o_vmapP,
                mesh.o_EToB,
                T,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                o_sgnq,
                o_Q,
                o_gradq,
                o_RHS);
  }


// this function reconstruct DG solution from FV cells and addd to DG storage,q
void lss_t::reconstruct_subcell(occa::memory& o_Q, occa::memory& o_sQ){

int all = 0; 
reconstructKernel(mesh.Nelements, 
                    all,
                    subcell->o_ElementList,
                    subcell->o_RMT,
                    o_sQ, 
                    o_Q);



}

void lss_t::Redistance_Subcell(occa::memory& o_Q, occa::memory & o_sQ, 
                                occa::memory& o_RHS, occa::memory& o_sRHS, const dfloat T){

    // extract q halo on DEVICE
    // traceHalo->ExchangeStart(o_Q, 1, ogs_dfloat);

    DetectTroubledCells(o_Q, subcell->o_ElementList); 

    mesh.halo->ExchangeStart(o_Q, mesh.Np, ogs_dfloat); 

    // Compute dq/dx for DG and DGFV elements
     partialRedistanceVolumeKernel(mesh.Nelements,
                                   subcell->o_ElementList,
                                   offset, 
                                   mesh.o_vgeo,
                                   mesh.o_DWmatrices, // o_Dmatrices
                                   T,
                                   mesh.o_x,
                                   mesh.o_y,
                                   mesh.o_z,
                                   o_Q,
                                   o_gradq);

     const int all = 0; // all elements->do not use element list
     // Project DG solution to FV face values directly for DGFV type
     projectDGKernel(mesh.Nelements,
                    all, 
                    subcell->o_ElementList, 
                    mesh.o_vmapM,
                    subcell->o_mFToE,
                    subcell->o_mFToF,
                    subcell->o_PFMT,
                    o_Q,
                    o_sface); 


    mesh.halo->ExchangeFinish(o_Q, mesh.Np, ogs_dfloat); 


     // Add the Contribution of DG Elements
     partialRedistanceSurfaceKernel(mesh.Nelements,
                offset, 
                subcell->o_ElementList,
                mesh.o_sgeo,
                mesh.o_LIFTT,
                mesh.o_vmapM,
                mesh.o_vmapP,
                mesh.o_EToB,
                T,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                o_sgnq,
                o_Q,
                o_gradq,
                o_RHS);



    // Compute Cell-Centered Subcell Values 
    projectKernel(mesh.Nelements + mesh.totalHaloPairs,
                  all,  
                  subcell->o_ElementList, 
                  subcell->o_PMT,
                  o_Q, 
                  o_sQ);

       
#if 0 
    // Reconstruct internal subcell FV face values using WENO2:FVFV elements
     reconstructInternalFaceKernel(mesh.Nelements, 
                                  subcell->o_ElementList,
                                  subcell->o_ielist,
                                  subcell->o_vgeo,
                                  subcell->o_sgeo, 
                                  subcell->o_emapP,
                                  subcell->o_fmapP, 
                                  o_Q, 
                                  o_sQ,
                                  o_sface);  


     // Reconstruct external subcell FV face values using WENO2
     reconstructExternalFaceKernel(mesh.Nelements, 
                                  subcell->o_ElementList,
                                  subcell->o_eelist,
                                  subcell->o_vgeo,
                                  subcell->o_sgeo, 
                                  subcell->o_emapP,
                                  subcell->o_fmapP, 
                                  o_Q, 
                                  o_sQ,
                                  o_sface);   
#else

 // Reconstruct face values for all subcells 
     reconstructFaceKernel(mesh.Nelements, 
                          subcell->o_ElementList,
                          subcell->o_vgeo,
                          subcell->o_sgeo, 
                          subcell->o_emapP,
                          subcell->o_fmapP, 
                          o_Q, 
                          o_sQ,
                          o_sface);   



#endif


     mesh.halo->Exchange(o_sface, subcell->Nsubcells*mesh.Nfaces, ogs_dfloat);

    // FV compute 
     subcellComputeKernel(mesh.Nelements, 
                              subcell->o_ElementList,
                              subcell->o_emapP,
                              subcell->o_fmapP, 
                              subcell->o_RMT,
                              subcell->o_vgeo,
                              subcell->o_sgeo, 
                              o_Q,
                              o_ssgnq,
                              o_sface, 
                              o_sRHS);  

      mixedRedistanceSurfaceKernel(mesh.Nelements,
                          offset, 
                          subcell->o_ElementList,
                          subcell->o_EToE, 
                          mesh.o_EToB,
                          mesh.o_sgeo,
                          subcell->o_sgeo, 
                          subcell->o_vgeo,
                          mesh.o_LIFTT,
                          subcell->o_SLIFTT,
                          subcell->o_RFMT,
                          mesh.o_vmapM,
                          mesh.o_vmapP,
                          subcell->o_emapP,
                          subcell->o_fmapP, 
                          subcell->o_mFToE,
                          subcell->o_mFToF,
                          subcell->o_mDGID,
                          T,
                          mesh.o_x,
                          mesh.o_y,
                          mesh.o_z,
                          o_sgnq,
                          o_Q,
                          o_sface,
                          o_gradq,
                          o_RHS);
}


void lss_t::DetectTroubledCells(occa::memory& o_Q, occa::memory& o_ElementList){

// Currently modal decay + skyline only
if(subcellStabilization){

  if(settings.compareSetting("INDICATOR TYPE", "MDA")){
  indicatorMDAKernel(mesh.Nelements,
                mesh.o_ggeo,
                subcell->o_ModMap, 
                mesh.o_MM,
                subcell->o_invVT,
                subcell->o_LSF,
                subcell->o_BLD, 
                o_Q,
                o_ElementList); 


  }else if(settings.compareSetting("INDICATOR TYPE", "MDH")){
  indicatorMDHKernel(mesh.Nelements,
                mesh.o_ggeo,
                mesh.o_MM,
                subcell->o_interpNM1T,
                o_Q,
                o_ElementList); 
  }else{
    printf("The detector type is not impmented\n");
  }

  // exchange element info before 
  mesh.halo->Exchange(o_ElementList, 1, ogs_dlong);

  // Find neighs i.e. DG neighs of FV elements....
  findNeighKernel(mesh.Nelements,
                 subcell->o_EToE, 
                 o_ElementList);     
}

}
