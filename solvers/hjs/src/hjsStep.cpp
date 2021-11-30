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

#include "hjs.hpp"
void hjs_t::rhsf(occa::memory& o_Q, occa::memory& o_RHS, const dfloat T){

// extract q halo on DEVICE
  traceHalo->ExchangeStart(o_Q, 1, ogs_dfloat);

   volumeKernel(mesh.Nelements,
               mesh.o_vgeo,
               mesh.o_DW,
               T,
               mesh.o_x,
               mesh.o_y,
               mesh.o_z,
               o_Q,
               o_gradq);

   traceHalo->ExchangeFinish(o_Q, 1, ogs_dfloat);


   surfaceKernel(mesh.Nelements,
                mesh.o_sgeo,
                mesh.o_LIFT,
                mesh.o_vmapM,
                mesh.o_vmapP,
                mesh.o_EToB,
                T,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                o_Q,
                o_gradq);

   // surfaceKernel(mesh.Nelements,
   //              mesh.o_sgeo,
   //              mesh.o_LIFT,
   //              mesh.o_vmapM,
   //              mesh.o_vmapP,
   //              mesh.o_EToB,
   //              T,
   //              mesh.o_x,
   //              mesh.o_y,
   //              mesh.o_z,
   //              o_Q,
   //              o_gradq, 
   //              o_RHS);


 hamiltonianKernel(mesh.Nelements,
               mesh.o_vgeo,
               T,
               mesh.o_x,
               mesh.o_y,
               mesh.o_z,
               o_gradq,
               o_RHS);


}


// //evaluate ODE rhs = f(q,t)
// void lss_t::rhsf(occa::memory& o_Q, occa::memory& o_RHS, const dfloat T){
//   // Solve Level-Set Advection  
//   if(advection){
//     Advection(o_Q, o_RHS, T);
//   }
//   // Solve Redistance Problem
//   if(redistance){
//    Redistance(o_Q, o_RHS, T); 
//   }

// }


// void lss_t::Advection(occa::memory& o_Q, occa::memory& o_RHS, const dfloat T){
 
// // extract q halo on DEVICE
//   traceHalo->ExchangeStart(o_Q, 1, ogs_dfloat);

//   setFlowFieldKernel(mesh.Nelements, 
//                     T, 
//                     mesh.o_x,
//                     mesh.o_y,
//                     mesh.o_z,
//                     o_U,
//                     o_Q);

//   if(cubature){
//      advectionVolumeKernel(mesh.Nelements,
//                          mesh.o_vgeo,
//                          mesh.o_cubvgeo,
//                          mesh.o_cubDWmatrices,
//                          mesh.o_cubInterpT,
//                          mesh.o_cubProjectT,
//                          mesh.o_x,
//                          mesh.o_y,
//                          mesh.o_z,
//                          T,
//                          o_U,
//                          o_Q,
//                          o_RHS);
//   }else{
//   advectionVolumeKernel(mesh.Nelements,
//                mesh.o_vgeo,
//                mesh.o_Dmatrices,
//                T,
//                mesh.o_x,
//                mesh.o_y,
//                mesh.o_z,
//                o_U,
//                o_Q,
//                o_RHS);
//  }

//   traceHalo->ExchangeFinish(o_Q, 1, ogs_dfloat);
  
//   if(cubature){
//     advectionSurfaceKernel(mesh.Nelements,
//                             mesh.o_vgeo,
//                             mesh.o_cubsgeo,
//                             mesh.o_vmapM,
//                             mesh.o_vmapP,
//                             mesh.o_EToB,
//                             mesh.o_intInterpT,
//                             mesh.o_intLIFTT,
//                             mesh.o_intx,
//                             mesh.o_inty,
//                             mesh.o_intz,
//                             T,
//                             o_U,
//                             o_Q,
//                             o_RHS);

//   }else{
//   advectionSurfaceKernel(mesh.Nelements,
//                 mesh.o_sgeo,
//                 mesh.o_LIFTT,
//                 mesh.o_vmapM,
//                 mesh.o_vmapP,
//                 mesh.o_EToB,
//                 T,
//                 mesh.o_x,
//                 mesh.o_y,
//                 mesh.o_z,
//                 o_U,
//                 o_Q,
//                 o_RHS);
//  }
// }



// void lss_t::rhsa(occa::memory& o_Q, const dfloat time, const dfloat dt){

//   if(redistance){
//     if(timeStepper->outputStep==0){

//    historyIndex +=1; 

//   // const dfloat dt = timeStepper->GetTimeStep(); 
//   // AK: Move to Device, note that time is not increased on timeStepper yet
//   rtime[shiftIndex] = time + dt; o_rtime.copyFrom(rtime); 


//   const dlong N = mesh.Nelements*mesh.Np*Nfields; 
//   o_phiH.copyFrom(o_Q, N*sizeof(dfloat), shiftIndex*N*sizeof(dfloat),0); 

//    if(historyIndex==(Nrecon/2)){
//     // Create Initial History
//     initialHistoryKernel(mesh.Nelements, 
//                     shiftIndex, 
//                     N, 
//                     o_phiH);
//    }



//   if(historyIndex>=(Nrecon/2))
//   reconstructENOKernel(mesh.Nelements, 
//                     shiftIndex, 
//                     N, 
//                     o_rtime,
//                     o_phiH,
//                     o_phi);

//   shiftIndex = (shiftIndex+Nrecon-1)%Nrecon;

//  DetectTroubledCells(o_Q, subcell->o_ElementList); 

// }


// }


// }



// //evaluate ODE rhs = f(q,t)
// void lss_t::rhsf_subcell(occa::memory& o_Q, occa::memory &o_sQ, occa::memory& o_RHS,occa::memory& o_sRHS, const dfloat T){
  
//    Redistance_Subcell(o_Q, o_sQ, o_RHS, o_sRHS, T); 
// }


// void lss_t::Redistance(occa::memory& o_Q, occa::memory& o_RHS, const dfloat T){

//     traceHalo->ExchangeStart(o_Q, 1, ogs_dfloat);

//     redistanceVolumeKernel(mesh.Nelements,
//                mesh.o_vgeo,
//                mesh.o_DWmatrices, // o_Dmatrices
//                T,
//                mesh.o_x,
//                mesh.o_y,
//                mesh.o_z,
//                o_Q,
//                o_gradq);


//   traceHalo->ExchangeFinish(o_Q, 1, ogs_dfloat);
  
//   redistanceSurfaceKernel(mesh.Nelements,
//                 mesh.o_sgeo,
//                 mesh.o_LIFTT,
//                 mesh.o_vmapM,
//                 mesh.o_vmapP,
//                 mesh.o_EToB,
//                 T,
//                 mesh.o_x,
//                 mesh.o_y,
//                 mesh.o_z,
//                 o_Q,
//                 o_gradq,
//                 o_RHS);
//   }


// // this function reconstruct DG solution from FV cells and addd to DG storage,q
// void lss_t::reconstruct_subcell(occa::memory& o_Q, occa::memory& o_sQ){

// int all = 0; 

// for(int fieldid=0; fieldid<Nfields; fieldid++){
// reconstructKernel(mesh.Nelements, 
//                    fieldid,
//                      all,
//                     subcell->o_ElementList,
//                     subcell->o_RMT,
//                     o_sQ, 
//                     o_Q);
// }
// }

// void lss_t::Redistance_Subcell(occa::memory& o_Q, occa::memory & o_sQ, 
//                                 occa::memory& o_RHS, occa::memory& o_sRHS, const dfloat T){

//     // extract q halo on DEVICE
//     // traceHalo->ExchangeStart(o_Q, 1, ogs_dfloat);

//    // DetectTroubledCells(o_Q, subcell->o_ElementList); 

//     // mesh.halo->ExchangeStart(o_Q, mesh.Np, ogs_dfloat); 

// #if 1
// const dfloat Nlocal= mesh.Nelements*mesh.Np*Nfields; 
//  // First gather-scatter then interpolate to the vertices of subcells
//  // o_gsq.copyFrom(o_Q, Nlocal*sizeof(dfloat),0, fieldid*Nlocal*sizeof(dfloat)); 
//  // o_gsq.copyFrom(o_Q, Nlocal*sizeof(dfloat),0,0); 
//  // mesh.ogs->GatherScatterVec(o_gsq, Nfields, ogs_dfloat, ogs_add, ogs_sym);
//  // mesh.ogs->GatherScatterMany(o_gsq, Nfields, mesh.Np, ogs_dfloat, ogs_add, ogs_sym);
// #endif


//     for(int fieldid=0; fieldid<Nfields; fieldid++){


//     // Compute dq/dx for DG and DGFV elements
//      partialRedistanceVolumeKernel(mesh.Nelements,
//                                     fieldid,
//                                    subcell->o_ElementList,
//                                    offset, 
//                                    mesh.o_vgeo,
//                                    mesh.o_DWmatrices, // o_Dmatrices
//                                    T,
//                                    mesh.o_x,
//                                    mesh.o_y,
//                                    mesh.o_z,
//                                    o_Q,
//                                    o_gradq);

//      const int all = 0; // all elements->do not use element list
//      // Project DG solution to FV face values directly for DGFV type
//      projectDGKernel(mesh.Nelements,
//                     fieldid, 
//                     all, 
//                     subcell->o_ElementList, 
//                     mesh.o_vmapM,
//                     subcell->o_mFToE,
//                     subcell->o_mFToF,
//                     subcell->o_PFMT,
//                     o_Q,
//                     o_sface); 


//     // mesh.halo->ExchangeFinish(o_Q, mesh.Np, ogs_dfloat); 


//      // Add the Contribution of DG Elements
//      partialRedistanceSurfaceKernel(mesh.Nelements,
//                 fieldid, 
//                 subcell->o_ElementList,
//                 mesh.o_sgeo,
//                 mesh.o_LIFTT,
//                 mesh.o_vmapM,
//                 mesh.o_vmapP,
//                 mesh.o_EToB,
//                 T,
//                 mesh.o_x,
//                 mesh.o_y,
//                 mesh.o_z,
//                 o_Q,
//                 o_gradq,
//                 o_RHS);



//     // Compute Cell-Centered Subcell Values 
//     projectKernel(mesh.Nelements + mesh.totalHaloPairs,
//                   fieldid,
//                   all,  
//                   subcell->o_ElementList, 
//                   subcell->o_PMT,
//                   o_Q, 
//                   o_sQ);

// #if 1       
//  // Reconstruct face values for all subcells 
//      reconstructFaceKernel(mesh.Nelements, 
//                           fieldid,
//                           subcell->o_ElementList,
//                           subcell->o_vgeo,
//                           subcell->o_sgeo, 
//                           subcell->o_emapP,
//                           subcell->o_fmapP, 
//                           o_Q, 
//                           o_sQ,
//                           o_sface);   


// #else

// // reconstructFaceKernel(mesh.Nelements, 
// //                           fieldid,
// //                           subcell->o_ElementList,
// //                           subcell->o_vgeo,
// //                           subcell->o_sgeo, 
// //                           subcell->o_emapP,
// //                           subcell->o_fmapP, 
// //                           o_Q, 
// //                           o_sQ,
// //                           o_sface);   

//  extractFieldKernel(mesh.Nelements, fieldid, o_Q, o_gsq); 

//  mesh.ogs->GatherScatter(o_gsq, ogs_dfloat, ogs_add, ogs_sym);
 
//  // Now project it
//   projectVertexKernel(mesh.Nelements,
//                      fieldid,
//                      1,  
//                      subcell->o_ElementList, 
//                      o_invDegree, 
//                      subcell->o_PVMT,
//                      o_gsq, 
//                      o_svq);



//  reconstructFace2Kernel(mesh.Nelements, 
//                           fieldid,
//                           subcell->o_mEToV, 
//                           subcell->o_ElementList,
//                           subcell->o_vgeo,
//                           subcell->o_sgeo, 
//                           subcell->o_emapP,
//                           subcell->o_fmapP, 
//                           o_Q, 
//                           o_sQ,
//                           o_svq,
//                           o_sface);   

// #endif


//     //  mesh.halo->Exchange(o_sface, subcell->Nsubcells*mesh.Nfaces, ogs_dfloat);

//     // FV compute 
//      subcellComputeKernel(mesh.Nelements, 
//                               fieldid, 
//                               subcell->o_ElementList,
//                               subcell->o_emapP,
//                               subcell->o_fmapP, 
//                               subcell->o_RMT,
//                               subcell->o_vgeo,
//                               subcell->o_sgeo, 
//                               o_Q,
//                               o_sface, 
//                               o_sRHS);  

//       // mixedRedistanceSurfaceKernel(mesh.Nelements,
//       //                     fieldid, 
//       //                     subcell->o_ElementList,
//       //                     subcell->o_EToE, 
//       //                     mesh.o_EToB,
//       //                     mesh.o_sgeo,
//       //                     subcell->o_sgeo, 
//       //                     subcell->o_vgeo,
//       //                     mesh.o_LIFTT,
//       //                     subcell->o_SLIFTT,
//       //                     subcell->o_RFMT,
//       //                     mesh.o_vmapM,
//       //                     mesh.o_vmapP,
//       //                     subcell->o_emapP,
//       //                     subcell->o_fmapP, 
//       //                     subcell->o_mFToE,
//       //                     subcell->o_mFToF,
//       //                     subcell->o_mDGID,
//       //                     T,
//       //                     mesh.o_x,
//       //                     mesh.o_y,
//       //                     mesh.o_z,
//       //                     o_Q,
//       //                     o_sface,
//       //                     o_gradq,
//       //                     o_RHS);

//     }
// }


// void lss_t::DetectTroubledCells(occa::memory& o_Q, occa::memory& o_ElementList){

// // Currently modal decay + skyline only
// if(subcellStabilization){

//   for (int i=0; i<Nfields; i++){
//     int fieldid  = i; 

//   if(settings.compareSetting("INDICATOR TYPE", "MDA")){
//   indicatorMDAKernel(mesh.Nelements,
//                 fieldid, 
//                 mesh.o_ggeo,
//                 subcell->o_ModMap, 
//                 mesh.o_MM,
//                 subcell->o_invVT,
//                 subcell->o_LSF,
//                 subcell->o_BLD, 
//                 o_Q,
//                 o_ElementList); 


//   }else if(settings.compareSetting("INDICATOR TYPE", "MDH")){
//   indicatorMDHKernel(mesh.Nelements,
//                 fieldid, 
//                 mesh.o_ggeo,
//                 mesh.o_MM,
//                 subcell->o_interpNM1T,
//                 o_Q,
//                 o_ElementList); 
//   }else{
//     printf("The detector type is not impmented\n");
//   }
// }

// mesh.halo->Exchange(o_ElementList, Nfields, ogs_dlong);

// for (int i=0; i<Nfields; i++){
//     int fieldid  = i;   
//   // Find neighs i.e. DG neighs of FV elements....
//   findNeighKernel(mesh.Nelements,
//                 fieldid,
//                  subcell->o_EToE, 
//                  o_ElementList);   
//   }  
// }

// }
