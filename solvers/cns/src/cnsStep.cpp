/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "cns.hpp"

dfloat cns_t::MaxWaveSpeed(deviceMemory<dfloat>& o_Q, const dfloat T){

  deviceMemory<dfloat> o_maxSpeed = platform.reserve<dfloat>(mesh.Nelements);

  maxWaveSpeedKernel(mesh.Nelements,
                     BCStateID,
                     mesh.o_vgeo,
                     mesh.o_sgeo,
                     mesh.o_vmapM,
                     mesh.o_EToB,
                     o_pCoeff, 
                     o_flowStates,
                     T,
                     mesh.o_x,
                     mesh.o_y,
                     mesh.o_z,
                     o_Q,
                     o_maxSpeed);

  const dfloat vmax = platform.linAlg().max(mesh.Nelements, o_maxSpeed, mesh.comm);
  return vmax;
}


//evaluate ODE rhs = f(q,t)
void cns_t::rhsf(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
  
  if(EulerSolve){
    if(stabType==Stab::ARTDIFF){rhsEulerArtDiff(o_Q, o_RHS, T);}
  }else{
    if(stabType==Stab::ARTDIFF){ rhsArtDiff(o_Q, o_RHS, T);}
    else if(stabType==Stab::LIMITER){ rhsLimiter(o_Q, o_RHS, T);}
    else if(stabType==Stab::NOSTAB ){ rhsNoStab(o_Q, o_RHS, T); }
  }
}


//evaluate ODE rhs = f(q,t)
void cns_t::rhsEulerArtDiff(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
  dlong NlocalGrads = mesh.Nelements*mesh.Np*Ngrads;
  dlong NhaloGrads  = mesh.totalHaloPairs*mesh.Np*Ngrads;
  deviceMemory<dfloat> o_gradq = platform.reserve<dfloat>(NlocalGrads+NhaloGrads);  

  // extract q trace halo and start exchange
  fieldTraceHalo.ExchangeStart(o_Q, 1);

  // compute volume contributions to gradients
  gradVolumeKernel(mesh.Nelements,
                 mesh.o_vgeo,
                 mesh.o_D,
                 o_Q,
                 o_gradq);    
  
  // complete trace halo exchange
  fieldTraceHalo.ExchangeFinish(o_Q, 1);

  // compute surface contributions to gradients
  gradSurfaceKernel(mesh.Nelements,
                    BCStateID, 
                    mesh.o_sgeo,
                    mesh.o_LIFT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    o_pCoeff, 
                    o_flowStates, 
                    T,
                    o_Q,
                    o_gradq);

  // extract viscousStresses trace halo and start exchange
  gradTraceHalo.ExchangeStart(o_gradq, 1);


  // const dfloat vmax = 0.0; 
   // platform.linAlg().scale(mesh.Nelements*mesh.Np, vmax, o_viscosity); 

  applyStab(o_Q, o_gradq, T);

  // compute volume contribution to cns RHS
  if (cubature) {
     cubatureVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubD,
                         mesh.o_cubPDT,
                         mesh.o_cubInterp,
                         mesh.o_cubProject,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         o_pCoeff, 
                         T,
                         o_viscosity,
                         o_Q,
                         o_gradq,
                         o_RHS);
  } else {
    volumeKernel(mesh.Nelements,
                 mesh.o_vgeo,
                 mesh.o_D,
                 mesh.o_x,
                 mesh.o_y,
                 mesh.o_z,
                 o_pCoeff, 
                 T,
                 o_viscosity,
                 o_Q,
                 o_gradq,
                 o_RHS);
  }

  // complete trace halo exchange
  gradTraceHalo.ExchangeFinish(o_gradq, 1);

  if (cubature) {
     // dfloat hmin = mesh.MinCharacteristicLength();
     cubatureSurfaceKernel(mesh.Nelements,
                            BCStateID, 
                            mesh.o_vgeo,
                            mesh.o_cubsgeo,
                            mesh.o_vmapM,
                            mesh.o_vmapP,
                            mesh.o_EToB,
                            mesh.o_intInterp,
                            mesh.o_intLIFT,
                            mesh.o_intx,
                            mesh.o_inty,
                            mesh.o_intz,
                            o_pCoeff, 
                            o_flowStates, 
                            T,
                            o_viscosity,
                            o_Q,
                            o_gradq,
                            o_RHS);
    } else {
      // dfloat hmin = mesh.MinCharacteristicLength();
     surfaceKernel(mesh.Nelements,
                        mesh.o_sgeo,
                        mesh.o_LIFT,
                        mesh.o_vmapM,
                        mesh.o_vmapP,
                        mesh.o_EToB,
                        mesh.o_x,
                        mesh.o_y,
                        mesh.o_z,
                        o_pCoeff, 
                        T,
                        o_viscosity,
                        o_Q,
                        o_gradq,
                        o_RHS);
  }
}





















// //evaluate ODE rhs = f(q,t)
// void cns_t::rhsArtDiff(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
//   dlong NlocalGrads = mesh.Nelements*mesh.Np*Ngrads;
//   dlong NhaloGrads  = mesh.totalHaloPairs*mesh.Np*Ngrads;
//   deviceMemory<dfloat> o_gradq = platform.reserve<dfloat>(NlocalGrads+NhaloGrads);  

//   // extract q trace halo and start exchange
//   fieldTraceHalo.ExchangeStart(o_Q, 1);
  
//   // compute volume contributions to gradients
//   gradVolumeKernel(mesh.Nelements,
//                  mesh.o_vgeo,
//                  mesh.o_D,
//                  o_Q,
//                  o_gradq);    

//   // complete trace halo exchange
//   fieldTraceHalo.ExchangeFinish(o_Q, 1);

//   // compute surface contributions to gradients
//   gradSurfaceKernel(mesh.Nelements,
//                     BCStateID, 
//                     mesh.o_sgeo,
//                     mesh.o_LIFT,
//                     mesh.o_vmapM,
//                     mesh.o_vmapP,
//                     mesh.o_EToB,
//                     mesh.o_x,
//                     mesh.o_y,
//                     mesh.o_z,
//                     o_pCoeff, 
//                     o_flowStates, 
//                     T,
//                     o_Q,
//                     o_gradq);


//   // Apply Stabilization: assumes gradQ is available but not communicated yet!
//   applyStab(o_Q, o_gradq, T);
   
//   // extract viscousStresses trace halo and start exchange
//   gradTraceHalo.ExchangeStart(o_gradq, 1);

//   // compute volume contribution to cns RHS
//   if (cubature) {
//      cubatureVolumeKernel(mesh.Nelements,
//                          mesh.o_vgeo,
//                          mesh.o_cubvgeo,
//                          mesh.o_cubD,
//                          mesh.o_cubPDT,
//                          mesh.o_cubInterp,
//                          mesh.o_cubProject,
//                          mesh.o_x,
//                          mesh.o_y,
//                          mesh.o_z,
//                          o_pCoeff, 
//                          T,
//                          o_viscosity,
//                          o_Q,
//                          o_gradq,
//                          o_RHS);
//   } else {
//     volumeKernel(mesh.Nelements,
//                  mesh.o_vgeo,
//                  mesh.o_D,
//                  mesh.o_x,
//                  mesh.o_y,
//                  mesh.o_z,
//                  o_pCoeff, 
//                  T,
//                  o_viscosity,
//                  o_Q,
//                  o_gradq,
//                  o_RHS);
//   }

//   // complete trace halo exchange
//   gradTraceHalo.ExchangeFinish(o_gradq, 1);

//   if (cubature) {
//      cubatureSurfaceKernel(mesh.Nelements,
//                             BCStateID, 
//                             mesh.o_vgeo,
//                             mesh.o_cubsgeo,
//                             mesh.o_vmapM,
//                             mesh.o_vmapP,
//                             mesh.o_EToB,
//                             mesh.o_intInterp,
//                             mesh.o_intLIFT,
//                             mesh.o_intx,
//                             mesh.o_inty,
//                             mesh.o_intz,
//                             o_pCoeff, 
//                             o_flowStates, 
//                             T,
//                             o_viscosity,
//                             o_Q,
//                             o_gradq,
//                             o_RHS);
//     } else {
//      surfaceKernel(mesh.Nelements,
//                         mesh.o_sgeo,
//                         mesh.o_LIFT,
//                         mesh.o_vmapM,
//                         mesh.o_vmapP,
//                         mesh.o_EToB,
//                         mesh.o_x,
//                         mesh.o_y,
//                         mesh.o_z,
//                         o_pCoeff, 
//                         T,
//                         o_viscosity,
//                         o_Q,
//                         o_gradq,
//                         o_RHS);
//   }
// }



//evaluate ODE rhs = f(q,t)
void cns_t::rhsArtDiff(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
  dlong NlocalGrads = mesh.Nelements*mesh.Np*Ngrads;
  dlong NhaloGrads  = mesh.totalHaloPairs*mesh.Np*Ngrads;
  deviceMemory<dfloat> o_gradq = platform.reserve<dfloat>(NlocalGrads+NhaloGrads);  

  // extract q trace halo and start exchange
  fieldTraceHalo.ExchangeStart(o_Q, 1);

  // if (cubature) {
  //   cubatureGradVolumeKernel(mesh.Nelements,
  //                    mesh.o_vgeo,
  //                    mesh.o_cubvgeo,
  //                    mesh.o_cubD,
  //                    mesh.o_cubPDT,
  //                    mesh.o_cubInterp,
  //                    mesh.o_cubProject,
  //                    o_Q,
  //                    o_gradq);    
  // }else{
    // compute volume contributions to gradients
    gradVolumeKernel(mesh.Nelements,
                   mesh.o_vgeo,
                   mesh.o_D,
                   o_Q,
                   o_gradq);    
  // }
  // complete trace halo exchange
  fieldTraceHalo.ExchangeFinish(o_Q, 1);


  // if (cubature) {
  //   cubatureGradSurfaceKernel(mesh.Nelements,
  //                           BCStateID, 
  //                           mesh.o_vgeo,
  //                           mesh.o_cubsgeo,
  //                           mesh.o_vmapM,
  //                           mesh.o_vmapP,
  //                           mesh.o_EToB,
  //                           mesh.o_intInterp,
  //                           mesh.o_intLIFT,
  //                           mesh.o_intx,
  //                           mesh.o_inty,
  //                           mesh.o_intz,
  //                           o_pCoeff, 
  //                           o_flowStates, 
  //                           T,
  //                           o_Q,
  //                           o_gradq);
  // }else{
   // compute surface contributions to gradients
  gradSurfaceKernel(mesh.Nelements,
                    BCStateID, 
                    mesh.o_sgeo,
                    mesh.o_LIFT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    o_pCoeff, 
                    o_flowStates, 
                    T,
                    o_Q,
                    o_gradq);
  // }

  // extract viscousStresses trace halo and start exchange
  gradTraceHalo.ExchangeStart(o_gradq, 1);


  applyStab(o_Q, o_gradq, T);
   // const dfloat vmax = 0.0; 
   // platform.linAlg().scale(mesh.Nelements*mesh.Np, vmax, o_viscosity); 


  // compute volume contribution to cns RHS
  if (cubature) {
     cubatureVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubD,
                         mesh.o_cubPDT,
                         mesh.o_cubInterp,
                         mesh.o_cubProject,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         o_pCoeff, 
                         T,
                         o_viscosity,
                         o_Q,
                         o_gradq,
                         o_RHS);
  } else {
    volumeKernel(mesh.Nelements,
                 mesh.o_vgeo,
                 mesh.o_D,
                 mesh.o_x,
                 mesh.o_y,
                 mesh.o_z,
                 o_pCoeff, 
                 T,
                 o_viscosity,
                 o_Q,
                 o_gradq,
                 o_RHS);
  }

  // complete trace halo exchange
  gradTraceHalo.ExchangeFinish(o_gradq, 1);

  if (cubature) {
     // dfloat hmin = mesh.MinCharacteristicLength();
     cubatureSurfaceKernel(mesh.Nelements,
                            BCStateID, 
                            mesh.o_vgeo,
                            mesh.o_cubsgeo,
                            mesh.o_vmapM,
                            mesh.o_vmapP,
                            mesh.o_EToB,
                            mesh.o_intInterp,
                            mesh.o_intLIFT,
                            mesh.o_intx,
                            mesh.o_inty,
                            mesh.o_intz,
                            o_pCoeff, 
                            o_flowStates, 
                            T,
                            o_viscosity,
                            o_Q,
                            o_gradq,
                            o_RHS);
    } else {
      // dfloat hmin = mesh.MinCharacteristicLength();
     surfaceKernel(mesh.Nelements,
                        mesh.o_sgeo,
                        mesh.o_LIFT,
                        mesh.o_vmapM,
                        mesh.o_vmapP,
                        mesh.o_EToB,
                        mesh.o_x,
                        mesh.o_y,
                        mesh.o_z,
                        o_pCoeff, 
                        T,
                        o_viscosity,
                        o_Q,
                        o_gradq,
                        o_RHS);
  }
}





// //evaluate ODE rhs = f(q,t)
// void cns_t::rhsArtDiff(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
//   dlong NlocalGrads = mesh.Nelements*mesh.Np*Ngrads;
//   dlong NhaloGrads  = mesh.totalHaloPairs*mesh.Np*Ngrads;
//   deviceMemory<dfloat> o_gradq = platform.reserve<dfloat>(NlocalGrads+NhaloGrads);  

//   // extract q trace halo and start exchange
//   fieldTraceHalo.ExchangeStart(o_Q, 1);

//   if (cubature) {
//     cubatureGradVolumeKernel(mesh.Nelements,
//                      mesh.o_vgeo,
//                      mesh.o_cubvgeo,
//                      mesh.o_cubD,
//                      mesh.o_cubPDT,
//                      mesh.o_cubInterp,
//                      mesh.o_cubProject,
//                      o_Q,
//                      o_gradq);    
//   }else{
//     // compute volume contributions to gradients
//     gradVolumeKernel(mesh.Nelements,
//                    mesh.o_vgeo,
//                    mesh.o_D,
//                    o_Q,
//                    o_gradq);    
//   }
//   // complete trace halo exchange
//   fieldTraceHalo.ExchangeFinish(o_Q, 1);




//   stab.Apply(o_Q, o_gradq, T); 

//   #if 0
//   // dfloat vmax = MaxWaveSpeed(o_Q, T);
//   dfloat vmax = 0.0;
//   // printf("vmax = %.f\n", vmax);
//   platform.linAlg().scale(mesh.Nelements*mesh.Nverts, vmax, stab.o_viscosity); 
//   #endif


//   if (cubature) {
//     cubatureGradSurfaceKernel(mesh.Nelements,
//                             mesh.o_vgeo,
//                             mesh.o_cubsgeo,
//                             mesh.o_vmapM,
//                             mesh.o_vmapP,
//                             mesh.o_EToB,
//                             mesh.o_intInterp,
//                             mesh.o_intLIFT,
//                             mesh.o_intx,
//                             mesh.o_inty,
//                             mesh.o_intz,
//                             o_pCoeff, 
//                             T,
//                             o_Q,
//                             o_gradq);
//   }else{
//    // compute surface contributions to gradients
//   gradSurfaceKernel(mesh.Nelements,
//                     mesh.o_sgeo,
//                     mesh.o_LIFT,
//                     mesh.o_vmapM,
//                     mesh.o_vmapP,
//                     mesh.o_EToB,
//                     mesh.o_x,
//                     mesh.o_y,
//                     mesh.o_z,
//                     o_pCoeff, 
//                     T,
//                     o_Q,
//                     o_gradq);
//   }



//   // extract viscousStresses trace halo and start exchange
//   gradTraceHalo.ExchangeStart(o_gradq, 1);

//   // compute volume contribution to cns RHS
//   if (cubature) {
//      cubatureVolumeKernel(mesh.Nelements,
//                          mesh.o_vgeo,
//                          mesh.o_cubvgeo,
//                          mesh.o_cubD,
//                          mesh.o_cubPDT,
//                          mesh.o_cubInterp,
//                          mesh.o_cubProject,
//                          mesh.o_x,
//                          mesh.o_y,
//                          mesh.o_z,
//                          o_pCoeff, 
//                          T,
//                          stab.o_viscosity,
//                          o_Q,
//                          o_gradq,
//                          o_RHS);
//   } else {
//     volumeKernel(mesh.Nelements,
//                  mesh.o_vgeo,
//                  mesh.o_D,
//                  mesh.o_x,
//                  mesh.o_y,
//                  mesh.o_z,
//                  o_pCoeff, 
//                  T,
//                  stab.o_viscosity,
//                  o_Q,
//                  o_gradq,
//                  o_RHS);
//   }

//   // complete trace halo exchange
//   gradTraceHalo.ExchangeFinish(o_gradq, 1);

//   if (cubature) {
//      // dfloat hmin = mesh.MinCharacteristicLength();
//      cubatureSurfaceKernel(mesh.Nelements,
//                             mesh.o_vgeo,
//                             mesh.o_cubsgeo,
//                             mesh.o_vmapM,
//                             mesh.o_vmapP,
//                             mesh.o_EToB,
//                             mesh.o_intInterp,
//                             mesh.o_intLIFT,
//                             mesh.o_intx,
//                             mesh.o_inty,
//                             mesh.o_intz,
//                             o_pCoeff, 
//                             T,
//                             stab.o_viscosity,
//                             o_Q,
//                             o_gradq,
//                             o_RHS);
//     } else {
//       // dfloat hmin = mesh.MinCharacteristicLength();
//      surfaceKernel(mesh.Nelements,
//                         mesh.o_sgeo,
//                         mesh.o_LIFT,
//                         mesh.o_vmapM,
//                         mesh.o_vmapP,
//                         mesh.o_EToB,
//                         mesh.o_x,
//                         mesh.o_y,
//                         mesh.o_z,
//                         o_pCoeff, 
//                         T,
//                         stab.o_viscosity,
//                         o_Q,
//                         o_gradq,
//                         o_RHS);
//   }
// }


//evaluate ODE rhs = f(q,t)
void cns_t::rhsLimiter(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

//   dlong NlocalGrads = mesh.Nelements*mesh.Np*Ngrads;
//   dlong NhaloGrads  = mesh.totalHaloPairs*mesh.Np*Ngrads;
//   deviceMemory<dfloat> o_gradq = platform.reserve<dfloat>(NlocalGrads+NhaloGrads);

// #if 1
// {
//   stab.Apply(o_Q, o_RHS, T);

//   // int alpha = 1; 
//   // platform.linAlg().set(mesh.Nelements, alpha, stab.o_elementList);
  
//   // Project Solution To Cell-Mean Values
//   limiterReconstructKernel(mesh.Nelements,
//                       mesh.o_vgeo,
//                       mesh.o_sgeo,
//                       stab.o_vgeo, 
//                       mesh.o_EToB, 
//                       mesh.o_vertexNodes, 
//                       stab.o_elementList, 
//                       mesh.o_vmapM, 
//                       mesh.o_vmapP, 
//                       mesh.o_x, 
//                       mesh.o_y, 
//                       mesh.o_z,
//                       o_pCoeff,
//                       T,
//                       o_Q, 
//                       stab.o_qc, 
//                       stab.o_qv, 
//                       stab.o_dq);

//   limiterGradientKernel(mesh.Nelements,
//                             mesh.o_vgeo,
//                             mesh.o_sgeo,
//                             stab.o_vgeo, 
//                             mesh.o_EToB, 
//                             mesh.o_vertexNodes, 
//                             stab.o_elementList, 
//                             mesh.o_vmapM, 
//                             mesh.o_vmapP, 
//                             mesh.o_x, 
//                             mesh.o_y, 
//                             mesh.o_z,
//                             o_pCoeff,
//                             T,
//                             o_Q, 
//                             stab.o_qc, 
//                             stab.o_qv, 
//                             stab.o_dq);
//  }
// #endif



//   // extract q trace halo and start exchange
//   fieldTraceHalo.ExchangeStart(o_Q, 1);

//   if (cubature) {
//     cubatureGradVolumeKernel(mesh.Nelements,
//                      mesh.o_vgeo,
//                      mesh.o_cubvgeo,
//                      mesh.o_cubD,
//                      mesh.o_cubPDT,
//                      mesh.o_cubInterp,
//                      mesh.o_cubProject,
//                      o_Q,
//                      o_gradq);    
//   }else{
//     // compute volume contributions to gradients
//     gradVolumeKernel(mesh.Nelements,
//                    mesh.o_vgeo,
//                    mesh.o_D,
//                    o_Q,
//                    o_gradq);    
//   }



//   // complete trace halo exchange
//   fieldTraceHalo.ExchangeFinish(o_Q, 1);


//   if (cubature) {
//     cubatureGradSurfaceKernel(mesh.Nelements,
//                               mesh.o_vgeo,
//                               mesh.o_cubsgeo,
//                               mesh.o_vmapM,
//                               mesh.o_vmapP,
//                               mesh.o_EToB,
//                               mesh.o_intInterp,
//                               mesh.o_intLIFT,
//                               mesh.o_intx,
//                               mesh.o_inty,
//                               mesh.o_intz,
//                               o_pCoeff, 
//                               T,
//                               o_Q,
//                               o_gradq);
//   }else{
//    // compute surface contributions to gradients
//   gradSurfaceKernel(mesh.Nelements,
//                     mesh.o_sgeo,
//                     mesh.o_LIFT,
//                     mesh.o_vmapM,
//                     mesh.o_vmapP,
//                     mesh.o_EToB,
//                     mesh.o_x,
//                     mesh.o_y,
//                     mesh.o_z,
//                     o_pCoeff, 
//                     T,
//                     o_Q,
//                     o_gradq);
//   }



//   // extract viscousStresses trace halo and start exchange
//   gradTraceHalo.ExchangeStart(o_gradq, 1);

//   // compute volume contribution to cns RHS
//   if (cubature) {
//      cubatureVolumeKernel(mesh.Nelements,
//                          mesh.o_vgeo,
//                          mesh.o_cubvgeo,
//                          mesh.o_cubD,
//                          mesh.o_cubPDT,
//                          mesh.o_cubInterp,
//                          mesh.o_cubProject,
//                          mesh.o_x,
//                          mesh.o_y,
//                          mesh.o_z,
//                          o_pCoeff, 
//                          T,
//                          o_Q,
//                          o_gradq,
//                          o_RHS);
//   } else {
//     volumeKernel(mesh.Nelements,
//                  mesh.o_vgeo,
//                  mesh.o_D,
//                  mesh.o_x,
//                  mesh.o_y,
//                  mesh.o_z,
//                  o_pCoeff, 
//                  T,
//                  o_Q,
//                  o_gradq,
//                  o_RHS);
//   }

//   // complete trace halo exchange
//   gradTraceHalo.ExchangeFinish(o_gradq, 1);

//   if (cubature) {
//      // dfloat hmin = mesh.MinCharacteristicLength();
//      cubatureSurfaceKernel(mesh.Nelements,
//                             mesh.o_vgeo,
//                             mesh.o_cubsgeo,
//                             mesh.o_vmapM,
//                             mesh.o_vmapP,
//                             mesh.o_EToB,
//                             mesh.o_intInterp,
//                             mesh.o_intLIFT,
//                             mesh.o_intx,
//                             mesh.o_inty,
//                             mesh.o_intz,
//                             o_pCoeff, 
//                             T,
//                             o_Q,
//                             o_gradq,
//                             o_RHS);
//     } else {
//       // dfloat hmin = mesh.MinCharacteristicLength();
//      surfaceKernel(mesh.Nelements,
//                         mesh.o_sgeo,
//                         mesh.o_LIFT,
//                         mesh.o_vmapM,
//                         mesh.o_vmapP,
//                         mesh.o_EToB,
//                         mesh.o_x,
//                         mesh.o_y,
//                         mesh.o_z,
//                         o_pCoeff, 
//                         T,
//                         o_Q,
//                         o_gradq,
//                         o_RHS);
//   }
}

//evaluate ODE rhs = f(q,t)
void cns_t::rhsNoStab(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

  dlong NlocalGrads = mesh.Nelements*mesh.Np*Ngrads;
  dlong NhaloGrads  = mesh.totalHaloPairs*mesh.Np*Ngrads;
  deviceMemory<dfloat> o_gradq = platform.reserve<dfloat>(NlocalGrads+NhaloGrads);

  // extract q trace halo and start exchange
  fieldTraceHalo.ExchangeStart(o_Q, 1);

 // compute volume contributions to gradients
  gradVolumeKernel(mesh.Nelements,
                   mesh.o_vgeo,
                   mesh.o_D,
                   o_Q,
                   o_gradq);

  // complete trace halo exchange
  fieldTraceHalo.ExchangeFinish(o_Q, 1);

  // compute surface contributions to gradients
  gradSurfaceKernel(mesh.Nelements,
                    mesh.o_sgeo,
                    mesh.o_LIFT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    o_pCoeff, 
                    T,
                    o_Q,
                    o_gradq);

  // extract viscousStresses trace halo and start exchange
  gradTraceHalo.ExchangeStart(o_gradq, 1);

  // compute volume contribution to cns RHS
  if (cubature) {
     cubatureVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubD,
                         mesh.o_cubPDT,
                         mesh.o_cubInterp,
                         mesh.o_cubProject,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         o_pCoeff, 
                         T,
                         o_Q,
                         o_gradq,
                         o_RHS);
  } else {
    volumeKernel(mesh.Nelements,
                 mesh.o_vgeo,
                 mesh.o_D,
                 mesh.o_x,
                 mesh.o_y,
                 mesh.o_z,
                 o_pCoeff, 
                 T,
                 o_Q,
                 o_gradq,
                 o_RHS);
  }

  // complete trace halo exchange
  gradTraceHalo.ExchangeFinish(o_gradq, 1);

  if (cubature) {
      cubatureSurfaceKernel(mesh.Nelements,
                            mesh.o_vgeo,
                            mesh.o_cubsgeo,
                            mesh.o_vmapM,
                            mesh.o_vmapP,
                            mesh.o_EToB,
                            mesh.o_intInterp,
                            mesh.o_intLIFT,
                            mesh.o_intx,
                            mesh.o_inty,
                            mesh.o_intz,
                            o_pCoeff, 
                            T,
                            o_Q,
                            o_gradq,
                            o_RHS);
    } else {
      surfaceKernel(mesh.Nelements,
                    mesh.o_sgeo,
                    mesh.o_LIFT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    o_pCoeff, 
                    T,
                    o_Q,
                    o_gradq,
                    o_RHS);
    }
}
