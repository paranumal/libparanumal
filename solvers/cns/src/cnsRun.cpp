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

void cns_t::Run(){

  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  initialConditionKernel(mesh.Nelements,
                         startTime,
                         o_pCoeff,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         o_q);

  dfloat cfl=1.0;
  settings.getSetting("CFL NUMBER", cfl);

  // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat vmax = MaxWaveSpeed(o_q, startTime);
  dfloat dtAdv  = cfl/(vmax*(mesh.N+1.)*(mesh.N+1.));

  dfloat visc   = 0.0; 

  if(stabType==Stab::NOSTAB){
    visc = pCoeff[MUID]; 
  }else if(stabType==Stab::ARTDIFF){
   visc = std::max(hmin/mesh.N, pCoeff[MUID]); 
  }

  dfloat dtVisc = mu > 1E-12 ? cfl*pow(hmin, 2)/(pow(mesh.N+1,4)*visc):1E12;
  dfloat dt = std::min(dtAdv, dtVisc);
  timeStepper.SetTimeStep(dt);

  if(mesh.rank==0)
  printf("time step size: %.4e startTime= %.4e FinalTime= %.4e \n", dt, startTime, finalTime);

  
   // Report(0,0); 
   // stab.Report(0.0, 0);

#if 0
   if(stabType!=Stab::NOSTAB){

    applyStab(o_q, o_q, 0);
    Report(0,0); 
      // stab.Apply(o_q, o_q, 0.0);

    // int alpha = 1; 
    // platform.linAlg().set(mesh.Nelements, alpha, stab.o_elementList);

  // for(int i=0; i<1; i++){

  // // Project Solution To Cell-Mean Values
  // limiterVertexBoundaryKernel(mesh.Nelements,
  //                     mesh.o_vgeo,
  //                     mesh.o_sgeo,
  //                     stab.o_vgeo, 
  //                     mesh.o_EToB, 
  //                     mesh.o_vertexNodes, 
  //                     stab.o_elementList, 
  //                     mesh.o_vmapM, 
  //                     mesh.o_vmapP, 
  //                     mesh.o_x, 
  //                     mesh.o_y, 
  //                     mesh.o_z,
  //                     o_pCoeff,
  //                     0.0,
  //                     o_q, 
  //                     stab.o_qc, 
  //                     stab.o_qv, 
  //                     stab.o_dq);

  // // Project Solution To Cell-Mean Values
  // limiterReconstructKernel(mesh.Nelements,
  //                     mesh.o_vgeo,
  //                     mesh.o_sgeo,
  //                     stab.o_vgeo, 
  //                     mesh.o_EToB, 
  //                     mesh.o_vertexNodes, 
  //                     stab.o_elementList, 
  //                     mesh.o_vmapM, 
  //                     mesh.o_vmapP, 
  //                     mesh.o_x, 
  //                     mesh.o_y, 
  //                     mesh.o_z,
  //                     o_pCoeff,
  //                     0.0,
  //                     o_q, 
  //                     stab.o_qc, 
  //                     stab.o_qv, 
  //                     stab.o_dq);

  // limiterGradientKernel(mesh.Nelements,
  //                           mesh.o_vgeo,
  //                           mesh.o_sgeo,
  //                           stab.o_vgeo, 
  //                           mesh.o_EToB, 
  //                           mesh.o_vertexNodes, 
  //                           stab.o_elementList, 
  //                           mesh.o_vmapM, 
  //                           mesh.o_vmapP, 
  //                           mesh.o_x, 
  //                           mesh.o_y, 
  //                           mesh.o_z,
  //                           o_pCoeff,
  //                           0.0,
  //                           o_q, 
  //                           stab.o_qc, 
  //                           stab.o_qv, 
  //                           stab.o_dq);

      // stab.Report(0.0, 0);
     // }
    }

#else
  timeStepper.Run(*this, o_q, startTime, finalTime);
// Report(0,0); 
#endif

  // output norm of final solution
  {
    //compute q.M*q
    dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
    deviceMemory<dfloat> o_Mq = platform.reserve<dfloat>(Nentries);
    mesh.MassMatrixApply(o_q, o_Mq);

    dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);
  }

}
