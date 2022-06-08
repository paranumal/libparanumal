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

#include "lbs.hpp"

void lbs_t::Run(){

  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  initialConditionKernel(mesh.Nelements,
                         nu,
                         startTime,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         o_U);

  
  dfloat cfl=1.0;
  settings.getSetting("CFL NUMBER", cfl);

  // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat vmax = MaxWaveSpeed();
  dfloat dtAdv  = hmin/(vmax*(mesh.N+1.)*(mesh.N+1.));
  // dfloat dtVisc = 1.0/tauInv;
  // dfloat dt = (semiAnalytic) ? cfl*dtAdv : cfl*mymin(dtAdv, dtVisc);

  dfloat dt = cfl*dtAdv;


  const dfloat gamma = alpha/dt; 
  phaseFieldKernel(mesh.Nelements,
                   startTime,
                   dt, 
                   gamma,
                   nu,
                   o_LBM,
                   mesh.o_x,
                   mesh.o_y,
                   mesh.o_z,
                   o_U,
                   o_q);

  /*
    Artificial warping of time step size for multirate testing
  */
#if 0
  if (settings.compareSetting("TIME INTEGRATOR","MRAB3") ||
      settings.compareSetting("TIME INTEGRATOR","MRSAAB3"))
    dt /= (1<<(mesh.mrNlevels-1));
#endif
  timeStepper.SetTimeStep(dt);

  timeStepper.Run(*this, o_q, startTime, finalTime);


  // output norm of final solution
  {
    // Compute velocity and density
    momentsKernel(mesh.Nelements, o_LBM, o_q, o_U); 
    //compute q.M*q
    mesh.MassMatrixApply(o_U, o_Mq);

    dlong Nentries = mesh.Nelements*mesh.Np*Nmacro;
    dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);
  }
}


