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

#include "ins.hpp"

void ins_t::Run(){

  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  initialConditionKernel(mesh.Nelements,
                         startTime,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         nu,
                         o_u,
                         o_p);

  dfloat cfl=1.0;
  settings.getSetting("CFL NUMBER", cfl);

  // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat vmax = MaxWaveSpeed(o_u, startTime);

  dfloat dtAdvc = cfl/(vmax*(mesh.N+1.)*(mesh.N+1.));
  dfloat dtDiff = nu>0.0 ? cfl*pow(hmin, 2)/(pow(mesh.N+1,4)*nu) : 1.0e9;

  dfloat dt = 0.0;
  if (settings.compareSetting("TIME INTEGRATOR","EXTBDF3")) {
    dt = dtAdvc;
  } else if (settings.compareSetting("TIME INTEGRATOR","SSBDF3")) {
    dt = Nsubcycles*dtAdvc;
    subStepper.SetTimeStep(dtAdvc);
  } else {
    dt = std::min(dtAdvc, dtDiff);
  }

  timeStepper.SetTimeStep(dt);

  timeStepper.Run(*this, o_u, startTime, finalTime);

  // output norm of final solution
  {
    //compute U.M*U
    mesh.MassMatrixApply(o_u, o_MU);

    dlong Nentries = mesh.Nelements*mesh.Np*NVfields;
    dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_u, o_MU, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);
  }

}
