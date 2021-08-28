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

#include "lbs.hpp"

void lbs_t::Report(dfloat time, int tstep){
  static int frame=0;
  // Compute velocity and density
  momentsKernel(mesh.Nelements, o_LBM, o_q, o_U); 

  //compute vorticity
  vorticityKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D, o_U, o_Vort);

  //compute q.M*q
  mesh.MassMatrixApply(o_U, o_Mq);

  dlong Nentries = mesh.Nelements*mesh.Np*Nmacro;
  dfloat norm2 = sqrt(platform.linAlg.innerProd(Nentries, o_q, o_Mq, mesh.comm));

  if(mesh.rank==0)
    printf("%5.2f (%d), %5.4f (time, timestep, norm)\n", time, tstep, norm2);

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

    // copy data back to host
    // o_q.copyTo(q);
    o_U.copyTo(U);
    o_Vort.copyTo(Vort);

    // output field files
    string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);

    // PlotFields(o_q, Vort, fname);
    PlotFields(U, Vort, fname);
  }
}
