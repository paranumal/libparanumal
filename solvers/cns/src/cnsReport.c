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

#include "cns.h"

void cnsReport(cns_t *cns, dfloat time, setupAide &options){

  mesh3D *mesh = cns->mesh;

  cns->vorticityKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_Dmatrices,
                       cns->o_q,
                       cns->o_Vort);

  // copy data back to host
  cns->o_q.copyTo(cns->q);
  cns->o_Vort.copyTo(cns->Vort);

  // do error stuff on host
  cnsError(cns, time);

  //  cnsForces(cns, time);

#ifdef RENDER
  if(options.compareArgs("OUTPUT FILE FORMAT","PPM")){

    // copy data back to host
    cns->o_q.copyTo(cns->q);
    cns->o_Vort.copyTo(cns->Vort);
   
    //
    char fname[BUFSIZ];
    string outName;
    options.getArgs("OUTPUT FILE NAME", outName);
    cnsRenderQuad3D(cns, (char*)outName.c_str(), cns->frame++);
  }
#endif

  if(options.compareArgs("OUTPUT FILE FORMAT","VTU")){  
    // output field files
    char fname[BUFSIZ];
    string outName;
    options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(), mesh->rank, cns->frame++);
    
    cnsPlotVTU(cns, fname);
  }

}
