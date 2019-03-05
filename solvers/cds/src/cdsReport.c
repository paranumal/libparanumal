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

#include "cds.h"

void cdsReport(cds_t *cds, dfloat time, int tstep){

  mesh_t *mesh = cds->mesh;
  
  // copy data back to host
  cds->o_S.copyTo(cds->S);
  cds->o_NS.copyTo(cds->NS);
  // do error stuff on host
  cdsError(cds, time);

#ifdef RENDER
  if(cds->options.compareArgs("OUTPUT FILE FORMAT","PPM")){

    // copy data back to host
    cds->o_P.copyTo(cds->P);
    cds->o_Vort.copyTo(cds->Vort);
   
    //
    char fname[BUFSIZ];
    string outName;
    cds->options.getArgs("OUTPUT FILE NAME", outName);
    cdsRenderQuad3D(cds, (char*)outName.c_str(), cds->frame++);
  }
#endif
  
  if(cds->options.compareArgs("OUTPUT TYPE","VTU")){ 
    // output field files
    char fname[BUFSIZ];
    string outName;
    cds->options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(), mesh->rank, cds->frame++);
    cdsPlotVTU(cds, fname);
  }
}

