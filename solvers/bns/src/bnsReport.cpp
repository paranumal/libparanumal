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

#include "bns.hpp"

void bns_t::Report(dfloat time, int tstep){

  static int frame=0;

  //compute vorticity
  vorticityKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_Dmatrices,
                  o_q, c, o_Vort);

  //compute q.M*q
  MassMatrixKernel(mesh.Nelements, mesh.o_ggeo, mesh.o_MM, o_q, o_Mq);

  dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
  dfloat norm2 = sqrt(linAlg.innerProd(Nentries, o_q, o_Mq, comm));

  if(mesh.rank==0)
    printf("%5.2f (%d), %5.2f (time, timestep, norm)\n", time, tstep, norm2);

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

    // copy data back to host
    o_q.copyTo(q);
    o_Vort.copyTo(Vort);

    // output field files
    string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);

    PlotFields(q, Vort, fname);
  }

  /*
  if(bns->dim==3){
    if(options.compareArgs("OUTPUT FILE FORMAT","ISO")){

      for (int gr=0; gr<bns->isoGNgroups; gr++){

        bns->isoNtris[0] = 0;
        bns->o_isoNtris.copyFrom(bns->isoNtris);
        if(mesh->nonPmlNelements){
        bns->isoSurfaceKernel(mesh->nonPmlNelements,    // Numner of elements
                              mesh->o_nonPmlElementIds,    // Element Ids
                              bns->isoField,               // which field to use for isosurfacing
                              bns->isoColorField,          // which field to use for isosurfacing
                              bns->isoGNlevels[gr],        // number of isosurface levels
                              bns->o_isoGLvalues[gr],      // array of isosurface levels
                              bns->isoMaxNtris,            // maximum number of generated triangles
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              bns->o_q,
                              bns->o_Vort,
                              bns->o_VortMag,
                              bns->o_plotInterp,
                              bns->o_plotEToV,
                              bns->o_isoNtris,             // output: number of generated triangles
                              bns->o_isoq);                // output: (p_dim+p_Nfields)*3*isoNtris[0] values (x,y,z,q0,q1..)

      }
        // find number of generated triangles
        bns->o_isoNtris.copyTo(bns->isoNtris);
        bns->isoNtris[0] = mymin(bns->isoNtris[0], bns->isoMaxNtris);

        //
        printf("Rank:%2d Group:%2d Triangles:%8d\n", mesh->rank, bns->isoNtris[0], gr);
        //
        int offset = 0;
        bns->o_isoq.copyTo(bns->isoq, bns->isoNtris[0]*(mesh->dim+bns->isoNfields)*3*sizeof(dfloat), offset);

        char fname[BUFSIZ];
        string outName;
        options.getArgs("OUTPUT FILE NAME", outName);


        if(options.compareArgs("OUTPUT FILE FORMAT", "WELD"))
        {
          int Ntris1 = bns->isoNtris[0];
          int Ntris2 = bnsWeldTriVerts(bns, Ntris1, bns->isoq);

          printf("Welding triangles:%8d to:%8d\n", Ntris1, Ntris2);
          sprintf(fname, "%s_%d_%d_ %04d_%04d.vtu",(char*)outName.c_str(), bns->isoField, gr, mesh->rank, bns->frame);
          bnsIsoWeldPlotVTU(bns,  fname);
        }
        else
        {
          sprintf(fname, "%s_%d_%d_ %04d_%04d.vtu",(char*)outName.c_str(), bns->isoField, gr, mesh->rank, bns->frame);
          bnsIsoPlotVTU(bns, bns->isoNtris[0], bns->isoq, fname);
        }
      }
      bns->frame++;
    }
  }
  */
}
