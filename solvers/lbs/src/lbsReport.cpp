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
  #if 1

  static int frame=0;
  // Compute velocity and density
  momentsKernel(mesh.Nelements, o_LBM, o_q, o_U); 

  // //compute vorticity
  // vorticityKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D, o_U, o_Vort);

  //compute q.M*q
  mesh.MassMatrixApply(o_U, o_Mq);

  dlong Nentries = mesh.Nelements*mesh.Np*Nmacro;
  dfloat norm2 = sqrt(platform.linAlg.innerProd(Nentries, o_q, o_Mq, mesh.comm));

  if(mesh.rank==0)
    printf("%5.2f (%d), %5.2f (time, timestep, norm)\n", time, tstep, norm2);

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
  #endif

  /*
  if(lbs->dim==3){
    if(options.compareArgs("OUTPUT FILE FORMAT","ISO")){

      for (int gr=0; gr<lbs->isoGNgroups; gr++){

        lbs->isoNtris[0] = 0;
        lbs->o_isoNtris.copyFrom(lbs->isoNtris);
        if(mesh->nonPmlNelements){
        lbs->isoSurfaceKernel(mesh->nonPmlNelements,    // Numner of elements
                              mesh->o_nonPmlElementIds,    // Element Ids
                              lbs->isoField,               // which field to use for isosurfacing
                              lbs->isoColorField,          // which field to use for isosurfacing
                              lbs->isoGNlevels[gr],        // number of isosurface levels
                              lbs->o_isoGLvalues[gr],      // array of isosurface levels
                              lbs->isoMaxNtris,            // maximum number of generated triangles
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              lbs->o_q,
                              lbs->o_Vort,
                              lbs->o_VortMag,
                              lbs->o_plotInterp,
                              lbs->o_plotEToV,
                              lbs->o_isoNtris,             // output: number of generated triangles
                              lbs->o_isoq);                // output: (p_dim+p_Nfields)*3*isoNtris[0] values (x,y,z,q0,q1..)

      }
        // find number of generated triangles
        lbs->o_isoNtris.copyTo(lbs->isoNtris);
        lbs->isoNtris[0] = mymin(lbs->isoNtris[0], lbs->isoMaxNtris);

        //
        printf("Rank:%2d Group:%2d Triangles:%8d\n", mesh->rank, lbs->isoNtris[0], gr);
        //
        int offset = 0;
        lbs->o_isoq.copyTo(lbs->isoq, lbs->isoNtris[0]*(mesh->dim+lbs->isoNfields)*3*sizeof(dfloat), offset);

        char fname[BUFSIZ];
        string outName;
        options.getArgs("OUTPUT FILE NAME", outName);


        if(options.compareArgs("OUTPUT FILE FORMAT", "WELD"))
        {
          int Ntris1 = lbs->isoNtris[0];
          int Ntris2 = lbsWeldTriVerts(lbs, Ntris1, lbs->isoq);

          printf("Welding triangles:%8d to:%8d\n", Ntris1, Ntris2);
          sprintf(fname, "%s_%d_%d_ %04d_%04d.vtu",(char*)outName.c_str(), lbs->isoField, gr, mesh->rank, lbs->frame);
          lbsIsoWeldPlotVTU(lbs,  fname);
        }
        else
        {
          sprintf(fname, "%s_%d_%d_ %04d_%04d.vtu",(char*)outName.c_str(), lbs->isoField, gr, mesh->rank, lbs->frame);
          lbsIsoPlotVTU(lbs, lbs->isoNtris[0], lbs->isoq, fname);
        }
      }
      lbs->frame++;
    }
  }
  */
}
