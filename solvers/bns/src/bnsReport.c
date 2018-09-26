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

#include "bns.h"

void bnsReport(bns_t *bns, dfloat time, setupAide &options){

mesh_t *mesh = bns->mesh; 

  bns->vorticityKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_Dmatrices,
                       bns->o_q,
                       bns->o_Vort,
                       bns->o_VortMag);

#if 0
  if(bns->dim==3){
    ogsGatherScatter(bns->o_VortMag, ogsDfloat, ogsAdd, mesh->ogs);  
    int Ntotal = mesh->Np*mesh->Nelements;
    bns->dotMultiplyKernel(Ntotal, bns->o_VortMag, mesh->ogs->o_invDegree); 
  }
#endif
  
  // report ramp function
  if(mesh->rank==0){
    dfloat fx, fy, fz, intfx, intfy, intfz;
    bnsBodyForce(time, &fx, &fy, &fz, &intfx, &intfy, &intfz);
    printf("t: %g (fx,fy,fz) = (%g,%g,%g), int(fx,fy,fz) = (%g,%g,%g)\n",
     time, fx,fy,fz, intfx, intfy, intfz);
  }

#ifdef RENDER
  if(options.compareArgs("OUTPUT FILE FORMAT","PPM")){

    // copy data back to host
    bns->o_q.copyTo(bns->q);
    bns->o_Vort.copyTo(bns->Vort);
    bns->o_VortMag.copyTo(bns->VortMag);
   
    //
    char fname[BUFSIZ];
    string outName;
    options.getArgs("OUTPUT FILE NAME", outName);
    bnsRenderQuad3D(bns, (char*)outName.c_str(), bns->frame++);
  }
#endif
  
  if(options.compareArgs("OUTPUT FILE FORMAT","VTU")){

    // copy data back to host
    bns->o_q.copyTo(bns->q);
    bns->o_Vort.copyTo(bns->Vort);
    bns->o_VortMag.copyTo(bns->VortMag);
   
    //
    char fname[BUFSIZ];
    string outName;
    options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(), mesh->rank, bns->frame++);
    bnsPlotVTU(bns, fname);
  }

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

          #if 0 // Currently we dont use gmsh postprocessing

          // int procid   = gethostid(); // Processor id for gmsh file. 
          int plotnum  = bns->frame;
          int N_offset = 0;          // Gmsh mpi node offset
          int E_offset = 0;          // Gmsh mpi element offset
          double plottime = 0.0;     // record time
          bool bBinary = false;      // toggle binary/ascii
          // bool bBinary = true;      // toggle binary/ascii
          int tstep = plotnum;       // dummy time-step
          sprintf(fname, "%s_%04d_%04d.msh",(char*)outName.c_str(), mesh->rank, bns->frame++);
          bnsIsoPlotGmsh(bns, Ntris2, fname, bns->tstep, N_offset, E_offset, plotnum, plottime, bBinary);
          #endif
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


  // if(options.compareArgs("OUTPUT FILE FORMAT","TEC")){ 
  //   // //boltzmannComputeVorticity2D(mesh, mesh->q,5, mesh->Nfields);
  //   // char fname[BUFSIZ];
  //   // sprintf(fname, "foo_v2_%04d.dat",rank);
  //   // bnsPlotTEC(bns, fname, t);
  // }
  
  
}
