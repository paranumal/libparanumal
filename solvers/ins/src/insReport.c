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

#include "ins.h"

void insReport(ins_t *ins, dfloat time, int tstep){

  mesh_t *mesh = ins->mesh;

  ins->vorticityKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_Dmatrices,
                       ins->fieldOffset,
                       ins->o_U,
                       ins->o_Vort);

  ins->divergenceVolumeKernel(mesh->Nelements,
                             mesh->o_vgeo,
                             mesh->o_Dmatrices,
                             ins->fieldOffset,
                             ins->o_U,
                             ins->o_Div);

  // gatherscatter vorticity field
  dlong Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
  for(int s=0; s<ins->dim; s++){
    ins->o_UH.copyFrom(ins->o_Vort,Ntotal*sizeof(dfloat),0,s*ins->fieldOffset*sizeof(dfloat));
  
    ogsGatherScatter(ins->o_UH, ogsDfloat, ogsAdd, mesh->ogs);  
    ins->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, ins->o_UH, ins->o_UH);
    
    ins->o_UH.copyTo(ins->o_Vort,Ntotal*sizeof(dfloat),s*ins->fieldOffset*sizeof(dfloat),0);
  }

  // gather-scatter divergence 
  ogsGatherScatter(ins->o_Div, ogsDfloat, ogsAdd, mesh->ogs);  
  ins->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, ins->o_Div, ins->o_Div);

  // copy data back to host
  ins->o_U.copyTo(ins->U);
  ins->o_P.copyTo(ins->P);
  if(ins->solveHeat)
    ins->o_T.copyTo(ins->T);

  ins->o_Vort.copyTo(ins->Vort);
  ins->o_Div.copyTo(ins->Div);

  // do error stuff on host
  insError(ins, time);

  if(ins->options.compareArgs("OUTPUT TYPE","VTU")){ 
    // output field files
    char fname[BUFSIZ];
    string outName;
    ins->options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(), mesh->rank, ins->frame++);

    insPlotVTU(ins, fname);
  }

  if(ins->options.compareArgs("OUTPUT TYPE","ISO") && (ins->dim==3)){ 

     for (int gr=0; gr<ins->isoGNgroups; gr++){

        ins->isoNtris[0] = 0; 
        ins->o_isoNtris.copyFrom(ins->isoNtris);

        ins->isoSurfaceKernel(mesh->Nelements,
                              ins->fieldOffset,    
                              ins->isoField,               
                              ins->isoColorField,          
                              ins->isoGNlevels[gr],        
                              ins->o_isoGLvalues[gr],      
                              ins->isoMaxNtris,            
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              ins->o_P, 
                              ins->o_U,
                              ins->o_Vort,
                              ins->o_plotInterp,
                              ins->o_plotEToV,
                              ins->o_isoNtris,             // output: number of generated triangles
                              ins->o_isoq);                // output: (p_dim+p_Nfields)*3*isoNtris[0] values (x,y,z,q0,q1..)
                  

        // find number of generated triangles
        ins->o_isoNtris.copyTo(ins->isoNtris);
        ins->isoNtris[0] = mymin(ins->isoNtris[0], ins->isoMaxNtris);
        // 
        printf("Rank:%2d Group:%2d Triangles:%8d\n", mesh->rank, gr, ins->isoNtris[0]);
        //
        int offset = 0;
        ins->o_isoq.copyTo(ins->isoq, ins->isoNtris[0]*(mesh->dim+ins->isoNfields)*3*sizeof(dfloat), offset);

        int Ntris1 = ins->isoNtris[0];
        int Ntris2 = insWeldTriVerts(ins, Ntris1, ins->isoq);

        printf("Welding triangles:%8d to:%8d\n", Ntris1, Ntris2);

        char fname[BUFSIZ];
        string outName;
        ins->options.getArgs("OUTPUT FILE NAME", outName);
        sprintf(fname, "%s_%d_%d_ %04d_%04d.vtu",(char*)outName.c_str(), ins->isoField, gr, mesh->rank, ins->frame++);
        insIsoPlotVTU(ins,  fname);
      }
  }

}

