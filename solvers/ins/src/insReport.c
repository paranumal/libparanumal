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

  // gather-scatter the vorticity
  ellipticParallelGatherScatter(mesh, mesh->ogs, ins->o_Vort, dfloatString, "add");  
  ins->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, ins->o_Vort, ins->o_Vort);
  
  // copy data back to host
  ins->o_U.copyTo(ins->U);
  ins->o_P.copyTo(ins->P);

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

