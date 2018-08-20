#include "mppf.h"

void mppfReport(mppf_t *mppf, dfloat time, int tstep){

  mesh_t *mesh = mppf->mesh;

  mppf->vorticityKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_Dmatrices,
                       mppf->fieldOffset,
                       mppf->o_U,
                       mppf->o_Vort);

  mppf->divergenceVolumeKernel(mesh->Nelements,
                             mesh->o_vgeo,
                             mesh->o_Dmatrices,
                             mppf->fieldOffset,
                             mppf->o_U,
                             mppf->o_Div);

  // // gather-scatter the vorticity
  // ellipticParallelGatherScatter(mesh, mesh->ogs, mppf->o_Vort, dfloatString, "add");  
  // mppf->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, mppf->o_Vort, mppf->o_Vort);
  
  // copy data back to host
  mppf->o_Phi.copyTo(mppf->Phi);
  mppf->o_U.copyTo(mppf->U); 
  mppf->o_GU.copyTo(mppf->GU); 
  mppf->o_P.copyTo(mppf->P);
  mppf->o_Rho.copyTo(mppf->Rho);

  mppf->o_Mu.copyTo(mppf->Mu);
  // mppf->o_Psi.copyTo(mppf->Mu);

  mppf->o_NPhi.copyTo(mppf->NPhi);

  mppf->o_Vort.copyTo(mppf->Vort);
  mppf->o_Div.copyTo(mppf->Div);

  // do error stuff on host
  mppfError(mppf, time);

  if(mppf->options.compareArgs("OUTPUT TYPE","VTU")){ 
    // output field files
    char fname[BUFSIZ];
    string outName;
    mppf->options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(), mesh->rank, mppf->frame++);

    mppfPlotVTU(mppf, fname);
  }

  // if(mppf->options.compareArgs("OUTPUT TYPE","ISO") && (mppf->dim==3)){ 

  //    for (int gr=0; gr<mppf->isoGNgroups; gr++){

  //       mppf->isoNtris[0] = 0; 
  //       mppf->o_isoNtris.copyFrom(mppf->isoNtris);

  //       mppf->isoSurfaceKernel(mesh->Nelements,
  //                             mppf->fieldOffset,    
  //                             mppf->isoField,               
  //                             mppf->isoColorField,          
  //                             mppf->isoGNlevels[gr],        
  //                             mppf->o_isoGLvalues[gr],      
  //                             mppf->isoMaxNtris,            
  //                             mesh->o_x,
  //                             mesh->o_y,
  //                             mesh->o_z,
  //                             mppf->o_P, 
  //                             mppf->o_U,
  //                             mppf->o_Vort,
  //                             mppf->o_plotInterp,
  //                             mppf->o_plotEToV,
  //                             mppf->o_isoNtris,             // output: number of generated triangles
  //                             mppf->o_isoq);                // output: (p_dim+p_Nfields)*3*isoNtris[0] values (x,y,z,q0,q1..)
                  

  //       // find number of generated triangles
  //       mppf->o_isoNtris.copyTo(mppf->isoNtris);
  //       mppf->isoNtris[0] = mymin(mppf->isoNtris[0], mppf->isoMaxNtris);
  //       // 
  //       printf("Rank:%2d Group:%2d Triangles:%8d\n", mesh->rank, gr, mppf->isoNtris[0]);
  //       //
  //       int offset = 0;
  //       mppf->o_isoq.copyTo(mppf->isoq, mppf->isoNtris[0]*(mesh->dim+mppf->isoNfields)*3*sizeof(dfloat), offset);

  //       int Ntris1 = mppf->isoNtris[0];
  //       int Ntris2 = mppfWeldTriVerts(mppf, Ntris1, mppf->isoq);

  //       printf("Welding triangles:%8d to:%8d\n", Ntris1, Ntris2);

  //       char fname[BUFSIZ];
  //       string outName;
  //       mppf->options.getArgs("OUTPUT FILE NAME", outName);
  //       sprintf(fname, "%s_%d_%d_ %04d_%04d.vtu",(char*)outName.c_str(), mppf->isoField, gr, mesh->rank, mppf->frame++);
  //       mppfIsoPlotVTU(mppf,  fname);
  //     }
  // }

}

