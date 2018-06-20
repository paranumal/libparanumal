#include "bns.h"

void bnsReport(bns_t *bns, dfloat time, setupAide &options){

mesh_t *mesh = bns->mesh; 

  bns->vorticityKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_Dmatrices,
                       bns->o_q,
                       bns->o_Vort);

  // report ramp function
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // copy data back to host
  bns->o_q.copyTo(bns->q);
  bns->o_Vort.copyTo(bns->Vort);

  if(rank==0){
    dfloat fx, fy, fz, intfx, intfy, intfz;
    bnsBodyForce(time, &fx, &fy, &fz, &intfx, &intfy, &intfz);
    printf("t: %g (fx,fy,fz) = (%g,%g,%g), int(fx,fy,fz) = (%g,%g,%g)\n",
     time, fx,fy,fz, intfx, intfy, intfz);
  }
  

  if(options.compareArgs("OUTPUT FILE FORMAT","VTU")){
   
    //
    char fname[BUFSIZ];
    string outName;
    options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(), rank, bns->frame++);
    bnsPlotVTU(bns, fname);
  }

  if(bns->dim==3){
    if(options.compareArgs("OUTPUT FILE FORMAT","ISO")){


      // for(int fld = 0; fld<bns->isoNfields; fld++){
      //   // Change Isocontour field
      //   bns->isoField = fld; 
         // Initialize
      bns->isoNtris[0] = 0; 
      bns->o_isoNtris.copyFrom(bns->isoNtris);


      bns->isoSurfaceKernel(mesh->nonPmlNelements,       // Numner of elements 
                            mesh->o_nonPmlElementIds,    // Element Ids
                            bns->isoField,               // which field to use for isosurfacing
                            bns->isoColorField,          // which field to use for isosurfacing
                            bns->isoNlevels,             // number of isosurface levels
                            bns->o_isoLevels,            // array of isosurface levels
                            bns->isoMaxNtris,            // maximum number of generated triangles
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            bns->o_q,
                            bns->o_Vort,
                            bns->o_plotInterp,
                            bns->o_plotEToV,
                            bns->o_isoNtris,  // output: number of generated triangles
                            bns->o_isoq       // output: (p_dim+p_Nfields)*3*isoNtris[0] values (x,y,z,q0,q1..)
                     );


      // find number of generated triangles
      bns->o_isoNtris.copyTo(bns->isoNtris);
      bns->isoNtris[0] = mymin(bns->isoNtris[0], bns->isoMaxNtris);

      printf("generated %d triangles\n", bns->isoNtris[0]);

      //
      int offset = 0;
      bns->o_isoq.copyTo(bns->isoq, bns->isoNtris[0]*(mesh->dim+bns->isoNfields)*3*sizeof(dfloat), offset);

      char fname[BUFSIZ];
      string outName;
      options.getArgs("OUTPUT FILE NAME", outName);

      if(options.compareArgs("OUTPUT FILE FORMAT", "WELD")){

        int Ntris1 = bns->isoNtris[0];
        int Ntris2 = bnsWeldTriVerts(bns, Ntris1, bns->isoq);

        printf("Ntri1: %d and Ntris2:%d \n", Ntris1, Ntris2);

        // // int procid   = gethostid(); // Processor id for gmsh file. 
        // int plotnum  = bns->frame;
        // int N_offset = 0;          // Gmsh mpi node offset
        // int E_offset = 0;          // Gmsh mpi element offset
        // double plottime = 0.0;     // record time
        // bool bBinary = false;      // toggle binary/ascii
        // // bool bBinary = true;      // toggle binary/ascii
        // int tstep = plotnum;       // dummy time-step

        // sprintf(fname, "%s_%04d_%04d.msh",(char*)outName.c_str(), rank, bns->frame++);
        // bnsIsoPlotGmsh(bns, Ntris2, fname, bns->tstep, N_offset, E_offset, plotnum, plottime, bBinary);
     
        sprintf(fname, "%s_%d_%04d_%04d.vtu",(char*)outName.c_str(), bns->isoField, rank, bns->frame++);
        bnsIsoWeldPlotVTU(bns,  fname);
      }else{
        sprintf(fname, "%s_%d_%04d_%04d.vtu",(char*)outName.c_str(), bns->isoField, rank, bns->frame++);
        bnsIsoPlotVTU(bns, bns->isoNtris[0], bns->isoq, fname);
      }



    // }
   


    }
  }




  if(options.compareArgs("OUTPUT FILE FORMAT","TEC")){ 
    // //boltzmannComputeVorticity2D(mesh, mesh->q,5, mesh->Nfields);
    // char fname[BUFSIZ];
    // sprintf(fname, "foo_v2_%04d.dat",rank);
    // bnsPlotTEC(bns, fname, t);
  }
  
  
}
