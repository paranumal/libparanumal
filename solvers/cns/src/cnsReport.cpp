/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "cns.hpp"

void cns_t::Report(dfloat time, int tstep){

  static int outFrame = 0;
  static int forceFrame = 0;     
  if(settings.compareSetting("REPORT FORCES","TRUE")){
    writeForces(time, tstep, forceFrame); 
    forceFrame++; 
  }

  //compute q.M*q
  dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
  deviceMemory<dfloat> o_Mq = platform.reserve<dfloat>(Nentries);
  mesh.MassMatrixApply(o_q, o_Mq);

  dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm));
  o_Mq.free();

  if(mesh.rank==0)
    printf("%5.2f (%d), %5.2f (time, timestep, norm)\n", time, tstep, norm2);

  
  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

   if(stabType!=Stab::NOSTAB){
      o_qdetect.copyTo(qdetect); 
      o_viscosity.copyTo(viscosity); 
    }
    
    // copy data back to host
    o_q.copyTo(q);
    
    //compute vorticity
    deviceMemory<dfloat> o_Vort = platform.reserve<dfloat>(mesh.dim*mesh.Nelements*mesh.Np);
    vorticityKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D, o_q, o_Vort);

    memory<dfloat> Vort(mesh.dim*mesh.Nelements*mesh.Np);

  
    o_Vort.copyTo(Vort);

    // output field files
    std::string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, outFrame++);
    PlotFields(q, Vort, std::string(fname));
  }
}


void cns_t::writeForces(dfloat time, int tstep, int frame){

  if(mesh.dim==2 && mesh.rank==0){
    printf("----------------------------------------------------------------------\n"); 
    if(reportComponent){
       printf("viscous_x - viscous_y - pressure_x - pressure_y - moment_z \n"); 
    }else{
       printf("force_x - force_y - moment_z \n"); 
    }
  }else if(mesh.dim==3 && mesh.rank==0){
    printf("----------------------------------------------------------------------\n"); 
    if(reportComponent){
       printf("viscous_x-viscous_y-viscous_z-pressure_x-pressure_y-pressure_z-moment_x-moment_y-moment_z\n"); 
    }else{
       printf("force_x - force_y - force_z - moment_x - moment_y - moment_z\n"); 
    }
  }

  // output field files
  std::string name;
  settings.getSetting("OUTPUT FILE NAME", name);
  name = name + "analysis.dat";
  // Open file
  FILE *fp; 
  fp = fopen(name.c_str(), "a");
  if(frame==0 && mesh.rank==0){
    if(reportComponent){
      if(mesh.dim==2){
       fprintf(fp, "time\tgroupID\tviscous_x\tviscous_y\tpressure_x\tpressure_y\tmoment_z\n"); 
      }else{
       fprintf(fp, "time\tgroupID\tviscous_x\tviscous_y\tviscous_z\tpressure_x\tpressure_y\tpressure_z\tmoment_x\tmoment_y\tmoment_z\n"); 
      }
    }else{
      if(mesh.dim==2){
       fprintf(fp, "time\tgroupID\tforce_x\tforce_y\tmoment_z\n"); 
      }else{
       fprintf(fp, "time\tgroupID\tforce_x\tforce_y\tforce_z\tmoment_x\tmoment_y\tmoment_z\n"); 
      }
    }
  }


   // Write out the integrated pressure and viscous forces
  dlong Nentries = mesh.dim==2 ? mesh.Nelements*mesh.Np*(mesh.dim+mesh.dim+1):
                                 mesh.Nelements*mesh.Np*(mesh.dim+mesh.dim+mesh.dim); 
  
  // Compute all forces on all boundaries
  deviceMemory<dfloat> o_F     = platform.reserve<dfloat>(Nentries);
#if 1
  deviceMemory<dfloat> o_gradq = platform.reserve<dfloat>(mesh.Nelements*mesh.Np*mesh.dim*mesh.dim);

  // compute volume contributions to gradients
  forcesVolumeKernel(mesh.Nelements,
                     mesh.o_vgeo,
                     mesh.o_D,
                     o_q,
                     o_gradq);



#else
dlong NlocalGrads = mesh.Nelements*mesh.Np*Ngrads;
dlong NhaloGrads  = mesh.totalHaloPairs*mesh.Np*Ngrads;
deviceMemory<dfloat> o_gradq = platform.reserve<dfloat>(NlocalGrads+NhaloGrads);  

// extract q trace halo and start exchange
  fieldTraceHalo.ExchangeStart(o_q, 1);
 gradVolumeKernel(mesh.Nelements,
                   mesh.o_vgeo,
                   mesh.o_D,
                   o_q,
                   o_gradq);    

fieldTraceHalo.ExchangeFinish(o_q, 1);
gradSurfaceKernel(mesh.Nelements,
                    BCStateID, 
                    mesh.o_sgeo,
                    mesh.o_LIFT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    o_pCoeff, 
                    o_flowStates, 
                    time,
                    o_q,
                    o_gradq);

#endif





  // Write out every force/moment components for every report group 
  for(int grp =0; grp<NreportGroups; grp++){
    // compute volume contributions to gradients
    forcesSurfaceKernel(mesh.Nelements,
                       grp, 
                       mesh.o_sgeo,
                       mesh.o_sM,
                       mesh.o_vmapM,
                       o_EToG,
                       o_reportGroups,
                       mesh.o_x,
                       mesh.o_y,
                       mesh.o_z,
                       o_momentCenter, 
                       o_pCoeff,
                       o_q,
                       o_gradq,
                       o_F);


     const dlong shift = mesh.Nelements*mesh.Np; 
    if(mesh.dim==2){
      const dfloat vFx = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+0*shift , mesh.comm); 
      const dfloat vFy = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+1*shift , mesh.comm); 
      const dfloat pFx = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+2*shift , mesh.comm); 
      const dfloat pFy = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+3*shift , mesh.comm);
      const dfloat Mz  = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+4*shift , mesh.comm);
      if(mesh.rank==0){
        if(reportComponent){
          printf("on report group %d : %.2e %.2e %.2e %.2e %.2e\n", grp, vFx, vFy, pFx, pFy, Mz);
          fprintf(fp,"%.6e %d %.6e %.6e %.6e %.6e %.6e\n", time, grp, vFx, vFy, pFx, pFy, Mz);
        }else{
          printf("on report group %d : %.2e %.2e %.2e\n", grp, vFx+pFx, vFy + pFy, Mz);
          fprintf(fp,"%.6e %d %.6e %.6e %.6e\n", time, grp, vFx+pFx, vFy+pFy, Mz);
        }
      } 
  }else{
    const dfloat vFx = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F + 0*shift, mesh.comm); 
    const dfloat vFy = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F + 1*shift, mesh.comm); 
    const dfloat vFz = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F + 2*shift, mesh.comm); 
    
    const dfloat pFx = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F + 3*shift, mesh.comm); 
    const dfloat pFy = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F + 4*shift, mesh.comm); 
    const dfloat pFz = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F + 5*shift, mesh.comm); 

    const dfloat Mx  = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F + 6*shift, mesh.comm); 
    const dfloat My  = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F + 7*shift, mesh.comm); 
    const dfloat Mz  = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F + 8*shift, mesh.comm); 
    if(mesh.rank==0){
      if(reportComponent){
          printf("on report group %d : %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e\n", grp,  vFx, vFy, vFz, pFx, pFy, pFz, Mx, My, Mz);
          fprintf(fp, "%.6e %d %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e\n", time, grp,  vFx, vFy, vFz, pFx, pFy, pFz, Mx, My, Mz);
        }else{
          printf("on report group %d : %.2e %.2e %.2e %.2e %.2e %.2e\n", grp, vFx+pFx, vFy+pFy, vFz+pFz, Mx, My, Mz);
          fprintf(fp, "%.6e %d %.2e %.2e %.2e %.2e %.2e %.2e\n", time, grp, vFx+pFx, vFy+pFy, vFz+pFz, Mx, My, Mz);
      }
    } 
  }

  }

  if(mesh.rank==0){
    printf("----------------------------------------------------------------------\n");     
  }
  
  fclose(fp); 
}