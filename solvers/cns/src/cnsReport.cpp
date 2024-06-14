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

  static int outFrame=0;
  static int forceFrame=0;
  if(settings.compareSetting("REPORT FORCES","TRUE")){

    if(forceFrame==0){
      momentCenter.malloc(mesh.dim, 0.0); 
      std::string stateStr;
      settings.getSetting("MOMENT CENTER", stateStr);
      tokenizer(mesh.dim, stateStr, momentCenter,  ',');
      o_momentCenter = platform.malloc<dfloat>(momentCenter); 
    }
    reportForces(time, tstep); 
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

    // reportSmoothFields(time, tstep); 

   if(stabType!=Stab::NOSTAB){
    o_qdetect.copyTo(qdetect); 
    o_viscosity.copyTo(viscosity); 
      // stab.Apply(o_q, o_gradq, 0.0);
      // stab.Report(time, tstep);
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
    // {
    // // sprintf(fname, ".dat");
    // FILE *fp;
    // fp = fopen("maxDensity.dat", "a");
    // dfloat maxr = 1000.0; 
    // for(int e=0; e<mesh.Nelements; e++){
    //   for(int n=0; n<mesh.Np; n++){
    //    maxr = maxr<q[e*Nfields*mesh.Np + n+ 0*mesh.Np] ? maxr: q[e*Nfields*mesh.Np + n+ 0*mesh.Np];
    //   }
    // }
    // fprintf(fp, " %g %g\n", time, maxr);
    // fclose(fp); 
    // }

    sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, outFrame++);
    PlotFields(q, Vort, std::string(fname));
  }
}

void cns_t::reportSmoothFields(dfloat time, int tstep){

  deviceMemory<dfloat> o_qv = platform.reserve<dfloat>(mesh.Nelements*mesh.Np*Nfields);
  deviceMemory<dfloat> o_qf = platform.reserve<dfloat>(mesh.Nelements*mesh.Np*Nfields);

  // arrange memory hit for vector-gather-scatter
  reportArrangeLayoutKernel(mesh.Nelements, 0, o_q, o_qv); 

  ogs.GatherScatter(o_qv, Nfields, ogs::Add, ogs::Sym);

  reportAverageKernel(mesh.Nelements, o_weight, o_qv); 

  reportArrangeLayoutKernel(mesh.Nelements, 1, o_qv, o_qf);

  o_qf.copyTo(q);  

}


void cns_t::reportForces(dfloat time, int tstep){

  // Write out the integrated pressure and viscous forces
  // Note that this is not the most efficient way but we can use the kernels
  // Viscous and pressure forces are seperated
  dlong Nentries = 0; 
  // Forces and moments
  Nentries = mesh.dim==2 ? mesh.Nelements*mesh.Np*(mesh.dim*mesh.dim+1):mesh.Nelements*mesh.Np*(mesh.dim*mesh.dim+mesh.dim); 
  deviceMemory<dfloat> o_F     = platform.reserve<dfloat>(Nentries);
  deviceMemory<dfloat> o_gradq = platform.reserve<dfloat>(mesh.Nelements*mesh.Np*mesh.dim*mesh.dim);

  // compute volume contributions to gradients
  forcesVolumeKernel(mesh.Nelements,
                     mesh.o_vgeo,
                     mesh.o_D,
                     o_q,
                     o_gradq);


  // compute volume contributions to gradients
  forcesSurfaceKernel(mesh.Nelements,
                   mesh.o_sgeo,
                   mesh.o_sM,
                   mesh.o_vmapM,
                   mesh.o_EToB,
                   mesh.o_x,
                   mesh.o_y,
                   mesh.o_z,
                   o_momentCenter, 
                   o_pCoeff,
                   o_q,
                   o_gradq,
                   o_F);
  // output field files
  std::string name;
  settings.getSetting("OUTPUT FILE NAME", name);
  name = name + "_forces.dat";
  // Open file
  FILE *fp; 
  fp = fopen(name.c_str(), "a");

  const dlong shift = mesh.Nelements*mesh.Np; 
  if(mesh.dim==2){
    const dfloat vFx = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+0*shift , mesh.comm); 
    const dfloat vFy = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+1*shift , mesh.comm); 
    const dfloat pFx = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+2*shift , mesh.comm); 
    const dfloat pFy = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+3*shift , mesh.comm);
    const dfloat Mz  = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+4*shift , mesh.comm);
    if(mesh.rank==0){
      printf("forces: %.4e %.4e %.4e %.4e %.4e\n", vFx, vFy, pFx, pFy, Mz);
      fprintf(fp,"%.6e %.6e %.6e %.6e %.6e %.6e\n", time, vFx, vFy, pFx, pFy, Mz);
    } 
  }else{
    const dfloat vFx = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+0*shift , mesh.comm); 
    const dfloat vFy = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+1*shift , mesh.comm); 
    const dfloat vFz = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+2*shift , mesh.comm); 
    
    const dfloat pFx = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+3*shift , mesh.comm); 
    const dfloat pFy = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+4*shift , mesh.comm); 
    const dfloat pFz = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+5*shift , mesh.comm); 

    const dfloat Mx = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+6*shift , mesh.comm); 
    const dfloat My = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+7*shift , mesh.comm); 
    const dfloat Mz = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+8*shift , mesh.comm); 
    if(mesh.rank==0){
      printf("forces: %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n", time, vFx, vFy, vFz, pFx, pFy, pFz, Mx, My, Mz);
      fprintf(fp,"%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n", time, vFx, vFy, vFz, pFx, pFy, pFz, Mx, My, Mz);

    } 

  }
  
  fclose(fp); 
}