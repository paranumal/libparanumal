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
  static int momentFrame=0;

  if(settings.compareSetting("REPORT FORCES","TRUE")){
    reportForces(time, tstep); forceFrame++; 
  }

  if(settings.compareSetting("REPORT MOMENTS","TRUE")){
    if(momentFrame==0){
      momentCenter.malloc(mesh.dim, 0.0); 
      std::string stateStr;
      settings.getSetting("MOMENT CENTER", stateStr);
      tokenizer(mesh.dim, stateStr, momentCenter,  ',');
      o_momentCenter = platform.malloc<dfloat>(momentCenter); 
    }
    reportMoments(time, tstep); momentFrame++; 
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
   // if(stab.type!=Stab::NOSTAB){
   //    stab.Apply(o_q, o_q, 0.0);
   //    stab.Report(0.0, tstep);}

    //compute vorticity
    deviceMemory<dfloat> o_Vort = platform.reserve<dfloat>(mesh.dim*mesh.Nelements*mesh.Np);
    vorticityKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D, o_q, o_Vort);

    memory<dfloat> Vort(mesh.dim*mesh.Nelements*mesh.Np);

    // copy data back to host
    o_q.copyTo(q);
    o_Vort.copyTo(Vort);

    // output field files
    std::string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, outFrame++);

    PlotFields(q, Vort, std::string(fname));
  }
}




void cns_t::reportForces(dfloat time, int tstep){

  // Write out the integrated pressure and viscous forces
  // Note that this is not the most efficient way but we can use the kernels
  // Viscous and pressure forces are seperated
  dlong Nentries = mesh.Nelements*mesh.Np*(2*mesh.dim);
  deviceMemory<dfloat> o_F = platform.reserve<dfloat>(Nentries);
  // compute volume contributions to gradients
  computeForcesKernel(mesh.Nelements,
                   mesh.o_vgeo,
                   mesh.o_sgeo,
                   mesh.o_D,
                   mesh.o_sM,
                   mesh.o_vmapM,
                   mesh.o_EToB,
                   mesh.o_x,
                   mesh.o_y,
                   mesh.o_z,
                   o_pCoeff,
                   o_q,
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
    if(mesh.rank==0) 
      printf("forces: %.4e %.4e %.4e %.4e \n", vFx, vFy, pFx, pFy);

    fprintf(fp,"%.6e %.6e %.6e %.6e %.6e\n", time, vFx, vFy, pFx, pFy);
  }else{
    const dfloat vFx = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+0*shift , mesh.comm); 
    const dfloat vFy = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+1*shift , mesh.comm); 
    const dfloat vFz = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+2*shift , mesh.comm); 
    
    const dfloat pFx = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+3*shift , mesh.comm); 
    const dfloat pFy = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+4*shift , mesh.comm); 
    const dfloat pFz = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_F+5*shift , mesh.comm); 
    if(mesh.rank==0) 
      printf("forces: %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n", time, vFx, vFy, vFz, pFx, pFy, pFz);

    fprintf(fp,"%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n", time, vFx, vFy, vFz, pFx, pFy, pFz);
  }
  
  fclose(fp); 
}


void cns_t::reportMoments(dfloat time, int tstep){

  // Write out the integrated pressure and viscous forces
  // Note that this is not the most efficient way but we can use the kernels
  dlong Nentries = mesh.dim==2 ? mesh.Nelements*mesh.Np:mesh.Nelements*mesh.Np*mesh.dim;
  deviceMemory<dfloat> o_M = platform.reserve<dfloat>(Nentries);
  // compute volume contributions to gradients
  computeMomentsKernel(mesh.Nelements,
                       mesh.o_vgeo,
                       mesh.o_sgeo,
                       mesh.o_D,
                       mesh.o_sM,
                       mesh.o_vmapM,
                       mesh.o_EToB,
                       mesh.o_x,
                       mesh.o_y,
                       mesh.o_z,
                       o_momentCenter, 
                       o_pCoeff,
                       o_q,
                       o_M);

  // output field files
  std::string name;
  settings.getSetting("OUTPUT FILE NAME", name);
  name = name + "_moments.dat";
  // Open file
  FILE *fp; fp = fopen(name.c_str(), "a");

  const dlong shift = mesh.Nelements*mesh.Np; 
  
  if(mesh.dim==2){
    const dfloat Mz = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_M+0*shift , mesh.comm); 

    if(mesh.rank==0){
      printf("moments:%.4e \n", Mz);
    }

    fprintf(fp,"%.6e %.6e\n", time, Mz);
    fclose(fp); 
  }else{
    const dfloat  Mx = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_M+0*shift , mesh.comm); 
    const dfloat  My = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_M+1*shift , mesh.comm); 
    const dfloat  Mz = platform.linAlg().sum(mesh.Nelements*mesh.Np, o_M+2*shift , mesh.comm); 
    if(mesh.rank==0)
      printf("moments: %.4e %.4e %.4e\n", Mx, My, Mz);

    // write to file
    fprintf(fp, "%.6e %.6e %.6e %.6e\n", time, Mx, My, Mz);
  }
  
  fclose(fp); 
}

