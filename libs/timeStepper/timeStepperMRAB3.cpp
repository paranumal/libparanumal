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

#include "core.hpp"
#include "timeStepper.hpp"

namespace libp {

namespace TimeStepper {

mrab3::mrab3(dlong Nelements, dlong NhaloElements,
               int Np, int _Nfields,
               platform_t& _platform, mesh_t& _mesh):
  timeStepperBase_t(Nelements, NhaloElements,
                    Np, _Nfields, _platform, _mesh.comm),
  mesh(_mesh),
  Nlevels(mesh.mrNlevels),
  Nfields(_Nfields) {

  Nstages = 3;

  memory<dfloat> rhsq0(N, 0.0);
  o_rhsq0 = platform.malloc<dfloat>(rhsq0);

  memory<dfloat> rhsq((Nstages-1)*N, 0.0);
  o_rhsq = platform.malloc<dfloat>(rhsq);

  o_fQM = platform.malloc<dfloat>((mesh.Nelements+mesh.totalHaloPairs)*mesh.Nfp
                                  *mesh.Nfaces*Nfields);

  properties_t kernelInfo = platform.props(); //copy base occa properties from solver

  const int blocksize=256;

  kernelInfo["defines/" "p_blockSize"] = blocksize;
  kernelInfo["defines/" "p_Nstages"] = Nstages;
  kernelInfo["defines/" "p_Np"] = mesh.Np;
  kernelInfo["defines/" "p_Nfp"] = mesh.Nfp;
  kernelInfo["defines/" "p_Nfaces"] = mesh.Nfaces;
  kernelInfo["defines/" "p_Nfields"] = Nfields;
  int maxNodes = std::max(mesh.Np, mesh.Nfp*mesh.Nfaces);
  kernelInfo["defines/" "p_maxNodes"] = maxNodes;

  updateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperMRAB.okl",
                                    "mrabUpdate",
                                    kernelInfo);
  traceUpdateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperMRAB.okl",
                                    "mrabTraceUpdate",
                                    kernelInfo);

  // initialize AB time stepping coefficients
  dfloat _ab_a[Nstages*Nstages] = {
                           1.0,      0.0,    0.0,
                         3./2.,   -1./2.,    0.0,
                       23./12., -16./12., 5./12.};
  dfloat _ab_b[Nstages*Nstages] = {
                         1./2.,      0.0,    0.0,
                         5./8.,   -1./8.,    0.0,
                       17./24.,  -7./24., 2./24.};

  ab_a.malloc(Nstages*Nstages);
  ab_b.malloc(Nstages*Nstages);
  ab_a.copyFrom(_ab_a);
  ab_b.copyFrom(_ab_b);

  h_shiftIndex = platform.hostMalloc<int>(Nlevels);
  o_shiftIndex = platform.malloc<int>(Nlevels);

  mrdt.malloc(Nlevels, 0.0);
  o_mrdt = platform.malloc<dfloat>(mrdt);

  o_ab_a = platform.malloc<dfloat>(ab_a);
  o_ab_b = platform.malloc<dfloat>(ab_b);
}

void mrab3::Run(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat start, dfloat end) {

  dfloat time = start;

  //set timesteps and shifting index
  for (int lev=0;lev<Nlevels;lev++) {
    mrdt[lev] = dt*(1 << lev);
    h_shiftIndex[lev] = 0;
  }
  o_mrdt.copyFrom(mrdt);
  h_shiftIndex.copyTo(o_shiftIndex);

  solver.Report(time,0);

  dfloat outputInterval=0.0;
  solver.settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat outputTime = time + outputInterval;

  // Populate Trace Buffer
  traceUpdateKernel(mesh.mrNelements[Nlevels-1],
                    mesh.o_mrElements[Nlevels-1],
                    mesh.o_mrLevel,
                    mesh.o_vmapM,
                    N,
                    o_shiftIndex,
                    o_mrdt,
                    o_ab_b,
                    o_rhsq0,
                    o_rhsq,
                    o_q,
                    o_fQM);

  dfloat DT = dt*(1 << (Nlevels-1));

  int tstep=0;
  int order=0;
  while (time < end) {
    Step(solver, o_q, time, dt, order);
    time += DT;
    tstep++;
    if (order<Nstages-1) order++;

    if (time>outputTime) {
      //report state
      solver.Report(outputTime,tstep);
      outputTime += outputInterval;
    }
  }
}

void mrab3::Step(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat time, dfloat _dt, int order) {

  deviceMemory<dfloat> o_A = o_ab_a+order*Nstages;
  deviceMemory<dfloat> o_B = o_ab_b+order*Nstages;

  for (int Ntick=0; Ntick < (1 << (Nlevels-1));Ntick++) {

    // intermediate stage time
    dfloat currentTime = time + dt*Ntick;

    int lev=0;
    for (;lev<Nlevels-1;lev++)
      if (Ntick % (1<<(lev+1)) != 0) break; //find the max lev to compute rhs

    //evaluate ODE rhs = f(q,t)
    solver.rhsf_MR(o_q, o_rhsq0, o_fQM, currentTime, lev);

    for (lev=0;lev<Nlevels-1;lev++)
      if ((Ntick+1) % (1<<(lev+1)) !=0) break; //find the max lev to update

    // update all elements of level <= lev
    if (mesh.mrNelements[lev])
      updateKernel(mesh.mrNelements[lev],
                   mesh.o_mrElements[lev],
                   mesh.o_mrLevel,
                   mesh.o_vmapM,
                   N,
                   o_shiftIndex,
                   o_mrdt,
                   o_A,
                   o_rhsq0,
                   o_rhsq,
                   o_fQM,
                   o_q);

    //rotate index
    if (Nstages>2)
      for (int l=0; l<=lev; l++)
        h_shiftIndex[l] = (h_shiftIndex[l]+Nstages-2)%(Nstages-1);

    //compute intermediate trace values on lev+1 / lev interface
    if (lev+1<Nlevels && mesh.mrInterfaceNelements[lev+1])
      traceUpdateKernel(mesh.mrInterfaceNelements[lev+1],
                        mesh.o_mrInterfaceElements[lev+1],
                        mesh.o_mrLevel,
                        mesh.o_vmapM,
                        N,
                        o_shiftIndex,
                        o_mrdt,
                        o_B,
                        o_rhsq0,
                        o_rhsq,
                        o_q,
                        o_fQM);

    // o_shiftIndex.copyFrom(h_shiftIndex, properties_t("async", true));
    h_shiftIndex.copyTo(o_shiftIndex); //Required to keep the update kernel overlapping the transfer,
                                       // but why does that happen?
  }
}

/**************************************************/
/* PML version                                    */
/**************************************************/

mrab3_pml::mrab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
                     int Np, int _Nfields, int _Npmlfields,
                     platform_t& _platform, mesh_t& _mesh):
  mrab3(Nelements, NhaloElements,
        Np, _Nfields, _platform, _mesh),
  Npml(NpmlElements*Np*_Npmlfields),
  Npmlfields(_Npmlfields) {

  if (Npml) {
    memory<dfloat> pmlq(Npml, 0.0);
    o_pmlq = platform.malloc<dfloat>(pmlq);

    memory<dfloat> rhspmlq0(Npml, 0.0);
    o_rhspmlq0 = platform.malloc<dfloat>(rhspmlq0);

    memory<dfloat> rhspmlq((Nstages-1)*Npml, 0.0);
    o_rhspmlq = platform.malloc<dfloat>(rhspmlq);

    properties_t kernelInfo = platform.props(); //copy base occa properties from solver

    const int blocksize=256;

    kernelInfo["defines/" "p_blockSize"] = blocksize;
    kernelInfo["defines/" "p_Nstages"] = Nstages;
    kernelInfo["defines/" "p_Np"] = mesh.Np;
    kernelInfo["defines/" "p_Nfp"] = mesh.Nfp;
    kernelInfo["defines/" "p_Nfaces"] = mesh.Nfaces;
    kernelInfo["defines/" "p_Nfields"] = Nfields;
    int maxNodes = std::max(mesh.Np, mesh.Nfp*mesh.Nfaces);
    kernelInfo["defines/" "p_maxNodes"] = maxNodes;

    pmlUpdateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                      "timeStepperMRAB.okl",
                                      "mrabPmlUpdate",
                                      kernelInfo);
  }
}

void mrab3_pml::Step(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat time, dfloat _dt, int order) {

  deviceMemory<dfloat> o_A = o_ab_a+order*Nstages;
  deviceMemory<dfloat> o_B = o_ab_b+order*Nstages;

  for (int Ntick=0; Ntick < (1 << (Nlevels-1));Ntick++) {

    // intermediate stage time
    dfloat currentTime = time + dt*Ntick;

    int lev=0;
    for (;lev<Nlevels-1;lev++)
      if (Ntick % (1<<(lev+1)) != 0) break; //find the max lev to compute rhs

    //evaluate ODE rhs = f(q,t)
    solver.rhsf_MR_pml(o_q, o_pmlq,
                       o_rhsq0, o_rhspmlq0,
                       o_fQM, currentTime, lev);

    for (lev=0;lev<Nlevels-1;lev++)
      if ((Ntick+1) % (1<<(lev+1)) !=0) break; //find the max lev to update

    // update all elements of level <= lev
    if (mesh.mrNelements[lev])
      updateKernel(mesh.mrNelements[lev],
                   mesh.o_mrElements[lev],
                   mesh.o_mrLevel,
                   mesh.o_vmapM,
                   N,
                   o_shiftIndex,
                   o_mrdt,
                   o_A,
                   o_rhsq0,
                   o_rhsq,
                   o_fQM,
                   o_q);

    if (mesh.mrNpmlElements[lev])
      pmlUpdateKernel(mesh.mrNpmlElements[lev],
                     mesh.o_mrPmlElements[lev],
                     mesh.o_mrPmlIds[lev],
                     mesh.o_mrLevel,
                     Npml,
                     Npmlfields,
                     o_shiftIndex,
                     o_mrdt,
                     o_A,
                     o_rhspmlq0,
                     o_rhspmlq,
                     o_pmlq);

    //rotate index
    if (Nstages>2)
      for (int l=0; l<=lev; l++)
        h_shiftIndex[l] = (h_shiftIndex[l]+Nstages-2)%(Nstages-1);

    //compute intermediate trace values on lev+1 / lev interface
    if (lev+1<Nlevels && mesh.mrInterfaceNelements[lev+1])
      traceUpdateKernel(mesh.mrInterfaceNelements[lev+1],
                        mesh.o_mrInterfaceElements[lev+1],
                        mesh.o_mrLevel,
                        mesh.o_vmapM,
                        N,
                        o_shiftIndex,
                        o_mrdt,
                        o_B,
                        o_rhsq0,
                        o_rhsq,
                        o_q,
                        o_fQM);

    // o_shiftIndex.copyFrom(h_shiftIndex, properties_t("async", true));
    h_shiftIndex.copyTo(o_shiftIndex); //Required to keep the update kernel overlapping the transfer,
                                       // but why does that happen?
  }
}

} //namespace TimeStepper

} //namespace libp
