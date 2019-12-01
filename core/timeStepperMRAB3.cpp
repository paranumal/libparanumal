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

#include "core.hpp"
#include "timeStepper.hpp"

namespace TimeStepper {

mrab3::mrab3(dlong Nelements, dlong NhaloElements,
               int Np, int _Nfields, solver_t& _solver):
  timeStepper_t(Nelements, NhaloElements, Np, _Nfields, _solver),
  mesh(solver.mesh),
  Nlevels(mesh.mrNlevels),
  Nfields(_Nfields) {

  Nstages = 3;

  dfloat *rhsq0 = (dfloat*) calloc(N, sizeof(dfloat));
  o_rhsq0 = device.malloc(N*sizeof(dfloat), rhsq0);
  free(rhsq0);

  dfloat *rhsq = (dfloat*) calloc((Nstages-1)*N, sizeof(dfloat));
  o_rhsq = device.malloc((Nstages-1)*N*sizeof(dfloat), rhsq);
  free(rhsq);

  o_fQM = device.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Nfp
                          *mesh.Nfaces*Nfields*sizeof(dfloat));

  occa::properties kernelInfo = props; //copy base occa properties from solver

  kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;
  kernelInfo["defines/" "p_Nstages"] = Nstages;
  kernelInfo["defines/" "p_Np"] = mesh.Np;
  kernelInfo["defines/" "p_Nfp"] = mesh.Nfp;
  kernelInfo["defines/" "p_Nfaces"] = mesh.Nfaces;
  kernelInfo["defines/" "p_Nfields"] = Nfields;
  int maxNodes = mymax(mesh.Np, mesh.Nfp*mesh.Nfaces);
  kernelInfo["defines/" "p_maxNodes"] = maxNodes;

  updateKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                    "timeStepperMRAB.okl",
                                    "mrabUpdate",
                                    kernelInfo, comm);
  traceUpdateKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                    "timeStepperMRAB.okl",
                                    "mrabTraceUpdate",
                                    kernelInfo, comm);

  // initialize AB time stepping coefficients
  dfloat _ab_a[Nstages*Nstages] = {
                           1.0,      0.0,    0.0,
                         3./2.,   -1./2.,    0.0,
                       23./12., -16./12., 5./12.};
  dfloat _ab_b[Nstages*Nstages] = {
                         1./2.,      0.0,    0.0,
                         5./8.,   -1./8.,    0.0,
                       17./24.,  -7./24., 2./24.};

  ab_a = (dfloat*) calloc(Nstages*Nstages, sizeof(dfloat));
  ab_b = (dfloat*) calloc(Nstages*Nstages, sizeof(dfloat));
  memcpy(ab_a, _ab_a, Nstages*Nstages*sizeof(dfloat));
  memcpy(ab_b, _ab_b, Nstages*Nstages*sizeof(dfloat));

  shiftIndex = (int*) occaHostMallocPinned(device, Nlevels*sizeof(int),
                                       NULL, o_shiftIndex, h_shiftIndex);

  mrdt = (dfloat*) calloc(Nlevels, sizeof(dfloat));
  o_mrdt = device.malloc(Nlevels*sizeof(dfloat), mrdt);

  o_ab_a = device.malloc(Nstages*Nstages*sizeof(dfloat), ab_a);
  o_ab_b = device.malloc(Nstages*Nstages*sizeof(dfloat), ab_b);
}

void mrab3::Run(occa::memory &o_q, dfloat start, dfloat end) {

  dfloat time = start;

  //set timesteps and shifting index
  for (int lev=0;lev<Nlevels;lev++) {
    mrdt[lev] = dt*(1 << lev);
    shiftIndex[lev] = 0;
  }
  o_mrdt.copyFrom(mrdt);
  o_shiftIndex.copyFrom(shiftIndex);

  solver.Report(time,0);

  dfloat outputInterval;
  settings.getSetting("OUTPUT INTERVAL", outputInterval);

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
    Step(o_q, time, dt, order);
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

void mrab3::Step(occa::memory &o_q, dfloat time, dfloat _dt, int order) {

  occa::memory o_A = o_ab_a+order*Nstages*sizeof(dfloat);
  occa::memory o_B = o_ab_b+order*Nstages*sizeof(dfloat);

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
        shiftIndex[l] = (shiftIndex[l]+Nstages-2)%(Nstages-1);

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

    // o_shiftIndex.copyFrom(shiftIndex, "async: true");
    o_shiftIndex.copyFrom(shiftIndex); //Required to keep the update kernel overlapping the transfer,
                                       // but why does that happen?
  }
}

mrab3::~mrab3() {
  if (o_rhsq.size()) o_rhsq.free();
  if (o_fQM.size()) o_fQM.free();

  if (ab_a) free(ab_a);
  if (ab_b) free(ab_b);

  updateKernel.free();
  traceUpdateKernel.free();
}

/**************************************************/
/* PML version                                    */
/**************************************************/

mrab3_pml::mrab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
               int Np, int _Nfields, int _Npmlfields, solver_t& _solver):
  mrab3(Nelements, NhaloElements, Np, _Nfields, _solver),
  Npml(NpmlElements*Np*_Npmlfields),
  Npmlfields(_Npmlfields) {

  if (Npml) {
    dfloat *pmlq = (dfloat*) calloc(Npml, sizeof(dfloat));
    o_pmlq = device.malloc(Npml*sizeof(dfloat), pmlq);
    free(pmlq);

    dfloat *rhspmlq0 = (dfloat*) calloc(Npml, sizeof(dfloat));
    o_rhspmlq0 = device.malloc(Npml*sizeof(dfloat), rhspmlq0);
    free(rhspmlq0);

    dfloat *rhspmlq = (dfloat*) calloc((Nstages-1)*Npml, sizeof(dfloat));
    o_rhspmlq = device.malloc((Nstages-1)*Npml*sizeof(dfloat), rhspmlq);
    free(rhspmlq);

    occa::properties kernelInfo = props; //copy base occa properties from solver

    kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;
    kernelInfo["defines/" "p_Nstages"] = Nstages;
    kernelInfo["defines/" "p_Np"] = mesh.Np;
    kernelInfo["defines/" "p_Nfp"] = mesh.Nfp;
    kernelInfo["defines/" "p_Nfaces"] = mesh.Nfaces;
    kernelInfo["defines/" "p_Nfields"] = Nfields;
    int maxNodes = mymax(mesh.Np, mesh.Nfp*mesh.Nfaces);
    kernelInfo["defines/" "p_maxNodes"] = maxNodes;

    pmlUpdateKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                      "timeStepperMRAB.okl",
                                      "mrabPmlUpdate",
                                      kernelInfo, comm);
  }
}

void mrab3_pml::Step(occa::memory &o_q, dfloat time, dfloat _dt, int order) {

  occa::memory o_A = o_ab_a+order*Nstages*sizeof(dfloat);
  occa::memory o_B = o_ab_b+order*Nstages*sizeof(dfloat);

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
        shiftIndex[l] = (shiftIndex[l]+Nstages-2)%(Nstages-1);

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

    // o_shiftIndex.copyFrom(shiftIndex, "async: true");
    o_shiftIndex.copyFrom(shiftIndex); //Required to keep the update kernel overlapping the transfer,
                                       // but why does that happen?
  }
}

mrab3_pml::~mrab3_pml() {
  if (o_pmlq.size()) o_pmlq.free();
  if (o_rhspmlq0.size()) o_rhspmlq0.free();
  if (o_rhspmlq.size()) o_rhspmlq.free();

  pmlUpdateKernel.free();
}

} //namespace TimeStepper
