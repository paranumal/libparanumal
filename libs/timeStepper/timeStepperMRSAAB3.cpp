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
#include <complex>

namespace libp {

namespace TimeStepper {

using std::complex;

mrsaab3::mrsaab3(dlong _Nelements, dlong _NhaloElements,
             int _Np, int _Nfields,
             memory<dfloat> _lambda,
             platform_t& _platform, mesh_t& _mesh):
  timeStepperBase_t(_Nelements, _NhaloElements, _Np, _Nfields,
                    _platform, _mesh.comm),
  mesh(_mesh),
  Nlevels(mesh.mrNlevels),
  Nfields(_Nfields) {

  lambda.malloc(Nfields);
  lambda.copyFrom(_lambda);

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
                                    "timeStepperMRSAAB.okl",
                                    "mrsaabUpdate",
                                    kernelInfo);
  traceUpdateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperMRSAAB.okl",
                                    "mrsaabTraceUpdate",
                                    kernelInfo);

  saab_x.malloc(Nlevels*Nfields);
  saab_a.malloc(Nlevels*Nfields*Nstages*Nstages);
  saab_b.malloc(Nlevels*Nfields*Nstages*Nstages);

  h_shiftIndex = platform.hostMalloc<int>(Nlevels);
  o_shiftIndex = platform.malloc<int>(Nlevels);

  mrdt.malloc(Nlevels, 0.0);
  o_mrdt = platform.malloc<dfloat>(mrdt);

  o_saab_x = platform.malloc<dfloat>(Nlevels*Nfields);
  o_saab_a = platform.malloc<dfloat>(Nlevels*Nfields*Nstages*Nstages);
  o_saab_b = platform.malloc<dfloat>(Nlevels*Nfields*Nstages*Nstages);
}

void mrsaab3::Run(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat start, dfloat end) {

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

  //Compute coefficients
  UpdateCoefficients();

  // Populate Trace Buffer
  traceUpdateKernel(mesh.mrNelements[Nlevels-1],
                    mesh.o_mrElements[Nlevels-1],
                    mesh.o_mrLevel,
                    mesh.o_vmapM,
                    N,
                    o_shiftIndex,
                    o_mrdt,
                    o_saab_x,
                    o_saab_b,
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

void mrsaab3::Step(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat time, dfloat _dt, int order) {

  deviceMemory<dfloat> o_A = o_saab_a+order*Nstages;
  deviceMemory<dfloat> o_B = o_saab_b+order*Nstages;

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
                   o_saab_x,
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
                        o_saab_x,
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


void mrsaab3::UpdateCoefficients() {

  const int Nr = 32;
  complex<dfloat> R[Nr];

  //contour integral integration points
  for(int ind=1; ind <= Nr; ++ind){
    const dfloat theta = (dfloat) (ind - 0.5) / (dfloat) Nr;
    complex<dfloat> z(0., M_PI* theta);
    R[ind-1] = exp(z);
  }

  for(int lev=0; lev<Nlevels; ++lev){

    for (int n=0;n<Nfields;n++) {

      if (lambda[n]==0) { //Zero exponential term, usual AB coefficients

        dfloat _saab_X[1]  = { 1.0 };
        dfloat _saab_A[Nstages*Nstages]
                      = {    1.0,      0.0,    0.0,
                           3./2.,   -1./2.,    0.0,
                         23./12., -16./12., 5./12.};
        dfloat _saab_B[Nstages*Nstages] = {
                           1./2.,      0.0,    0.0,
                           5./8.,   -1./8.,    0.0,
                         17./24.,  -7./24., 2./24.};

        saab_x.copyFrom(_saab_X,               1, n                +lev*Nfields);
        saab_a.copyFrom(_saab_A, Nstages*Nstages, n*Nstages*Nstages+lev*Nfields*Nstages*Nstages);
        saab_b.copyFrom(_saab_B, Nstages*Nstages, n*Nstages*Nstages+lev*Nfields*Nstages*Nstages);

      } else {

        //Compute coefficients via contour integral
        dfloat alpha = lambda[n]*dt*(1<<lev);

        // Initialize complex variables for contour integral
        complex<double> a11(0.,0.);
        complex<double> a21(0.,0.);
        complex<double> a22(0.,0.);
        complex<double> a31(0.,0.);
        complex<double> a32(0.,0.);
        complex<double> a33(0.,0.);

        complex<double> b11(0.,0.);
        complex<double> b21(0.,0.);
        complex<double> b22(0.,0.);
        complex<double> b31(0.,0.);
        complex<double> b32(0.,0.);
        complex<double> b33(0.,0.);

        for(int i = 0; i<Nr; ++i ){
          complex<double> lr = alpha  + R[i];

          a11 +=  (exp(lr) - 1.)/lr;
          b11 +=  (exp(lr/2.) - 1.)/lr;

          a21 +=  (-2.*lr + (1.+lr)*exp(lr) - 1.)/pow(lr,2);
          a22 +=  (lr - exp(lr) + 1.)/pow(lr,2);
          b21 +=  (-1.5*lr + (1.+lr)*exp(lr/2.) - 1.)/pow(lr,2);
          b22 +=  (0.5*lr - exp(lr/2.) + 1.)/pow(lr,2);

          a31 += (-2.5*lr - 3.*pow(lr,2) + (1.+pow(lr,2)+1.5*lr)*exp(lr) - 1.)/pow(lr,3);
          a32 += (4.*lr + 3.*pow(lr,2)- (2.*lr + 2.0)*exp(lr) + 2.)/pow(lr,3);
          a33 +=-(1.5*lr + pow(lr,2)- (0.5*lr + 1.)*exp(lr) + 1.)/pow(lr,3);
          b31 += (exp(lr/2.)- 2.*lr - (15.*pow(lr,2))/8. + (pow(lr,2) + 1.5*lr)*exp(lr/2.) - 1.)/pow(lr,3);
          b32 += (3.*lr - 2.*exp(lr/2.0) + 1.25*pow(lr,2) - 2.*lr*exp(lr/2.) + 2.)/pow(lr,3);
          b33 +=-(lr - exp(lr/2.) + 0.375*pow(lr,2) - 0.5*lr*exp(lr/2.) + 1.)/pow(lr,3);
        }


        dfloat aa11=real(a11)/ (double) Nr;
        dfloat bb11=real(b11)/ (double) Nr;

        dfloat aa21=real(a21)/ (double) Nr;
        dfloat aa22=real(a22)/ (double) Nr;
        dfloat bb21=real(b21)/ (double) Nr;
        dfloat bb22=real(b22)/ (double) Nr;

        dfloat aa31=real(a31)/ (double) Nr;
        dfloat aa32=real(a32)/ (double) Nr;
        dfloat aa33=real(a33)/ (double) Nr;

        dfloat bb31=real(b31)/ (double) Nr;
        dfloat bb32=real(b32)/ (double) Nr;
        dfloat bb33=real(b33)/ (double) Nr;

        dfloat _saab_X[1]  = { std::exp(alpha) };
        dfloat _saab_A[Nstages*Nstages]
                        ={   aa11,   0.0,   0.0,
                             aa21,  aa22,   0.0,
                             aa31,  aa32,  aa33 };
        dfloat _saab_B[Nstages*Nstages]
                        ={   bb11,   0.0,   0.0,
                             bb21,  bb22,   0.0,
                             bb31,  bb32,  bb33 };

        saab_x.copyFrom(_saab_X,               1, n                +lev*Nfields);
        saab_a.copyFrom(_saab_A, Nstages*Nstages, n*Nstages*Nstages+lev*Nfields*Nstages*Nstages);
        saab_b.copyFrom(_saab_B, Nstages*Nstages, n*Nstages*Nstages+lev*Nfields*Nstages*Nstages);
      }
    }

    // move data to platform
    o_saab_x.copyFrom(saab_x);
    o_saab_a.copyFrom(saab_a);
    o_saab_b.copyFrom(saab_b);
  }
}

/**************************************************/
/* PML version                                    */
/**************************************************/

mrsaab3_pml::mrsaab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
                         int Np, int _Nfields, int _Npmlfields,
                         memory<dfloat> _lambda,
                         platform_t& _platform, mesh_t& _mesh):
  mrsaab3(Nelements, NhaloElements,
          Np, _Nfields, _lambda, _platform, _mesh),
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
                                      "timeStepperMRSAAB.okl",
                                      "mrsaabPmlUpdate",
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

    pmlsaab_a.malloc(Nstages*Nstages);
    pmlsaab_b.malloc(Nstages*Nstages);
    pmlsaab_a.copyFrom(_ab_a);
    pmlsaab_b.copyFrom(_ab_b);

    o_pmlsaab_a = platform.malloc<dfloat>(pmlsaab_a);
    o_pmlsaab_b = platform.malloc<dfloat>(pmlsaab_b);
  }
}

void mrsaab3_pml::Step(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat time, dfloat _dt, int order) {

  deviceMemory<dfloat> o_A = o_saab_a+order*Nstages;
  deviceMemory<dfloat> o_B = o_saab_b+order*Nstages;

  deviceMemory<dfloat> o_pmlA;
  if (Npml) o_pmlA = o_pmlsaab_a+order*Nstages;

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
                   o_saab_x,
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
                     o_pmlA,
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
                        o_saab_x,
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
