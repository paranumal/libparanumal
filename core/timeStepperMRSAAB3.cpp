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
#include <complex>

namespace TimeStepper {

using std::complex;

mrsaab3::mrsaab3(dlong _Nelements, dlong _NhaloElements,
             int _Np, int _Nfields,
             dfloat *_lambda, solver_t& _solver):
  timeStepper_t(_Nelements, _NhaloElements, _Np, _Nfields, _solver),
  mesh(solver.mesh),
  Nlevels(mesh.mrNlevels),
  Nfields(_Nfields) {

  lambda = (dfloat *) malloc(Nfields*sizeof(dfloat));
  memcpy(lambda, _lambda, Nfields*sizeof(dfloat));

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
                                    "timeStepperMRSAAB.okl",
                                    "mrsaabUpdate",
                                    kernelInfo, comm);
  traceUpdateKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                    "timeStepperMRSAAB.okl",
                                    "mrsaabTraceUpdate",
                                    kernelInfo, comm);

  saab_x = (dfloat*) calloc(Nlevels*Nfields, sizeof(dfloat));
  saab_a = (dfloat*) calloc(Nlevels*Nfields*Nstages*Nstages, sizeof(dfloat));
  saab_b = (dfloat*) calloc(Nlevels*Nfields*Nstages*Nstages, sizeof(dfloat));

  shiftIndex = (int*) occaHostMallocPinned(device, Nlevels*sizeof(int),
                                       NULL, o_shiftIndex, h_shiftIndex);

  mrdt = (dfloat*) calloc(Nlevels, sizeof(dfloat));
  o_mrdt = device.malloc(Nlevels*sizeof(dfloat), mrdt);

  o_saab_x = device.malloc(Nlevels*Nfields*sizeof(dfloat));
  o_saab_a = device.malloc(Nlevels*Nfields*Nstages*Nstages*sizeof(dfloat));
  o_saab_b = device.malloc(Nlevels*Nfields*Nstages*Nstages*sizeof(dfloat));
}

void mrsaab3::Run(occa::memory &o_q, dfloat start, dfloat end) {

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

void mrsaab3::Step(occa::memory &o_q, dfloat time, dfloat _dt, int order) {

  occa::memory o_A = o_saab_a+order*Nstages*sizeof(dfloat);
  occa::memory o_B = o_saab_b+order*Nstages*sizeof(dfloat);

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
                        o_saab_x,
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

        memcpy(saab_x+n                +lev*Nfields,
              _saab_X, 1*sizeof(dfloat));
        memcpy(saab_a+n*Nstages*Nstages+lev*Nfields*Nstages*Nstages,
              _saab_A,Nstages*Nstages*sizeof(dfloat));
        memcpy(saab_b+n*Nstages*Nstages+lev*Nfields*Nstages*Nstages,
              _saab_B,Nstages*Nstages*sizeof(dfloat));

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

        dfloat _saab_X[1]  = { exp(alpha) };
        dfloat _saab_A[Nstages*Nstages]
                        ={   aa11,   0.0,   0.0,
                             aa21,  aa22,   0.0,
                             aa31,  aa32,  aa33 };
        dfloat _saab_B[Nstages*Nstages]
                        ={   bb11,   0.0,   0.0,
                             bb21,  bb22,   0.0,
                             bb31,  bb32,  bb33 };

        memcpy(saab_x+n                +lev*Nfields,
              _saab_X, 1*sizeof(dfloat));
        memcpy(saab_a+n*Nstages*Nstages+lev*Nfields*Nstages*Nstages,
              _saab_A,Nstages*Nstages*sizeof(dfloat));
        memcpy(saab_b+n*Nstages*Nstages+lev*Nfields*Nstages*Nstages,
              _saab_B,Nstages*Nstages*sizeof(dfloat));
      }
    }

    // move data to device
    o_saab_x.copyFrom(saab_x);
    o_saab_a.copyFrom(saab_a);
    o_saab_b.copyFrom(saab_b);
  }
}

mrsaab3::~mrsaab3() {
  if (o_rhsq.size()) o_rhsq.free();
  if (o_fQM.size()) o_fQM.free();

  if (saab_x) free(saab_x);
  if (saab_a) free(saab_a);
  if (saab_b) free(saab_b);

  if (o_saab_x.size()) o_saab_x.free();
  if (o_saab_a.size()) o_saab_a.free();
  if (o_saab_b.size()) o_saab_b.free();

  updateKernel.free();
  traceUpdateKernel.free();
}

/**************************************************/
/* PML version                                    */
/**************************************************/

mrsaab3_pml::mrsaab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
                         int Np, int _Nfields, int _Npmlfields,
                         dfloat *_lambda, solver_t& _solver):
  mrsaab3(Nelements, NhaloElements, Np, _Nfields, _lambda, _solver),
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
                                      "timeStepperMRSAAB.okl",
                                      "mrsaabPmlUpdate",
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

    pmlsaab_a = (dfloat*) calloc(Nstages*Nstages, sizeof(dfloat));
    pmlsaab_b = (dfloat*) calloc(Nstages*Nstages, sizeof(dfloat));
    memcpy(pmlsaab_a, _ab_a, Nstages*Nstages*sizeof(dfloat));
    memcpy(pmlsaab_b, _ab_b, Nstages*Nstages*sizeof(dfloat));

    o_pmlsaab_a = device.malloc(Nstages*Nstages*sizeof(dfloat), pmlsaab_a);
    o_pmlsaab_b = device.malloc(Nstages*Nstages*sizeof(dfloat), pmlsaab_b);
  }
}

void mrsaab3_pml::Step(occa::memory &o_q, dfloat time, dfloat _dt, int order) {

  occa::memory o_A = o_saab_a+order*Nstages*sizeof(dfloat);
  occa::memory o_B = o_saab_b+order*Nstages*sizeof(dfloat);

  occa::memory o_pmlA;
  if (Npml) o_pmlA = o_pmlsaab_a+order*Nstages*sizeof(dfloat);

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
                        o_saab_x,
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

mrsaab3_pml::~mrsaab3_pml() {
  if (o_pmlq.size()) o_pmlq.free();
  if (o_rhspmlq0.size()) o_rhspmlq0.free();
  if (o_rhspmlq.size()) o_rhspmlq.free();

  pmlUpdateKernel.free();
}

} //namespace TimeStepper