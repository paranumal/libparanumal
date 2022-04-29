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

#ifndef TIMESTEPPER_HPP
#define TIMESTEPPER_HPP

#include "core.hpp"
#include "settings.hpp"
#include "mesh.hpp"
#include "solver.hpp"

namespace libp {

//forward declare
namespace TimeStepper { class timeStepperBase_t; }

/* General TimeStepper object*/
class timeStepper_t {
 public:
  timeStepper_t() = default;

  /*Generic setup. Create a Stepper object and wrap it in a shared_ptr*/
  template<class Stepper, class... Args>
  void Setup(Args&& ... args) {
    ts = std::make_shared<Stepper>(args...);
  }

  void Run(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat start, dfloat end);

  void SetTimeStep(dfloat dt_);

  dfloat GetTimeStep();

  dfloat GetGamma();

 private:
  std::shared_ptr<TimeStepper::timeStepperBase_t> ts=nullptr;

  void assertInitialized();
};

namespace TimeStepper {

//base time stepper class
class timeStepperBase_t {
public:
  platform_t platform;
  comm_t comm;

  dlong N;
  dlong Nhalo;

  dfloat dt;

  timeStepperBase_t(dlong Nelements, dlong NhaloElements,
                    int Np, int Nfields,
                    platform_t& _platform, comm_t _comm):
    platform(_platform),
    comm(_comm),
    N(Nelements*Np*Nfields),
    Nhalo(NhaloElements*Np*Nfields) {}

  virtual void Run(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat start, dfloat end)=0;

  void SetTimeStep(dfloat dt_) {dt = dt_;};

  dfloat GetTimeStep() {return dt;};

  virtual dfloat GetGamma() {
    LIBP_FORCE_ABORT("GetGamma() not available in this Timestepper");
    return 0.0;
  }
};

/* Adams Bashforth, order 3 */
class ab3: public timeStepperBase_t {
protected:
  int Nstages;
  int shiftIndex;

  memory<dfloat> ab_a;
  deviceMemory<dfloat> o_ab_a;

  deviceMemory<dfloat> o_rhsq;

  kernel_t updateKernel;

  virtual void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

public:
  ab3(dlong Nelements, dlong NhaloElements,
      int Np, int Nfields,
      platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat start, dfloat end);
};

/* Low-Storage Explicit Runge-Kutta, order 4 */
class lserk4: public timeStepperBase_t {
protected:
  int Nrk;
  memory<dfloat> rka, rkb, rkc;

  deviceMemory<dfloat> o_rhsq;
  deviceMemory<dfloat> o_resq;

  deviceMemory<dfloat> o_saveq;

  kernel_t updateKernel;

  virtual void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt);

public:
  lserk4(dlong Nelements, dlong NhaloElements,
         int Np, int Nfields,
         platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat start, dfloat end);
};

/* Dormand-Prince method */
/* Explict Runge-Kutta, order 5 with embedded order 4 and adaptive time-stepping */
class dopri5: public timeStepperBase_t {
protected:
  int Nrk;

  dlong Nblock;

  memory<dfloat> rkC, rkA, rkE;
  deviceMemory<dfloat> o_rkA, o_rkE;

  deviceMemory<dfloat> o_errtmp;
  pinnedMemory<dfloat> h_errtmp;

  deviceMemory<dfloat> o_rhsq;
  deviceMemory<dfloat> o_rkq;
  deviceMemory<dfloat> o_rkrhsq;
  deviceMemory<dfloat> o_rkerr;

  deviceMemory<dfloat> o_saveq;


  kernel_t rkUpdateKernel;
  kernel_t rkStageKernel;
  kernel_t rkErrorEstimateKernel;

  dfloat dtMIN; //minumum allowed timestep
  dfloat ATOL;  //absolute error tolerance
  dfloat RTOL;  //relative error tolerance
  dfloat safe;   //safety factor

  //error control parameters
  dfloat beta;
  dfloat factor1;
  dfloat factor2;

  dfloat exp1;
  dfloat invfactor1;
  dfloat invfactor2;
  dfloat facold;
  dfloat sqrtinvNtotal;

  virtual void Backup(deviceMemory<dfloat> &o_Q);
  virtual void Restore(deviceMemory<dfloat> &o_Q);
  virtual void AcceptStep(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rq);

  virtual void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt);

  virtual dfloat Estimater(deviceMemory<dfloat>& o_q);

public:
  dopri5(dlong Nelements, dlong NhaloElements,
         int Np, int Nfields,
         platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat start, dfloat end);
};

/* Semi-Analytic Adams-Bashforth, order 3 */
class saab3: public timeStepperBase_t {
protected:
  int Nstages;
  int shiftIndex;

  int Np, Nfields;
  dlong Nblock, Nelements, NhaloElements;

  memory<dfloat> lambda;

  pinnedMemory<dfloat> h_saab_x, h_saab_a;
  deviceMemory<dfloat> o_saab_x, o_saab_a;

  deviceMemory<dfloat> o_rhsq;

  kernel_t updateKernel;

  virtual void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

  virtual void UpdateCoefficients();

public:
  saab3(dlong _Nelements, dlong _NhaloElements,
        int _Np, int _Nfields,
        memory<dfloat> _lambda,
        platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat start, dfloat end);
};

/* Semi-Analytic Explict Runge-Kutta, order 4 with embedded order 3 and adaptive time-stepping */
class sark4: public timeStepperBase_t {
protected:
  int Nrk;
  int order, embeddedOrder;

  int Np, Nfields;
  dlong Nblock, Nelements, NhaloElements;

  memory<dfloat> lambda;

  memory<dfloat> rkC;
  deviceMemory<dfloat> o_rkX, o_rkA, o_rkE;
  pinnedMemory<dfloat> h_rkX, h_rkA, h_rkE;

  deviceMemory<dfloat> o_rhsq;
  deviceMemory<dfloat> o_rkq;
  deviceMemory<dfloat> o_rkrhsq;
  deviceMemory<dfloat> o_rkerr;

  deviceMemory<dfloat> o_saveq;

  deviceMemory<dfloat> o_errtmp;
  pinnedMemory<dfloat> h_errtmp;

  kernel_t rkUpdateKernel;
  kernel_t rkStageKernel;
  kernel_t rkErrorEstimateKernel;

  dfloat dtMIN; //minumum allowed timestep
  dfloat ATOL;  //absolute error tolerance
  dfloat RTOL;  //relative error tolerance
  dfloat safe;   //safety factor

  //error control parameters
  dfloat beta;
  dfloat factor1;
  dfloat factor2;

  dfloat exp1;
  dfloat invfactor1;
  dfloat invfactor2;
  dfloat facold;
  dfloat sqrtinvNtotal;

  virtual void Backup(deviceMemory<dfloat> &o_Q);
  virtual void Restore(deviceMemory<dfloat> &o_Q);
  virtual void AcceptStep(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rq);

  virtual void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt);

  dfloat Estimater(deviceMemory<dfloat>& o_q);

  void UpdateCoefficients();

public:
  sark4(dlong _Nelements, dlong _NhaloElements,
        int _Np, int _Nfields,
        memory<dfloat> _lambda,
        platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat start, dfloat end);
};

/* Semi-Analytic Explict Runge-Kutta, order 5 with embedded order 4 and adaptive time-stepping */
class sark5: public timeStepperBase_t {
protected:
  int Nrk;
  int order, embeddedOrder;

  int Np, Nfields;
  dlong Nblock, Nelements, NhaloElements;

  memory<dfloat> lambda;

  memory<dfloat> rkC;
  deviceMemory<dfloat> o_rkX, o_rkA, o_rkE;
  pinnedMemory<dfloat> h_rkX, h_rkA, h_rkE;


  deviceMemory<dfloat> o_rhsq;
  deviceMemory<dfloat> o_rkq;
  deviceMemory<dfloat> o_rkrhsq;
  deviceMemory<dfloat> o_rkerr;

  deviceMemory<dfloat> o_saveq;

  deviceMemory<dfloat> o_errtmp;
  pinnedMemory<dfloat> h_errtmp;

  kernel_t rkUpdateKernel;
  kernel_t rkStageKernel;
  kernel_t rkErrorEstimateKernel;

  dfloat dtMIN; //minumum allowed timestep
  dfloat ATOL;  //absolute error tolerance
  dfloat RTOL;  //relative error tolerance
  dfloat safe;   //safety factor

  //error control parameters
  dfloat beta;
  dfloat factor1;
  dfloat factor2;

  dfloat exp1;
  dfloat invfactor1;
  dfloat invfactor2;
  dfloat facold;
  dfloat sqrtinvNtotal;

  virtual void Backup(deviceMemory<dfloat> &o_Q);
  virtual void Restore(deviceMemory<dfloat> &o_Q);
  virtual void AcceptStep(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rq);

  virtual void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt);

  dfloat Estimater(deviceMemory<dfloat>& o_q);

  void UpdateCoefficients();

public:
  sark5(dlong _Nelements, dlong _NhaloElements,
        int _Np, int _Nfields,
        memory<dfloat> _lambda,
        platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat start, dfloat end);
};

/* Backward Difference Formula, order 3, with extrapolation */
class extbdf3: public timeStepperBase_t {
protected:
  int Nstages;
  int shiftIndex;

  memory<dfloat> extbdf_a;
  memory<dfloat> extbdf_b;
  deviceMemory<dfloat> o_extbdf_a;
  deviceMemory<dfloat> o_extbdf_b;

  deviceMemory<dfloat> o_rhs;
  deviceMemory<dfloat> o_qn;
  deviceMemory<dfloat> o_F;

  kernel_t rhsKernel;

  virtual void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

public:
  extbdf3(dlong Nelements, dlong NhaloElements,
      int Np, int Nfields,
      platform_t& _platform, comm_t _comm);

  dfloat GetGamma();

  void Run(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat start, dfloat end);
};

/* Backward Difference Formula, order 3, with subcycling */
class ssbdf3: public timeStepperBase_t {
protected:
  int Nstages;
  int shiftIndex;

  memory<dfloat> ssbdf_b;
  deviceMemory<dfloat> o_ssbdf_b;

  deviceMemory<dfloat> o_rhs;
  deviceMemory<dfloat> o_qn;
  deviceMemory<dfloat> o_qhat;

  kernel_t rhsKernel;

  virtual void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

public:
  ssbdf3(dlong Nelements, dlong NhaloElements,
      int Np, int Nfields,
      platform_t& _platform, comm_t _comm);

  dfloat GetGamma();

  void Run(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat start, dfloat end);
};

/* Multi-rate Adams-Bashforth, order 3 */
class mrab3: public timeStepperBase_t {
protected:
  mesh_t mesh;

  int Nstages;
  int Nlevels;
  int Nfields;

  deviceMemory<int> o_shiftIndex;
  deviceMemory<int> h_shiftIndex;

  memory<dfloat> mrdt;
  deviceMemory<dfloat> o_mrdt;

  memory<dfloat> ab_a, ab_b;
  deviceMemory<dfloat> o_ab_a, o_ab_b;

  deviceMemory<dfloat> o_rhsq0, o_rhsq, o_fQM;

  kernel_t updateKernel;
  kernel_t traceUpdateKernel;

  virtual void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

public:
  mrab3(dlong _Nelements, dlong _NhaloElements,
         int _Np, int _Nfields,
         platform_t& _platform, mesh_t& _mesh);

  void Run(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat start, dfloat end);
};

/* Multi-rate Semi-Analytic Adams-Bashforth, order 3 */
class mrsaab3: public timeStepperBase_t {
protected:
  mesh_t mesh;

  int Nstages;
  int Nlevels;
  int Nfields;

  memory<dfloat> lambda;

  deviceMemory<int> o_shiftIndex;
  pinnedMemory<int> h_shiftIndex;

  memory<dfloat> mrdt;
  deviceMemory<dfloat> o_mrdt;

  memory<dfloat> saab_x, saab_a, saab_b;
  deviceMemory<dfloat> o_saab_x, o_saab_a, o_saab_b;

  deviceMemory<dfloat> o_rhsq0, o_rhsq, o_fQM;

  kernel_t updateKernel;
  kernel_t traceUpdateKernel;

  virtual void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

  void UpdateCoefficients();

public:
  mrsaab3(dlong _Nelements, dlong _NhaloElements,
         int _Np, int _Nfields,
         memory<dfloat> _lambda,
         platform_t& _platform, mesh_t& _mesh);

  void Init();
  void Run(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat start, dfloat end);
};


/**************************************************/
/* Derived Time Integrators which step PML fields */
/**************************************************/

/* Adams Bashforth, order 3 */
class ab3_pml: public ab3 {
private:
  dlong Npml;

  deviceMemory<dfloat> o_pmlq;
  deviceMemory<dfloat> o_rhspmlq;

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

public:
  ab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
          int Np, int Nfields, int Npmlfields,
          platform_t& _platform, comm_t _comm);
};

/* Low-Storage Explicit Runge-Kutta, order 4 */
class lserk4_pml: public lserk4 {
private:
  dlong Npml;

  deviceMemory<dfloat> o_pmlq;
  deviceMemory<dfloat> o_rhspmlq;
  deviceMemory<dfloat> o_respmlq;

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt);

public:
  lserk4_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int Npmlfields,
            platform_t& _platform, comm_t _comm);
};

/* Dormand-Prince method */
/* Explict Runge-Kutta, order 5 with embedded order 4 and adaptive time-stepping */
class dopri5_pml: public dopri5 {
private:
  dlong Npml;

  deviceMemory<dfloat> o_pmlq;
  deviceMemory<dfloat> o_rhspmlq;
  deviceMemory<dfloat> o_rkpmlq;
  deviceMemory<dfloat> o_rkrhspmlq;

  deviceMemory<dfloat> o_savepmlq;

  kernel_t rkPmlUpdateKernel;

  void Backup(deviceMemory<dfloat> &o_Q);
  void Restore(deviceMemory<dfloat> &o_Q);
  void AcceptStep(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rq);

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt);

public:
  dopri5_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int Npmlfields,
            platform_t& _platform, comm_t _comm);
};

/* Semi-Analytic Adams-Bashforth, order 3 */
// Note: PML fields are not stepped Semi-Analytically
class saab3_pml: public saab3 {
private:
  dlong Npml;

  memory<dfloat> pmlsaab_x, pmlsaab_a;
  deviceMemory<dfloat> o_pmlsaab_x, o_pmlsaab_a;

  deviceMemory<dfloat> o_pmlq;
  deviceMemory<dfloat> o_rhspmlq;

  kernel_t pmlUpdateKernel;

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

public:
  saab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields,
            memory<dfloat> _lambda,
            platform_t& _platform, comm_t _comm);
};

/* Semi-Analytic Explict Runge-Kutta, order 4 with embedded order 3 and adaptive time-stepping */
// Note: PML fields are not stepped Semi-Analytically
class sark4_pml: public sark4 {
private:
  dlong Npml;

  memory<dfloat> pmlrkA;
  deviceMemory<dfloat> o_pmlrkA;

  deviceMemory<dfloat> o_pmlq;
  deviceMemory<dfloat> o_rhspmlq;
  deviceMemory<dfloat> o_rkpmlq;
  deviceMemory<dfloat> o_rkrhspmlq;

  deviceMemory<dfloat> o_savepmlq;

  kernel_t rkPmlUpdateKernel;
  kernel_t rkPmlStageKernel;

  void Backup(deviceMemory<dfloat> &o_Q);
  void Restore(deviceMemory<dfloat> &o_Q);
  void AcceptStep(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rq);

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt);

public:
  sark4_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields,
            memory<dfloat> _lambda, platform_t& _platform, comm_t _comm);
};

/* Semi-Analytic Explict Runge-Kutta, order 5 with embedded order 4 and adaptive time-stepping */
// Note: PML fields are not stepped Semi-Analytically
class sark5_pml: public sark5 {
private:
  dlong Npml;

  memory<dfloat> pmlrkA;
  deviceMemory<dfloat> o_pmlrkA;

  deviceMemory<dfloat> o_pmlq;
  deviceMemory<dfloat> o_rhspmlq;
  deviceMemory<dfloat> o_rkpmlq;
  deviceMemory<dfloat> o_rkrhspmlq;

  deviceMemory<dfloat> o_savepmlq;

  kernel_t rkPmlUpdateKernel;
  kernel_t rkPmlStageKernel;

  void Backup(deviceMemory<dfloat> &o_Q);
  void Restore(deviceMemory<dfloat> &o_Q);
  void AcceptStep(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rq);

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt);

public:
  sark5_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields,
            memory<dfloat> _lambda, platform_t& _platform, comm_t _comm);
};


/* Multi-rate Adams-Bashforth, order 3 */
class mrab3_pml: public mrab3 {
private:
  dlong Npml;
  int Npmlfields;

  deviceMemory<dfloat> o_pmlq;
  deviceMemory<dfloat> o_rhspmlq0, o_rhspmlq;

  kernel_t pmlUpdateKernel;

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

public:
  mrab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields, platform_t& _platform, mesh_t& _mesh);
};

/* Multi-rate Semi-Analytic Adams-Bashforth, order 3 */
// Note: PML fields are not stepped Semi-Analytically
class mrsaab3_pml: public mrsaab3 {
private:
  dlong Npml;
  int Npmlfields;

  deviceMemory<dfloat> o_pmlq;

  memory<dfloat> pmlsaab_a, pmlsaab_b;
  deviceMemory<dfloat> o_pmlsaab_a, o_pmlsaab_b;

  deviceMemory<dfloat> o_rhspmlq0, o_rhspmlq;

  kernel_t pmlUpdateKernel;

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

public:
  mrsaab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields,
            memory<dfloat> _lambda, platform_t& _platform, mesh_t& _mesh);
};

} //namespace TimeStepper

} //namespace libp

#endif
