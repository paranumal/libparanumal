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
namespace TimeStepper {
  class timeStepperBase_t;
  class pmlTimeStepperBase_t;
}

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

/* General PML TimeStepper object*/
class pmlTimeStepper_t {
 public:
  pmlTimeStepper_t() = default;

  /*Generic setup. Create a Stepper object and wrap it in a shared_ptr*/
  template<class Stepper, class... Args>
  void Setup(Args&& ... args) {
    ts = std::make_shared<Stepper>(args...);
  }

  void Run(solver_t& solver,
           deviceMemory<dfloat>& o_q,
           deviceMemory<dfloat>& o_pmlq,
           dfloat start, dfloat end);

  void SetTimeStep(dfloat dt_);

  dfloat GetTimeStep();

  dfloat GetGamma();

 private:
  std::shared_ptr<TimeStepper::pmlTimeStepperBase_t> ts=nullptr;

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
  static constexpr int Nstages{3};
  int shiftIndex;

  memory<dfloat> ab_a;
  deviceMemory<dfloat> o_ab_a;

  deviceMemory<dfloat> o_rhsq;

  kernel_t updateKernel;

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

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

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt);

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

  void Backup(deviceMemory<dfloat> &o_Q);
  void Restore(deviceMemory<dfloat> &o_Q);
  void AcceptStep(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rq);

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt);

  dfloat Estimater(deviceMemory<dfloat>& o_q);

public:
  dopri5(dlong Nelements, dlong NhaloElements,
         int Np, int Nfields,
         platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat start, dfloat end);
};

/* Semi-Analytic Adams-Bashforth, order 3 */
class saab3: public timeStepperBase_t {
protected:
  static constexpr int Nstages{3};
  int shiftIndex;

  int Np, Nfields;
  dlong Nblock, Nelements, NhaloElements;

  memory<dfloat> lambda;

  pinnedMemory<dfloat> h_saab_x, h_saab_a;
  deviceMemory<dfloat> o_saab_x, o_saab_a;

  deviceMemory<dfloat> o_rhsq;

  kernel_t updateKernel;

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

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
  static constexpr int Nrk{5};
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

  void Backup(deviceMemory<dfloat> &o_Q);
  void Restore(deviceMemory<dfloat> &o_Q);
  void AcceptStep(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rq);

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt);

  dfloat Estimater(deviceMemory<dfloat>& o_q);

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
  static constexpr int Nrk{7};
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

  void Backup(deviceMemory<dfloat> &o_Q);
  void Restore(deviceMemory<dfloat> &o_Q);
  void AcceptStep(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rq);

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt);

  dfloat Estimater(deviceMemory<dfloat>& o_q);

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
  static constexpr int Nstages{3};
  int shiftIndex;

  memory<dfloat> extbdf_a;
  memory<dfloat> extbdf_b;
  deviceMemory<dfloat> o_extbdf_a;
  deviceMemory<dfloat> o_extbdf_b;

  deviceMemory<dfloat> o_rhs;
  deviceMemory<dfloat> o_qn;
  deviceMemory<dfloat> o_F;

  kernel_t rhsKernel;

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

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
  static constexpr int Nstages{3};
  int shiftIndex;

  memory<dfloat> ssbdf_b;
  deviceMemory<dfloat> o_ssbdf_b;

  deviceMemory<dfloat> o_rhs;
  deviceMemory<dfloat> o_qn;
  deviceMemory<dfloat> o_qhat;

  kernel_t rhsKernel;

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

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

  static constexpr int Nstages{3};
  int Nlevels;
  int Nfields;

  deviceMemory<int> o_shiftIndex;
  pinnedMemory<int> h_shiftIndex;

  memory<dfloat> mrdt;
  deviceMemory<dfloat> o_mrdt;
  deviceMemory<dfloat> o_zeros;

  memory<dfloat> ab_a, ab_b;
  deviceMemory<dfloat> o_ab_a, o_ab_b;

  deviceMemory<dfloat> o_rhsq0, o_rhsq, o_fQM;

  kernel_t updateKernel;
  kernel_t traceUpdateKernel;

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

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

  static constexpr int Nstages{3};
  int Nlevels;
  int Nfields;

  memory<dfloat> lambda;

  deviceMemory<int> o_shiftIndex;
  pinnedMemory<int> h_shiftIndex;

  memory<dfloat> mrdt;
  deviceMemory<dfloat> o_mrdt;
  deviceMemory<dfloat> o_zeros;

  pinnedMemory<dfloat> h_saab_x, h_saab_a, h_saab_b;
  deviceMemory<dfloat> o_saab_x, o_saab_a, o_saab_b;

  deviceMemory<dfloat> o_rhsq0, o_rhsq, o_fQM;

  kernel_t updateKernel;
  kernel_t traceUpdateKernel;

  void Step(solver_t& solver, deviceMemory<dfloat>& o_q, dfloat time, dfloat dt, int order);

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

//base pml time stepper class
class pmlTimeStepperBase_t {
public:
  platform_t platform;
  comm_t comm;

  dlong N;
  dlong Nhalo;
  dlong Npml;

  dfloat dt;

  pmlTimeStepperBase_t(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
                    int Np, int Nfields, int Npmlfields,
                    platform_t& _platform, comm_t _comm):
    platform(_platform),
    comm(_comm),
    N(Nelements*Np*Nfields),
    Nhalo(NhaloElements*Np*Nfields),
    Npml(NpmlElements*Np*Npmlfields) {}

  virtual void Run(solver_t& solver,
                   deviceMemory<dfloat>& o_q,
                   deviceMemory<dfloat>& o_pmlq,
                   dfloat start, dfloat end) = 0;

  void SetTimeStep(dfloat dt_) {dt = dt_;};

  dfloat GetTimeStep() {return dt;};

  virtual dfloat GetGamma() {
    LIBP_FORCE_ABORT("GetGamma() not available in this Timestepper");
    return 0.0;
  }
};

/* Adams Bashforth, order 3 */
class ab3_pml: public pmlTimeStepperBase_t {
private:
  static constexpr int Nstages{3};
  int shiftIndex;

  memory<dfloat> ab_a;
  deviceMemory<dfloat> o_ab_a;

  deviceMemory<dfloat> o_rhsq;
  deviceMemory<dfloat> o_rhspmlq;

  kernel_t updateKernel;

  void Step(solver_t& solver,
            deviceMemory<dfloat>& o_q,
            deviceMemory<dfloat>& o_pmlq,
            dfloat time, dfloat dt, int order);

public:
  ab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
          int Np, int Nfields, int Npmlfields,
          platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver,
           deviceMemory<dfloat>& o_q,
           deviceMemory<dfloat>& o_pmlq,
           dfloat start, dfloat end);
};

/* Low-Storage Explicit Runge-Kutta, order 4 */
class lserk4_pml: public pmlTimeStepperBase_t {
private:
  int Nrk;
  memory<dfloat> rka, rkb, rkc;

  deviceMemory<dfloat> o_rhsq;
  deviceMemory<dfloat> o_resq;
  deviceMemory<dfloat> o_rhspmlq;
  deviceMemory<dfloat> o_respmlq;

  deviceMemory<dfloat> o_saveq;
  deviceMemory<dfloat> o_savepmlq;

  kernel_t updateKernel;

  void Step(solver_t& solver,
            deviceMemory<dfloat>& o_q,
            deviceMemory<dfloat>& o_pmlq,
            dfloat time, dfloat dt);

public:
  lserk4_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int Npmlfields,
            platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver,
           deviceMemory<dfloat>& o_q,
           deviceMemory<dfloat>& o_pmlq,
           dfloat start, dfloat end);
};

/* Dormand-Prince method */
/* Explict Runge-Kutta, order 5 with embedded order 4 and adaptive time-stepping */
class dopri5_pml: public pmlTimeStepperBase_t {
private:
  int Nrk;

  dlong Nblock;

  memory<dfloat> rkC, rkA, rkE;
  deviceMemory<dfloat> o_rkA, o_rkE;

  deviceMemory<dfloat> o_errtmp;
  pinnedMemory<dfloat> h_errtmp;

  deviceMemory<dfloat> o_rhsq;
  deviceMemory<dfloat> o_rkq;
  deviceMemory<dfloat> o_rkrhsq;
  deviceMemory<dfloat> o_rhspmlq;
  deviceMemory<dfloat> o_rkpmlq;
  deviceMemory<dfloat> o_rkrhspmlq;

  deviceMemory<dfloat> o_rkerr;

  deviceMemory<dfloat> o_saveq;
  deviceMemory<dfloat> o_savepmlq;


  kernel_t rkUpdateKernel;
  kernel_t rkPmlUpdateKernel;
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

  void Backup(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_pmlq);
  void Restore(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_pmlq);
  void AcceptStep(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rq,
                  deviceMemory<dfloat> &o_pmlq, deviceMemory<dfloat> &o_rpmlq);

  void Step(solver_t& solver,
            deviceMemory<dfloat>& o_q,
            deviceMemory<dfloat>& o_pmlq,
            dfloat time, dfloat dt);

  dfloat Estimater(deviceMemory<dfloat>& o_q);

public:
  dopri5_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int Npmlfields,
            platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver,
           deviceMemory<dfloat>& o_q,
           deviceMemory<dfloat>& o_pmlq,
           dfloat start, dfloat end);
};

/* Semi-Analytic Adams-Bashforth, order 3 */
// Note: PML fields are not stepped Semi-Analytically
class saab3_pml: public pmlTimeStepperBase_t {
private:
  static constexpr int Nstages{3};
  int shiftIndex;

  int Np, Nfields;
  dlong Nblock, Nelements, NhaloElements;

  memory<dfloat> lambda;

  pinnedMemory<dfloat> h_saab_x, h_saab_a;
  deviceMemory<dfloat> o_saab_x, o_saab_a;

  memory<dfloat> pmlsaab_x, pmlsaab_a;
  deviceMemory<dfloat> o_pmlsaab_x, o_pmlsaab_a;

  deviceMemory<dfloat> o_rhsq;
  deviceMemory<dfloat> o_rhspmlq;

  kernel_t updateKernel;
  kernel_t pmlUpdateKernel;

  void Step(solver_t& solver,
            deviceMemory<dfloat>& o_q,
            deviceMemory<dfloat>& o_pmlq,
            dfloat time, dfloat dt, int order);

public:
  saab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields,
            memory<dfloat> _lambda,
            platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver,
           deviceMemory<dfloat>& o_q,
           deviceMemory<dfloat>& o_pmlq,
           dfloat start, dfloat end);
};

/* Semi-Analytic Explict Runge-Kutta, order 4 with embedded order 3 and adaptive time-stepping */
// Note: PML fields are not stepped Semi-Analytically
class sark4_pml: public pmlTimeStepperBase_t {
private:
  static constexpr int Nrk{5};
  int order, embeddedOrder;

  int Np, Nfields;
  dlong Nblock, Nelements, NhaloElements;

  memory<dfloat> lambda;

  memory<dfloat> rkC;
  deviceMemory<dfloat> o_rkX, o_rkA, o_rkE;
  pinnedMemory<dfloat> h_rkX, h_rkA, h_rkE;
  memory<dfloat> pmlrkA;
  deviceMemory<dfloat> o_pmlrkA;

  deviceMemory<dfloat> o_rhsq;
  deviceMemory<dfloat> o_rkq;
  deviceMemory<dfloat> o_rkrhsq;
  deviceMemory<dfloat> o_rhspmlq;
  deviceMemory<dfloat> o_rkpmlq;
  deviceMemory<dfloat> o_rkrhspmlq;

  deviceMemory<dfloat> o_rkerr;

  deviceMemory<dfloat> o_saveq;
  deviceMemory<dfloat> o_savepmlq;

  deviceMemory<dfloat> o_errtmp;
  pinnedMemory<dfloat> h_errtmp;

  kernel_t rkUpdateKernel;
  kernel_t rkPmlUpdateKernel;
  kernel_t rkStageKernel;
  kernel_t rkPmlStageKernel;
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

  void Backup(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_pmlq);
  void Restore(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_pmlq);
  void AcceptStep(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rq,
                  deviceMemory<dfloat> &o_pmlq, deviceMemory<dfloat> &o_rpmlq);

  void Step(solver_t& solver,
            deviceMemory<dfloat>& o_q,
            deviceMemory<dfloat>& o_pmlq,
            dfloat time, dfloat dt);

  dfloat Estimater(deviceMemory<dfloat>& o_q);

public:
  sark4_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields,
            memory<dfloat> _lambda, platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver,
           deviceMemory<dfloat>& o_q,
           deviceMemory<dfloat>& o_pmlq,
           dfloat start, dfloat end);
};

/* Semi-Analytic Explict Runge-Kutta, order 5 with embedded order 4 and adaptive time-stepping */
// Note: PML fields are not stepped Semi-Analytically
class sark5_pml: public pmlTimeStepperBase_t {
private:
  static constexpr int Nrk{7};
  int order, embeddedOrder;

  int Np, Nfields;
  dlong Nblock, Nelements, NhaloElements;

  memory<dfloat> lambda;

  memory<dfloat> rkC;
  deviceMemory<dfloat> o_rkX, o_rkA, o_rkE;
  pinnedMemory<dfloat> h_rkX, h_rkA, h_rkE;
  memory<dfloat> pmlrkA;
  deviceMemory<dfloat> o_pmlrkA;


  deviceMemory<dfloat> o_rhsq;
  deviceMemory<dfloat> o_rkq;
  deviceMemory<dfloat> o_rkrhsq;
  deviceMemory<dfloat> o_rhspmlq;
  deviceMemory<dfloat> o_rkpmlq;
  deviceMemory<dfloat> o_rkrhspmlq;

  deviceMemory<dfloat> o_rkerr;

  deviceMemory<dfloat> o_saveq;
  deviceMemory<dfloat> o_savepmlq;

  deviceMemory<dfloat> o_errtmp;
  pinnedMemory<dfloat> h_errtmp;

  kernel_t rkUpdateKernel;
  kernel_t rkPmlUpdateKernel;
  kernel_t rkStageKernel;
  kernel_t rkPmlStageKernel;
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

  void Backup(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_pmlq);
  void Restore(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_pmlq);
  void AcceptStep(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rq,
                  deviceMemory<dfloat> &o_pmlq, deviceMemory<dfloat> &o_rpmlq);

  void Step(solver_t& solver,
            deviceMemory<dfloat>& o_q,
            deviceMemory<dfloat>& o_pmlq,
            dfloat time, dfloat dt);

  dfloat Estimater(deviceMemory<dfloat>& o_q);

public:
  sark5_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields,
            memory<dfloat> _lambda, platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver,
           deviceMemory<dfloat>& o_q,
           deviceMemory<dfloat>& o_pmlq,
           dfloat start, dfloat end);
};

/* Multi-rate Adams-Bashforth, order 3 */
class mrab3_pml: public pmlTimeStepperBase_t {
private:
  mesh_t mesh;

  static constexpr int Nstages{3};
  int Nlevels;
  int Nfields;
  int Npmlfields;

  deviceMemory<int> o_shiftIndex;
  pinnedMemory<int> h_shiftIndex;

  memory<dfloat> mrdt;
  deviceMemory<dfloat> o_mrdt;
  deviceMemory<dfloat> o_zeros;

  memory<dfloat> ab_a, ab_b;
  deviceMemory<dfloat> o_ab_a, o_ab_b;

  deviceMemory<dfloat> o_rhsq0, o_rhsq, o_fQM;
  deviceMemory<dfloat> o_rhspmlq0, o_rhspmlq;

  kernel_t updateKernel;
  kernel_t pmlUpdateKernel;
  kernel_t traceUpdateKernel;

  void Step(solver_t& solver,
            deviceMemory<dfloat>& o_q,
            deviceMemory<dfloat>& o_pmlq,
            dfloat time, dfloat dt, int order);

public:
  mrab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields, platform_t& _platform, mesh_t& _mesh);

  void Run(solver_t& solver,
           deviceMemory<dfloat>& o_q,
           deviceMemory<dfloat>& o_pmlq,
           dfloat start, dfloat end);
};

/* Multi-rate Semi-Analytic Adams-Bashforth, order 3 */
// Note: PML fields are not stepped Semi-Analytically
class mrsaab3_pml: public pmlTimeStepperBase_t {
private:
  mesh_t mesh;

  static constexpr int Nstages{3};
  int Nlevels;
  int Nfields;
  int Npmlfields;

  memory<dfloat> lambda;

  deviceMemory<int> o_shiftIndex;
  pinnedMemory<int> h_shiftIndex;

  memory<dfloat> mrdt;
  deviceMemory<dfloat> o_mrdt;
  deviceMemory<dfloat> o_zeros;

  pinnedMemory<dfloat> h_saab_x, h_saab_a, h_saab_b;
  deviceMemory<dfloat> o_saab_x, o_saab_a, o_saab_b;
  memory<dfloat> pmlsaab_a, pmlsaab_b;
  deviceMemory<dfloat> o_pmlsaab_a, o_pmlsaab_b;

  deviceMemory<dfloat> o_rhsq0, o_rhsq, o_fQM;
  deviceMemory<dfloat> o_rhspmlq0, o_rhspmlq;

  kernel_t updateKernel;
  kernel_t pmlUpdateKernel;
  kernel_t traceUpdateKernel;

  void Step(solver_t& solver,
            deviceMemory<dfloat>& o_q,
            deviceMemory<dfloat>& o_pmlq,
            dfloat time, dfloat dt, int order);

public:
  mrsaab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields,
            memory<dfloat> _lambda, platform_t& _platform, mesh_t& _mesh);

  void Run(solver_t& solver,
           deviceMemory<dfloat>& o_q,
           deviceMemory<dfloat>& o_pmlq,
           dfloat start, dfloat end);
};

} //namespace TimeStepper

} //namespace libp

#endif
