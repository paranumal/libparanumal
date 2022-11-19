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

#include <optional>
#include "core.hpp"
#include "settings.hpp"
#include "mesh.hpp"
#include "solver.hpp"

namespace libp {

//forward declare
namespace TimeStepper {
  class timeStepperBase_t;
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

  void Run(solver_t& solver,
           deviceMemory<dfloat>& o_q,
           dfloat start, dfloat end);

  void RunWithPml(solver_t& solver,
                  deviceMemory<dfloat>& o_q,
                  deviceMemory<dfloat>& o_pmlq,
                  dfloat start, dfloat end);

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
  dlong Npml;

  dfloat dt;

  timeStepperBase_t(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
                    int Np, int Nfields, int Npmlfields,
                    platform_t& _platform, comm_t _comm):
    platform(_platform),
    comm(_comm),
    N(Nelements*Np*Nfields),
    Nhalo(NhaloElements*Np*Nfields),
    Npml(NpmlElements*Np*Npmlfields) {}

  virtual void Run(solver_t& solver,
                   deviceMemory<dfloat> o_q,
                   std::optional<deviceMemory<dfloat>> o_pmlq,
                   dfloat start, dfloat end)=0;

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

  kernel_t updateKernel;

  void Step(solver_t& solver,
            deviceMemory<dfloat> o_q,
            deviceMemory<dfloat> o_rhsq,
            std::optional<deviceMemory<dfloat>> o_pmlq,
            deviceMemory<dfloat> o_rhspmlq,
            dfloat time, dfloat _dt, int order);

public:
  ab3(dlong Nelements, dlong NhaloElements,
      int Np, int Nfields,
      platform_t& _platform, comm_t _comm);
  ab3(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
      int Np, int Nfields, int Npmlfields,
      platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver,
           deviceMemory<dfloat> o_q,
           std::optional<deviceMemory<dfloat>> o_pmlq,
           dfloat start, dfloat end);
};

/* Low-Storage Explicit Runge-Kutta, order 4 */
class lserk4: public timeStepperBase_t {
protected:
  int Nrk;
  memory<dfloat> rka, rkb, rkc;

  kernel_t updateKernel;

  void Step(solver_t& solver,
            deviceMemory<dfloat> o_q,
            std::optional<deviceMemory<dfloat>> o_pmlq,
            dfloat time, dfloat dt);

public:
  lserk4(dlong Nelements, dlong NhaloElements,
         int Np, int Nfields,
         platform_t& _platform, comm_t _comm);
  lserk4(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
         int Np, int Nfields, int Npmlfields,
         platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver,
           deviceMemory<dfloat> o_q,
           std::optional<deviceMemory<dfloat>> o_pmlq,
           dfloat start, dfloat end);
};

/* Dormand-Prince method */
/* Explict Runge-Kutta, order 5 with embedded order 4 and adaptive time-stepping */
class dopri5: public timeStepperBase_t {
protected:
  int Nrk;

  memory<dfloat> rkC, rkA, rkE;
  deviceMemory<dfloat> o_rkA, o_rkE;

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

  void Step(solver_t& solver,
            deviceMemory<dfloat> o_q,
            std::optional<deviceMemory<dfloat>> o_pmlq,
            deviceMemory<dfloat> o_rkq,
            deviceMemory<dfloat> o_rkpmlq,
            deviceMemory<dfloat> o_rkerr,
            dfloat time, dfloat _dt);

  dfloat Estimater(deviceMemory<dfloat>& o_q,
                   deviceMemory<dfloat>& o_rkq,
                   deviceMemory<dfloat>& o_rkerr);

public:
  dopri5(dlong Nelements, dlong NhaloElements,
         int Np, int Nfields,
         platform_t& _platform, comm_t _comm);
  dopri5(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
         int Np, int Nfields, int Npmlfields,
         platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver,
           deviceMemory<dfloat> o_q,
           std::optional<deviceMemory<dfloat>> o_pmlq,
           dfloat start, dfloat end);
};

/* Semi-Analytic Adams-Bashforth, order 3 */
class saab3: public timeStepperBase_t {
protected:
  static constexpr int Nstages{3};
  int shiftIndex;

  int Np, Nfields;
  dlong Nelements, NhaloElements;

  memory<dfloat> lambda;

  pinnedMemory<dfloat> h_saab_x, h_saab_a;
  deviceMemory<dfloat> o_saab_x, o_saab_a;

  memory<dfloat> pmlsaab_x, pmlsaab_a;
  deviceMemory<dfloat> o_pmlsaab_x, o_pmlsaab_a;

  kernel_t updateKernel;
  kernel_t pmlUpdateKernel;

  void Step(solver_t& solver,
            deviceMemory<dfloat> o_q,
            deviceMemory<dfloat> o_rhsq,
            std::optional<deviceMemory<dfloat>> o_pmlq,
            deviceMemory<dfloat> o_rhspmlq,
            dfloat time, dfloat _dt, int order);

  void UpdateCoefficients();

public:
  saab3(dlong _Nelements, dlong _NhaloElements,
        int _Np, int _Nfields,
        memory<dfloat> _lambda,
        platform_t& _platform, comm_t _comm);
  saab3(dlong _Nelements, dlong NpmlElements, dlong _NhaloElements,
        int _Np, int _Nfields, int Npmlfields,
        memory<dfloat> _lambda,
        platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver,
           deviceMemory<dfloat> o_q,
           std::optional<deviceMemory<dfloat>> o_pmlq,
           dfloat start, dfloat end);
};

/* Semi-Analytic Explict Runge-Kutta, order 4 with embedded order 3 and adaptive time-stepping */
class sark4: public timeStepperBase_t {
protected:
  static constexpr int Nrk{5};
  int order, embeddedOrder;

  int Np, Nfields;
  dlong Nelements, NhaloElements;

  memory<dfloat> lambda;

  memory<dfloat> rkC;
  deviceMemory<dfloat> o_rkX, o_rkA, o_rkE;
  pinnedMemory<dfloat> h_rkX, h_rkA, h_rkE;
  memory<dfloat> pmlrkA;
  deviceMemory<dfloat> o_pmlrkA;

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


  void Step(solver_t& solver,
            deviceMemory<dfloat> o_q,
            std::optional<deviceMemory<dfloat>> o_pmlq,
            deviceMemory<dfloat> o_rkq,
            deviceMemory<dfloat> o_rkpmlq,
            deviceMemory<dfloat> o_rkerr,
            dfloat time, dfloat _dt);

  dfloat Estimater(deviceMemory<dfloat>& o_q,
                   deviceMemory<dfloat>& o_rkq,
                   deviceMemory<dfloat>& o_rkerr);

  void UpdateCoefficients();

public:
  sark4(dlong Nelements, dlong NhaloElements,
        int Np, int Nfields,
        memory<dfloat> _lambda,
        platform_t& _platform, comm_t _comm);
  sark4(dlong _Nelements, dlong NpmlElements, dlong _NhaloElements,
        int _Np, int _Nfields, int _Npmlfields,
        memory<dfloat> _lambda,
        platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver,
           deviceMemory<dfloat> o_q,
           std::optional<deviceMemory<dfloat>> o_pmlq,
           dfloat start, dfloat end);
};

/* Semi-Analytic Explict Runge-Kutta, order 5 with embedded order 4 and adaptive time-stepping */
class sark5: public timeStepperBase_t {
protected:
  static constexpr int Nrk{7};
  int order, embeddedOrder;

  int Np, Nfields;
  dlong Nelements, NhaloElements;

  memory<dfloat> lambda;

  memory<dfloat> rkC;
  deviceMemory<dfloat> o_rkX, o_rkA, o_rkE;
  pinnedMemory<dfloat> h_rkX, h_rkA, h_rkE;
  memory<dfloat> pmlrkA;
  deviceMemory<dfloat> o_pmlrkA;

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

   void Step(solver_t& solver,
             deviceMemory<dfloat> o_q,
             std::optional<deviceMemory<dfloat>> o_pmlq,
             deviceMemory<dfloat> o_rkq,
             deviceMemory<dfloat> o_rkpmlq,
             deviceMemory<dfloat> o_rkerr,
             dfloat time, dfloat _dt);

  dfloat Estimater(deviceMemory<dfloat>& o_q,
                   deviceMemory<dfloat>& o_rkq,
                   deviceMemory<dfloat>& o_rkerr);

  void UpdateCoefficients();

public:
  sark5(dlong Nelements, dlong NhaloElements,
        int Np, int Nfields,
        memory<dfloat> _lambda,
        platform_t& _platform, comm_t _comm);
  sark5(dlong _Nelements, dlong NpmlElements, dlong _NhaloElements,
        int _Np, int _Nfields, int Npmlfields,
        memory<dfloat> _lambda,
        platform_t& _platform, comm_t _comm);

  void Run(solver_t& solver,
           deviceMemory<dfloat> o_q,
           std::optional<deviceMemory<dfloat>> o_pmlq,
           dfloat start, dfloat end);
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

  kernel_t rhsKernel;

  void Step(solver_t& solver,
            deviceMemory<dfloat> o_q,
            deviceMemory<dfloat> o_qn,
            deviceMemory<dfloat> o_F,
            dfloat time, dfloat _dt, int order);

public:
  extbdf3(dlong Nelements, dlong NhaloElements,
          int Np, int Nfields,
          platform_t& _platform, comm_t _comm);

  dfloat GetGamma();

  void Run(solver_t& solver,
           deviceMemory<dfloat> o_q,
           std::optional<deviceMemory<dfloat>> o_pmlq,
           dfloat start, dfloat end);
};

/* Backward Difference Formula, order 3, with subcycling */
class ssbdf3: public timeStepperBase_t {
protected:
  static constexpr int Nstages{3};
  int shiftIndex;

  memory<dfloat> ssbdf_b;
  deviceMemory<dfloat> o_ssbdf_b;

  kernel_t rhsKernel;

  void Step(solver_t& solver,
            deviceMemory<dfloat> o_q,
            deviceMemory<dfloat> o_qn,
            dfloat time, dfloat _dt, int order);

public:
  ssbdf3(dlong Nelements, dlong NhaloElements,
      int Np, int Nfields,
      platform_t& _platform, comm_t _comm);

  dfloat GetGamma();

  void Run(solver_t& solver,
           deviceMemory<dfloat> o_q,
           std::optional<deviceMemory<dfloat>> o_pmlq,
           dfloat start, dfloat end);
};

/* Multi-rate Adams-Bashforth, order 3 */
class mrab3: public timeStepperBase_t {
protected:
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

  kernel_t updateKernel;
  kernel_t pmlUpdateKernel;
  kernel_t traceUpdateKernel;

  void Step(solver_t& solver,
            deviceMemory<dfloat> o_q,
            std::optional<deviceMemory<dfloat>> o_pmlq,
            deviceMemory<dfloat> o_rhsq0,
            deviceMemory<dfloat> o_rhsq,
            deviceMemory<dfloat> o_rhspmlq0,
            deviceMemory<dfloat> o_rhspmlq,
            deviceMemory<dfloat> o_fQM,
            dfloat time, dfloat _dt, int order);

public:
  mrab3(dlong Nelements, dlong NhaloElements,
        int _Np, int _Nfields,
        platform_t& _platform, mesh_t& _mesh);
  mrab3(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
        int _Np, int _Nfields, int _Npmlfields,
        platform_t& _platform, mesh_t& _mesh);

  void Run(solver_t& solver,
           deviceMemory<dfloat> o_q,
           std::optional<deviceMemory<dfloat>> o_pmlq,
           dfloat start, dfloat end);
};

/* Multi-rate Semi-Analytic Adams-Bashforth, order 3 */
class mrsaab3: public timeStepperBase_t {
protected:
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

  kernel_t updateKernel;
  kernel_t pmlUpdateKernel;
  kernel_t traceUpdateKernel;

  void Step(solver_t& solver,
            deviceMemory<dfloat> o_q,
            std::optional<deviceMemory<dfloat>> o_pmlq,
            deviceMemory<dfloat> o_rhsq0,
            deviceMemory<dfloat> o_rhsq,
            deviceMemory<dfloat> o_rhspmlq0,
            deviceMemory<dfloat> o_rhspmlq,
            deviceMemory<dfloat> o_fQM,
            dfloat time, dfloat _dt, int order);

  void UpdateCoefficients();

public:
  mrsaab3(dlong _Nelements, dlong _NhaloElements,
          int _Np, int _Nfields,
          memory<dfloat> _lambda,
          platform_t& _platform, mesh_t& _mesh);
  mrsaab3(dlong _Nelements, dlong NpmlElements, dlong _NhaloElements,
          int _Np, int _Nfields, int _Npmlfields,
          memory<dfloat> _lambda,
          platform_t& _platform, mesh_t& _mesh);

  void Run(solver_t& solver,
           deviceMemory<dfloat> o_q,
           std::optional<deviceMemory<dfloat>> o_pmlq,
           dfloat start, dfloat end);
};

} //namespace TimeStepper

} //namespace libp

#endif
