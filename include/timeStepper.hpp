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

#ifndef TIMESTEPPER_HPP
#define TIMESTEPPER_HPP

#include <occa.hpp>
#include "types.h"
#include "utils.hpp"
#include "solver.hpp"
#include "settings.hpp"

namespace TimeStepper {

//base time stepper class
class timeStepper_t {
public:
  dlong N;
  dlong Nhalo;

  solver_t& solver;

  MPI_Comm& comm;
  occa::device& device;
  settings_t& settings;
  occa::properties& props;

  dfloat dt;

  timeStepper_t(dlong Nelements, dlong NhaloElements,
                 int Np, int Nfields, solver_t& _solver):
    N(Nelements*Np*Nfields),
    Nhalo(NhaloElements*Np*Nfields),
    solver(_solver),
    comm(_solver.comm),
    device(_solver.device),
    settings(_solver.settings),
    props(_solver.props) {}

  virtual ~timeStepper_t() {};
  virtual void Run(occa::memory& o_q, dfloat start, dfloat end)=0;

  void SetTimeStep(dfloat dt_) {dt = dt_;};
};

/* Adams Bashforth, order 3 */
class ab3: public timeStepper_t {
protected:
  int Nstages;
  int shiftIndex;

  dfloat *ab_a;
  occa::memory o_ab_a;

  occa::memory o_rhsq;

  occa::kernel updateKernel;

  virtual void Step(occa::memory& o_q, dfloat time, dfloat dt, int order);

public:
  ab3(dlong Nelements, dlong NhaloElements,
      int Np, int Nfields, solver_t& solver);
  ~ab3();

  void Run(occa::memory& o_q, dfloat start, dfloat end);
};

/* Low-Storage Explicit Runge-Kutta, order 4 */
class lserk4: public timeStepper_t {
protected:
  int Nrk;
  dfloat *rka, *rkb, *rkc;

  occa::memory o_rhsq;
  occa::memory o_resq;

  occa::kernel updateKernel;

  virtual void Step(occa::memory& o_q, dfloat time, dfloat dt);

public:
  lserk4(dlong Nelements, dlong NhaloElements,
         int Np, int Nfields, solver_t& solver);
  ~lserk4();

  void Run(occa::memory& o_q, dfloat start, dfloat end);
};

/* Dormand-Prince method */
/* Explict Runge-Kutta, order 5 with embedded order 4 and adaptive time-stepping */
class dopri5: public timeStepper_t {
protected:
  int Nrk;

  dlong Nblock;

  dfloat *rkC, *rkA, *rkE;
  occa::memory o_rkA, o_rkE;

  dfloat *errtmp;
  occa::memory o_errtmp, h_errtmp;

  occa::memory o_rhsq;
  occa::memory o_rkq;
  occa::memory o_rkrhsq;
  occa::memory o_rkerr;

  occa::memory o_saveq;


  occa::kernel rkUpdateKernel;
  occa::kernel rkStageKernel;
  occa::kernel rkErrorEstimateKernel;

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

  virtual void Backup(occa::memory &o_Q);
  virtual void Restore(occa::memory &o_Q);
  virtual void AcceptStep(occa::memory &o_q, occa::memory &o_rq);

  virtual void Step(occa::memory& o_q, dfloat time, dfloat dt);

  virtual dfloat Estimater(occa::memory& o_q);

public:
  dopri5(dlong Nelements, dlong NhaloElements,
         int Np, int Nfields, solver_t& solver);
  ~dopri5();

  void Run(occa::memory& o_q, dfloat start, dfloat end);
};

/* Semi-Analytic Adams-Bashforth, order 3 */
class saab3: public timeStepper_t {
protected:
  int Nstages;
  int shiftIndex;

  int Np, Nfields;
  dlong Nblock, Nelements, NhaloElements;

  dfloat *lambda;

  dfloat *saab_x, *saab_a;
  occa::memory o_saab_x, o_saab_a;

  occa::memory o_rhsq;

  occa::kernel updateKernel;

  virtual void Step(occa::memory& o_q, dfloat time, dfloat dt, int order);

  virtual void UpdateCoefficients();

public:
  saab3(dlong _Nelements, dlong _NhaloElements,
        int _Np, int _Nfields,
        dfloat *_lambda,
        solver_t& _solver);
  ~saab3();

  void Run(occa::memory& o_q, dfloat start, dfloat end);
};

/* Semi-Analytic Explict Runge-Kutta, order 4 with embedded order 3 and adaptive time-stepping */
class sark4: public timeStepper_t {
protected:
  int Nrk;
  int order, embeddedOrder;

  int Np, Nfields;
  dlong Nblock, Nelements, NhaloElements;

  dfloat *lambda;

  dfloat *rkC, *rkX, *rkA, *rkE;
  occa::memory h_rkX, h_rkA, h_rkE;
  occa::memory o_rkX, o_rkA, o_rkE;

  dfloat *errtmp;

  occa::memory o_rhsq;
  occa::memory o_rkq;
  occa::memory o_rkrhsq;
  occa::memory o_rkerr;

  occa::memory o_saveq;

  occa::memory o_errtmp;

  occa::kernel rkUpdateKernel;
  occa::kernel rkStageKernel;
  occa::kernel rkErrorEstimateKernel;

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

  virtual void Backup(occa::memory &o_Q);
  virtual void Restore(occa::memory &o_Q);
  virtual void AcceptStep(occa::memory &o_q, occa::memory &o_rq);

  virtual void Step(occa::memory& o_q, dfloat time, dfloat dt);

  dfloat Estimater(occa::memory& o_q);

  void UpdateCoefficients();

public:
  sark4(dlong _Nelements, dlong _NhaloElements,
        int _Np, int _Nfields,
        dfloat *_lambda,
        solver_t& _solver);
  ~sark4();

  void Run(occa::memory& o_q, dfloat start, dfloat end);
};

/* Semi-Analytic Explict Runge-Kutta, order 5 with embedded order 4 and adaptive time-stepping */
class sark5: public timeStepper_t {
protected:
  int Nrk;
  int order, embeddedOrder;

  int Np, Nfields;
  dlong Nblock, Nelements, NhaloElements;

  dfloat *lambda;

  dfloat *rkC, *rkX, *rkA, *rkE;
  occa::memory h_rkX, h_rkA, h_rkE;
  occa::memory o_rkX, o_rkA, o_rkE;

  dfloat *errtmp;

  occa::memory o_rhsq;
  occa::memory o_rkq;
  occa::memory o_rkrhsq;
  occa::memory o_rkerr;

  occa::memory o_saveq;

  occa::memory o_errtmp;

  occa::kernel rkUpdateKernel;
  occa::kernel rkStageKernel;
  occa::kernel rkErrorEstimateKernel;

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

  virtual void Backup(occa::memory &o_Q);
  virtual void Restore(occa::memory &o_Q);
  virtual void AcceptStep(occa::memory &o_q, occa::memory &o_rq);

  virtual void Step(occa::memory& o_q, dfloat time, dfloat dt);

  dfloat Estimater(occa::memory& o_q);

  void UpdateCoefficients();

public:
  sark5(dlong _Nelements, dlong _NhaloElements,
        int _Np, int _Nfields,
        dfloat *_lambda,
        solver_t& _solver);
  ~sark5();

  void Run(occa::memory& o_q, dfloat start, dfloat end);
};


/* Multi-rate Adams-Bashforth, order 3 */
class mrab3: public timeStepper_t {
protected:
  mesh_t &mesh;

  int Nstages;
  int Nlevels;
  int Nfields;

  int* shiftIndex;
  occa::memory o_shiftIndex, h_shiftIndex;

  dfloat *mrdt;
  occa::memory o_mrdt;

  dfloat *ab_a, *ab_b;
  occa::memory o_ab_a, o_ab_b;

  occa::memory o_rhsq0, o_rhsq, o_fQM;

  occa::kernel updateKernel;
  occa::kernel traceUpdateKernel;

  virtual void Step(occa::memory& o_q, dfloat time, dfloat dt, int order);

public:
  mrab3(dlong _Nelements, dlong _NhaloElements,
         int _Np, int _Nfields,
         solver_t& _solver);
  ~mrab3();

  void Run(occa::memory& o_q, dfloat start, dfloat end);
};

/* Multi-rate Semi-Analytic Adams-Bashforth, order 3 */
class mrsaab3: public timeStepper_t {
protected:
  mesh_t &mesh;

  int Nstages;
  int Nlevels;
  int Nfields;

  dfloat *lambda;

  int* shiftIndex;
  occa::memory o_shiftIndex, h_shiftIndex;

  dfloat *mrdt;
  occa::memory o_mrdt;

  dfloat *saab_x, *saab_a, *saab_b;
  occa::memory o_saab_x, o_saab_a, o_saab_b;

  occa::memory o_rhsq0, o_rhsq, o_fQM;

  occa::kernel updateKernel;
  occa::kernel traceUpdateKernel;

  virtual void Step(occa::memory& o_q, dfloat time, dfloat dt, int order);

  void UpdateCoefficients();

public:
  mrsaab3(dlong _Nelements, dlong _NhaloElements,
         int _Np, int _Nfields,
         dfloat *_lambda,
         solver_t& _solver);
  ~mrsaab3();

  void Init();
  void Run(occa::memory& o_q, dfloat start, dfloat end);
};


/**************************************************/
/* Derived Time Integrators which step PML fields */
/**************************************************/

/* Adams Bashforth, order 3 */
class ab3_pml: public ab3 {
private:
  dlong Npml;

  occa::memory o_pmlq;
  occa::memory o_rhspmlq;

  void Step(occa::memory& o_q, dfloat time, dfloat dt, int order);

public:
  ab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
          int Np, int Nfields, int Npmlfields, solver_t& solver);
  ~ab3_pml();
};

/* Low-Storage Explicit Runge-Kutta, order 4 */
class lserk4_pml: public lserk4 {
private:
  dlong Npml;

  occa::memory o_pmlq;
  occa::memory o_rhspmlq;
  occa::memory o_respmlq;

  void Step(occa::memory& o_q, dfloat time, dfloat dt);

public:
  lserk4_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int Npmlfields, solver_t& solver);
  ~lserk4_pml();
};

/* Dormand-Prince method */
/* Explict Runge-Kutta, order 5 with embedded order 4 and adaptive time-stepping */
class dopri5_pml: public dopri5 {
private:
  dlong Npml;

  occa::memory o_pmlq;
  occa::memory o_rhspmlq;
  occa::memory o_rkpmlq;
  occa::memory o_rkrhspmlq;

  occa::memory o_savepmlq;

  occa::kernel rkPmlUpdateKernel;

  void Backup(occa::memory &o_Q);
  void Restore(occa::memory &o_Q);
  void AcceptStep(occa::memory &o_q, occa::memory &o_rq);

  void Step(occa::memory& o_q, dfloat time, dfloat dt);

public:
  dopri5_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int Npmlfields, solver_t& solver);
  ~dopri5_pml();
};

/* Semi-Analytic Adams-Bashforth, order 3 */
// Note: PML fields are not stepped Semi-Analytically
class saab3_pml: public saab3 {
private:
  dlong Npml;

  dfloat *pmlsaab_x, *pmlsaab_a;
  occa::memory o_pmlsaab_x, o_pmlsaab_a;

  occa::memory o_pmlq;
  occa::memory o_rhspmlq;

  occa::kernel pmlUpdateKernel;

  void Step(occa::memory& o_q, dfloat time, dfloat dt, int order);

public:
  saab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields,
            dfloat *_lambda, solver_t& solver);
  ~saab3_pml();
};

/* Semi-Analytic Explict Runge-Kutta, order 4 with embedded order 3 and adaptive time-stepping */
// Note: PML fields are not stepped Semi-Analytically
class sark4_pml: public sark4 {
private:
  dlong Npml;

  dfloat *pmlrkA;
  occa::memory o_pmlrkA;

  occa::memory o_pmlq;
  occa::memory o_rhspmlq;
  occa::memory o_rkpmlq;
  occa::memory o_rkrhspmlq;

  occa::memory o_savepmlq;

  occa::kernel rkPmlUpdateKernel;
  occa::kernel rkPmlStageKernel;

  void Backup(occa::memory &o_Q);
  void Restore(occa::memory &o_Q);
  void AcceptStep(occa::memory &o_q, occa::memory &o_rq);

  void Step(occa::memory& o_q, dfloat time, dfloat dt);

public:
  sark4_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields,
            dfloat *_lambda, solver_t& solver);
  ~sark4_pml();
};

/* Semi-Analytic Explict Runge-Kutta, order 5 with embedded order 4 and adaptive time-stepping */
// Note: PML fields are not stepped Semi-Analytically
class sark5_pml: public sark5 {
private:
  dlong Npml;

  dfloat *pmlrkA;
  occa::memory o_pmlrkA;

  occa::memory o_pmlq;
  occa::memory o_rhspmlq;
  occa::memory o_rkpmlq;
  occa::memory o_rkrhspmlq;

  occa::memory o_savepmlq;

  occa::kernel rkPmlUpdateKernel;
  occa::kernel rkPmlStageKernel;

  void Backup(occa::memory &o_Q);
  void Restore(occa::memory &o_Q);
  void AcceptStep(occa::memory &o_q, occa::memory &o_rq);

  void Step(occa::memory& o_q, dfloat time, dfloat dt);

public:
  sark5_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields,
            dfloat *_lambda, solver_t& solver);
  ~sark5_pml();
};


/* Multi-rate Adams-Bashforth, order 3 */
class mrab3_pml: public mrab3 {
private:
  dlong Npml;
  int Npmlfields;

  occa::memory o_pmlq;
  occa::memory o_rhspmlq0, o_rhspmlq;

  occa::kernel pmlUpdateKernel;

  void Step(occa::memory& o_q, dfloat time, dfloat dt, int order);

public:
  mrab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields, solver_t& solver);
  ~mrab3_pml();
};

/* Multi-rate Semi-Analytic Adams-Bashforth, order 3 */
// Note: PML fields are not stepped Semi-Analytically
class mrsaab3_pml: public mrsaab3 {
private:
  dlong Npml;
  int Npmlfields;

  occa::memory o_pmlq;
  dfloat *pmlsaab_a, *pmlsaab_b;
  occa::memory o_pmlsaab_a, o_pmlsaab_b;

  occa::memory o_rhspmlq0, o_rhspmlq;

  occa::kernel pmlUpdateKernel;

  void Step(occa::memory& o_q, dfloat time, dfloat dt, int order);

public:
  mrsaab3_pml(dlong Nelements, dlong NpmlElements, dlong NhaloElements,
            int Np, int Nfields, int _Npmlfields,
            dfloat *_lambda, solver_t& solver);
  ~mrsaab3_pml();
};

} //namespace TimeStepper

#endif
