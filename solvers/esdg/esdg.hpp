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

#ifndef ESDG_HPP
#define ESDG_HPP 1

#define mymin(a,b)  (((a)<(b))  ? (a) : (b))
#define mymax(a,b)  (((a)>(b))  ? (a) : (b))



#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"

#define DESDG LIBP_DIR"/solvers/esdg/"

using namespace libp;

class esdgSettings_t: public settings_t {
public:
  esdgSettings_t(comm_t _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename);
};

class esdg_t: public solver_t {
public:

  mesh_t mesh;
  
  int Nfields;
  int Ngrads;

  dfloat gamma;

  int cubature;
  int entropyStable;

  timeStepper_t timeStepper;

  ogs::halo_t fieldTraceHalo;
  ogs::halo_t gradTraceHalo;

  memory<dfloat> q, qexact;

  deviceMemory<dfloat> o_q;
  deviceMemory<dfloat> o_qexact;

  memory<dfloat> Vort;
  deviceMemory<dfloat> o_Vort;

  memory<dfloat> l2Error, linfError;
  deviceMemory<dfloat> o_l2Error, o_linfError;

  deviceMemory<dfloat> o_cubFilterT, o_cubDmatricesT, o_cubw; // cub2cub differentiation matrices, cubw missing ?
  
  deviceMemory<dfloat> o_Mq;

  // filter recon degree
  int reconN, reconNp;

  // cell centers (needed for Riemann probem initialization)
  deviceMemory<dfloat> o_cx, o_cy, o_cz;
  
  deviceMemory<dfloat> o_filterT;    // top mode filter on regular volume nodes

  kernel_t vorticityKernel;
  kernel_t dgVorticityKernel;

  kernel_t initialConditionKernel;

  kernel_t MassMatrixKernel;
  kernel_t filterKernel;
  kernel_t errorKernel;

  // =======>
  // entropy stable

  int esNp;

  deviceMemory<dfloat> o_esMu;   // variable viscosity (accomodates viscosity ramp at outflow)
  
  deviceMemory<dfloat> o_esSurfRHS;
  
  deviceMemory<dfloat> o_esX; // x coordinates of vol and surf es nodes
  deviceMemory<dfloat> o_esY; // y coordinates of vol and surf es nodes
  deviceMemory<dfloat> o_esZ; // z coordinates of vol and surf es nodes

  deviceMemory<dfloat> o_esRelaxOpT; // relaxation filter terms for entropy vars
  deviceMemory<dfloat> o_esIqT; // interpolate to volume quadratures
  deviceMemory<dfloat> o_esIqfT; // interpolate to combined quadratures
  deviceMemory<dfloat> o_esQNrT;  // quadrature to quadrature 'r' derivative
  deviceMemory<dfloat> o_esQNsT;  // quadrature to quadrature 's' derivative
  deviceMemory<dfloat> o_esPqT;   // quadrature to PN nodes
  deviceMemory<dfloat> o_esFqT;   // quadrature to quadrature through PN
  deviceMemory<dfloat> o_esPqLfT; // project from volume quad and lift from surface quadrature to PN nodes
  deviceMemory<dfloat> o_esItMT; // interpolate and integrate
  deviceMemory<dfloat> o_esLfT;   // lift from surface to PN nodes

  deviceMemory<dfloat> o_esSrrT;
  deviceMemory<dfloat> o_esSrsT;
  deviceMemory<dfloat> o_esSsrT;
  deviceMemory<dfloat> o_esSssT;
  deviceMemory<dfloat> o_esLumpedMassMatrixInv;
  deviceMemory<dfloat> o_esLocalMassMatrixInv;
  
  deviceMemory<dfloat> o_esIqfDrPqT; // diffferentiate from vol quadrature to all quadrature via PN
  deviceMemory<dfloat> o_esIqfDsPqT; // diffferentiate from vol quadrature to all quadrature via PN
  deviceMemory<dfloat> o_esIqfLfT; // interpolate from surface quadrature to all quadature via PN
  
  deviceMemory<dfloat> o_esQc;      // storage for conserved variables at combined quadrature
  deviceMemory<dfloat> o_esQe;      // storage for entropy variables at combined quadrature
  deviceMemory<dfloat> o_esQp;      // storage for primitive variables at combined quadrature
  deviceMemory<dfloat> o_esQcrecon; // storage for reconstructed conservative variables at combined quadrature

  deviceMemory<dfloat> o_esdQedx;  // 'x' derivative of entropy variables at combined quadrature
  deviceMemory<dfloat> o_esdQedy;  // 'y' derivative of entropy variables at combined quadrature
  deviceMemory<dfloat> o_esdQedz;  // 'z' derivative of entropy variables at combined quadrature

  deviceMemory<dfloat> o_esDQe;  // art diff P1 subcell
  
  memory<dfloat> entropyChange;
  deviceMemory<dfloat> o_entropyChange;

  
  deviceMemory<dlong> o_esVmapM;
  deviceMemory<dlong> o_esVmapP;

  deviceMemory<dfloat> o_esR;
  deviceMemory<dfloat> o_esS;
  deviceMemory<dfloat> o_esRq;
  deviceMemory<dfloat> o_esSq;
  deviceMemory<dfloat> o_esWq;
  deviceMemory<dfloat> o_esWf;

  deviceMemory<dlong> o_EToE;
  deviceMemory<dlong> o_EToB;
  
  deviceMemory<dfloat> o_esTotalEntropy;
  memory<dfloat> esTotalEntropy;
  
  kernel_t esInterpolateKernel;
  kernel_t esVolumeCubatureKernel;
  kernel_t esSurfaceCubatureKernel;

  kernel_t esVolumeKernel;
  kernel_t esSurfaceKernel;

  kernel_t esIntegrateEntropyChangeKernel;
  kernel_t esIntegrateEntropyKernel;

  kernel_t esVolumeGradientKernel;
  kernel_t esSurfaceGradientKernel;

  kernel_t esDiffusionFluxesKernel;
  
  kernel_t esVolumeDivergenceKernel;
  kernel_t esSurfaceDivergenceKernel;

  
  //  kernel_t esProjectKernel;
  // <=======
  



  esdg_t() = delete;
  esdg_t(platform_t &_platform, mesh_t& _mesh, esdgSettings_t& _settings){
    Setup(_platform, _mesh, _settings);
  }

  ~esdg_t();

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh, 
	     esdgSettings_t& _settings);

  void SetupTri2D();
  void SetupQuad2D();

  
  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(memory<dfloat> &Q, memory<dfloat> &V, char *fileName);

  void rhsf(deviceMemory<dfloat> &o_q, deviceMemory<dfloat> &o_rhs, const dfloat time);

  void saveCheckpoint(memory<dfloat> &outq, dfloat time);

  void loadCheckpoint(memory<dfloat> &inq, dfloat &time);

  dfloat integrateEntropy(deviceMemory<dfloat> &o_Q);
  dfloat integrateEntropyChange(deviceMemory<dfloat> &o_Q, deviceMemory<dfloat>& o_RHS);
};

#endif
