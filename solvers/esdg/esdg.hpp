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

#ifndef ESDG_HPP
#define ESDG_HPP 1

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

  dfloat mu;
  dfloat gamma;

  int cubature;
  int isothermal;

  timeStepper_t timeStepper;

  ogs::halo_t fieldTraceHalo;
  ogs::halo_t gradTraceHalo;

  memory<dfloat> q;
  deviceMemory<dfloat> o_q;

  kernel_t volumeKernel;
  kernel_t surfaceKernel;
  kernel_t cubatureVolumeKernel;
  kernel_t cubatureSurfaceKernel;

  kernel_t gradVolumeKernel;
  kernel_t gradSurfaceKernel;

  kernel_t vorticityKernel;

  kernel_t constrainKernel;

  kernel_t initialConditionKernel;
  kernel_t maxWaveSpeedKernel;



  esdg_t() = default;
  esdg_t(platform_t &_platform, mesh_t &_mesh,
              esdgSettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh,
             esdgSettings_t& _settings);

  void Run();

  void Report(dfloat time, int tstep) override;

  void PlotFields(memory<dfloat> Q, memory<dfloat> V, std::string fileName);

  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);

  dfloat MaxWaveSpeed(deviceMemory<dfloat>& o_Q, const dfloat T);

  // ESDG specifics
  int esNp;
  dlong maxLVToE;
  dlong maxEToVToE;
  deviceMemory<dlong> o_LVToE;
  deviceMemory<dlong> o_NLVToE;

  deviceMemory<dlong> o_NEToVToE;
  deviceMemory<dlong> o_EToVToE;
  
  deviceMemory<dlong> o_EToE;
  deviceMemory<dlong> o_EToB;
  
  deviceMemory<dlong>  o_esVmapM;
  deviceMemory<dlong>  o_esVmapP; 
  
  deviceMemory<dfloat> o_esIqfLfT;
  deviceMemory<dfloat> o_esIqfDrPqT;
  deviceMemory<dfloat> o_esIqfDsPqT;
  
  deviceMemory<dfloat> o_esRelaxOpT;
  deviceMemory<dfloat> o_esMu;
  deviceMemory<dfloat> o_esSurfRHS;
  deviceMemory<dfloat> o_esX;
  deviceMemory<dfloat> o_esY;
  deviceMemory<dfloat> o_esZ;
  
  deviceMemory<dfloat> o_esIqT;
  deviceMemory<dfloat> o_esIqfT;
  deviceMemory<dfloat> o_esQNrT;
  deviceMemory<dfloat> o_esQNsT;
  deviceMemory<dfloat> o_esPqT;
  deviceMemory<dfloat> o_esPqLfT;
  deviceMemory<dfloat> o_esFqT;
  deviceMemory<dfloat> o_esLfT;
  deviceMemory<dfloat> o_esItMT;
  deviceMemory<dfloat> o_esQc;
  deviceMemory<dfloat> o_esQe;
  deviceMemory<dfloat> o_esQp;
  
  deviceMemory<dfloat> o_esQcrecon;
  
  deviceMemory<dfloat> o_esR;
  deviceMemory<dfloat> o_esS;
  
  deviceMemory<dfloat> o_esRq;
  deviceMemory<dfloat> o_esSq;
  deviceMemory<dfloat> o_esWq;
  deviceMemory<dfloat> o_esWf;
  
  deviceMemory<dfloat> o_esdQedx;
  deviceMemory<dfloat> o_esdQedy;
  deviceMemory<dfloat> o_esdQedz;
  
  deviceMemory<dfloat> o_esDQe;

  deviceMemory<dfloat> o_cx;
  deviceMemory<dfloat> o_cy;
  deviceMemory<dfloat> o_cz;
  
  kernel_t esInterpolateKernel;
  kernel_t esVolumeCubatureKernel;
  kernel_t esSurfaceCubatureKernel;

  kernel_t esVolumeKernel;
  kernel_t esSurfaceKernel;

  kernel_t esIntegrateEntropyChangeKernel;
  kernel_t esIntegrateEntropyKernel;

  kernel_t esRelaxationKernel;
  
  kernel_t esVolumeGradientKernel;
  kernel_t esSurfaceGradientKernel;

  kernel_t esDiffusionFluxesKernel;
  
  kernel_t esVolumeDivergenceKernel;
  kernel_t esSurfaceDivergenceKernel;

  
};

#endif
