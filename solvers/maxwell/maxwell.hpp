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

#ifndef MAXWELL_HPP
#define MAXWELL_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"

#define DMAXWELL LIBP_DIR"/solvers/maxwell/"

using namespace libp;

class maxwellSettings_t: public settings_t {
public:
  maxwellSettings_t(comm_t _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename);
};

typedef enum{
  ISOTROPIC=1,
  HETEROGENEOUS=2
}materialType_e;

class maxwell_t: public solver_t {
public:
  mesh_t mesh;

  int Nfields;
  int Npmlfields;

  int materialType;
  
  timeStepper_t timeStepper;

  ogs::halo_t traceHalo;

  memory<dfloat> q;
  deviceMemory<dfloat> o_q;

  // Pml
  int pmlOrder;
  dfloat  sigmaXmax, sigmaYmax, sigmaZmax;
  memory<dfloat> pmlSigma;
  deviceMemory<dfloat> o_pmlSigma;
  dfloat pmlAlpha;

  memory<dfloat> pmlq;
  deviceMemory<dfloat> o_pmlq;
  
  // Flag for using cubature integration for sigma terms in pml
  int pmlcubature;
  
  int   materialNfields;
  dlong materialNlocal;
  dlong materialNhalo;
  dlong materialNtotal;
  memory<dfloat> materialCoefficients;
  memory<dfloat> materialInverseWeights;
  memory<dfloat> materialUpwindWeights;
  
  deviceMemory<dfloat> o_materialInverseWeights;
  deviceMemory<dfloat> o_materialUpwindWeights;
  
  ogs::halo_t materialHalo;
  
  kernel_t volumeKernel;
  kernel_t surfaceKernel;

  kernel_t pmlVolumeKernel;
  kernel_t pmlSurfaceKernel;
  kernel_t pmlCubatureTermsKernel;
  
  kernel_t heterogeneousSurfaceKernel;
  kernel_t heterogeneousProjectKernel;

  kernel_t initialConditionKernel;
  kernel_t errorKernel;

  maxwell_t() = default;
  maxwell_t(platform_t &_platform, mesh_t &_mesh,
              maxwellSettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh,
             maxwellSettings_t& _settings);

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(memory<dfloat> Q, const std::string fileName);

  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);
  void rhsf_pml(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_pmlQ,
                deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_pmlRHS, const dfloat T);

  void rhsVolume(dlong N, deviceMemory<dlong>& o_ids,
		 deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS);

  void rhsSurface(dlong N, deviceMemory<dlong>& o_ids,
		  deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T);


  void rhsPmlVolume(dlong N, deviceMemory<dlong>& o_ids, deviceMemory<dlong>& o_pmlids,
		    deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_pmlQ,
		    deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_pmlRHS);

  
  void rhsPmlSurface(dlong N, deviceMemory<dlong>& o_ids, deviceMemory<dlong>& o_pmlids,
		     deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_pmlQ,
		     deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_pmlRHS, const dfloat T);


  dfloat MaxWaveSpeed();

  void PmlSetup();
};

#endif
