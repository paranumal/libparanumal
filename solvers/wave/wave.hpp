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

#ifndef WAVE_HPP
#define WAVE_HPP

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"
#include "elliptic.hpp"
#include "initialGuess.hpp"


#define DWAVE LIBP_DIR"/solvers/wave/"

using namespace libp;

class waveSettings_t: public settings_t {
public:
  waveSettings_t(comm_t& _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename);

   ellipticSettings_t extractEllipticSettings();
};

class wave_t: public solver_t {
public:
   
   mesh_t mesh;
   
   ellipticSettings_t ellipticSettings;
   elliptic_t elliptic;
   
   linearSolver_t<dfloat>  linearSolver;

   int Nfields=1;
   int disc_c0=0;
   
   int Niter;

   dlong Nall;
   dlong NglobalDofs;
   
   dfloat ellipticTOL;
   dfloat tau;

   int Nstages;
   int embedded;
   dfloat gamma;
   dfloat lambdaSolve;
   
   linAlgMatrix_t<dfloat> alpha, beta, betahat;
   
   deviceMemory<dfloat> o_alphatilde;
   deviceMemory<dfloat> o_gammatilde;
   deviceMemory<dfloat> o_betatilde;
   deviceMemory<dfloat> o_betahattilde;
   deviceMemory<dfloat> o_alpha;
   deviceMemory<dfloat> o_beta;
   deviceMemory<dfloat> o_betaAlpha;
   deviceMemory<dfloat> o_betahatAlpha;
   deviceMemory<dfloat> o_betahat;
   deviceMemory<dfloat> o_gamma;
   
   memory<dfloat> DL;
   memory<dfloat> PL;
   memory<dfloat> DrhsL;
   memory<dfloat> DhatL;
   memory<dfloat> PhatL;

   memory<dfloat> WJ;
   memory<dfloat> invWJ;

   deviceMemory<dfloat> o_DL;
   deviceMemory<dfloat> o_PL;
   deviceMemory<dfloat> o_DrhsL;
   deviceMemory<dfloat> o_DhatL;
   deviceMemory<dfloat> o_PhatL;

   deviceMemory<dfloat> o_DtildeL;
   deviceMemory<dfloat> o_Dtilde;

   deviceMemory<dfloat> o_Drhs;
   deviceMemory<dfloat> o_scratch1;
   deviceMemory<dfloat> o_scratch2;
   deviceMemory<dfloat> o_scratch1L;
   deviceMemory<dfloat> o_scratch2L;
   
   deviceMemory<dfloat> o_invMM;
   deviceMemory<dfloat> o_MM;
   deviceMemory<dfloat> o_invWJ;
   deviceMemory<dfloat> o_WJ;

   kernel_t waveStageUpdateKernel;
   kernel_t waveCombineKernel;
   kernel_t waveErrorEstimateKernel;
   kernel_t waveStepInitializeKernel;
   kernel_t waveStepFinalizeKernel;
   kernel_t waveStageInitializeKernel;
   kernel_t waveStageFinalizeKernel;
   kernel_t waveInitialConditionsKernel;
   
   wave_t() = default;
   wave_t(platform_t &_platform,
          mesh_t &_mesh,
          waveSettings_t& _settings) {
     Setup(_platform, _mesh, _settings);
   }
   
   //setup
   void Setup(platform_t& _platform, mesh_t& _mesh,
              waveSettings_t& _settings);
   
   void Solve(deviceMemory<dfloat> &_DL,
              deviceMemory<dfloat> &_PL,
              const dfloat T);
   
   void Run();
   
   void Report(dfloat time, int tstep);
   
   void PlotFields(libp::memory<dfloat>& DL,
                   libp::memory<dfloat>& PL,
                   std::string fileName);

};

#endif
