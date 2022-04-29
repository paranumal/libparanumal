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

#ifndef ELLIPTIC_HPP
#define ELLIPTIC_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "linAlg.hpp"
#include "precon.hpp"
#include "linearSolver.hpp"
#include "parAlmond.hpp"

#define DELLIPTIC LIBP_DIR"/solvers/elliptic/"

using namespace libp;

class ellipticSettings_t: public settings_t {
public:
  ellipticSettings_t() = default;
  ellipticSettings_t(const comm_t& _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename);
};
void ellipticAddRunSettings(settings_t& settings);
void ellipticAddSettings(settings_t& settings,
                         const std::string prefix="");

class elliptic_t: public solver_t {
public:
  mesh_t mesh;

  dlong Ndofs, Nhalo;
  int Nfields;

  dfloat lambda;
  dfloat tau;

  int disc_ipdg, disc_c0;

  deviceMemory<dfloat> o_AqL;

  ogs::halo_t traceHalo;

  precon_t precon;

  memory<dfloat> grad;
  deviceMemory<dfloat> o_grad;

  memory<dfloat> weight, weightG;
  deviceMemory<dfloat> o_weight, o_weightG;

  //C0-FEM mask data
  ogs::ogs_t ogsMasked;
  ogs::halo_t gHalo;
  memory<int> mapB;      // boundary flag of face nodes
  deviceMemory<int> o_mapB;

  dlong Nmasked;
  memory<dlong> maskIds;
  memory<hlong> maskedGlobalIds;
  memory<hlong> maskedGlobalNumbering;
  memory<dlong> GlobalToLocal;

  deviceMemory<dlong> o_maskIds;
  deviceMemory<dlong> o_GlobalToLocal;

  int NBCTypes;
  memory<int> BCType;
  memory<int> EToB;
  deviceMemory<int> o_EToB;

  int allNeumann;
  dfloat allNeumannPenalty;
  dfloat allNeumannScale;

  kernel_t maskKernel;
  kernel_t partialAxKernel;
  kernel_t partialGradientKernel;
  kernel_t partialIpdgKernel;

  elliptic_t() = default;
  elliptic_t(platform_t &_platform, mesh_t &_mesh,
              settings_t& _settings, dfloat _lambda,
              const int _NBCTypes, const memory<int> _BCType) {
    Setup(_platform, _mesh, _settings, _lambda, _NBCTypes, _BCType);
  }

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh,
             settings_t& _settings, dfloat _lambda,
             const int _NBCTypes, const memory<int> _BCType);

  void BoundarySetup();

  void Run();

  int Solve(linearSolver_t& linearSolver, deviceMemory<dfloat> &o_x, deviceMemory<dfloat> &o_r,
            const dfloat tol, const int MAXIT, const int verbose);

  void PlotFields(memory<dfloat>& Q, std::string fileName);

  void Operator(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_Aq);

  void BuildOperatorMatrixIpdg(parAlmond::parCOO& A);
  void BuildOperatorMatrixContinuous(parAlmond::parCOO& A);

  void BuildOperatorMatrixContinuousTri2D(parAlmond::parCOO& A);
  void BuildOperatorMatrixContinuousTri3D(parAlmond::parCOO& A);
  void BuildOperatorMatrixContinuousQuad2D(parAlmond::parCOO& A);
  void BuildOperatorMatrixContinuousQuad3D(parAlmond::parCOO& A);
  void BuildOperatorMatrixContinuousTet3D(parAlmond::parCOO& A);
  void BuildOperatorMatrixContinuousHex3D(parAlmond::parCOO& A);

  void BuildOperatorMatrixIpdgTri2D(parAlmond::parCOO& A);
  void BuildOperatorMatrixIpdgTri3D(parAlmond::parCOO& A);
  void BuildOperatorMatrixIpdgQuad2D(parAlmond::parCOO& A);
  void BuildOperatorMatrixIpdgQuad3D(parAlmond::parCOO& A);
  void BuildOperatorMatrixIpdgTet3D(parAlmond::parCOO& A);
  void BuildOperatorMatrixIpdgHex3D(parAlmond::parCOO& A);

  void BuildOperatorDiagonal(memory<dfloat>& diagA);

  void BuildOperatorDiagonalContinuousTri2D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalContinuousTri3D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalContinuousQuad2D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalContinuousQuad3D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalContinuousTet3D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalContinuousHex3D(memory<dfloat>& diagA);

  void BuildOperatorDiagonalIpdgTri2D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalIpdgTri3D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalIpdgQuad2D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalIpdgQuad3D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalIpdgTet3D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalIpdgHex3D(memory<dfloat>& diagA);

  elliptic_t SetupNewDegree(mesh_t& meshF);

  elliptic_t SetupRingPatch(mesh_t& meshPatch);

  void ZeroMean(deviceMemory<dfloat> &o_q);
};


#endif

