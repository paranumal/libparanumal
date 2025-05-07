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

#ifndef MASS_HPP
#define MASS_HPP 1

#define DINS LIBP_DIR"/solvers/ins/"

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "linAlg.hpp"
#include "precon.hpp"
#include "linearSolver.hpp"
#include "parAlmond.hpp"

#define DMASS LIBP_DIR"/solvers/mass/"

using namespace libp;

class mass_t: public solver_t {
public:
  mesh_t mesh;

  dlong Ndofs, Nhalo;
  int Nfields;

  ogs::halo_t traceHalo;

  precon_t precon;

  // NOTE pfloat
  memory<pfloat> weight, weightG;
  deviceMemory<pfloat> o_weight, o_weightG;

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

  kernel_t reciprocalKernel;
  kernel_t maskKernel;
  kernel_t massAxKernel;
  kernel_t floatMassAxKernel;
  kernel_t massPartialAxKernel;
  kernel_t floatMassPartialAxKernel;
  kernel_t massScatterKernel;

  kernel_t buildOperatorDiagonalKernel;
  
  memory<dfloat> nut;
  deviceMemory<dfloat> o_nut;
  
  mass_t() = default;
  mass_t(platform_t &_platform, mesh_t &_mesh,
	 const int _NBCTypes, const memory<int> _BCType) {
    Setup(_platform, _mesh, _NBCTypes, _BCType);
  }

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh, 
             const int _NBCTypes, const memory<int> _BCType);

  void BoundarySetup();

  void Run();

  int Solve(linearSolver_t<dfloat>& linearSolver, deviceMemory<dfloat> &o_x, deviceMemory<dfloat> &o_r,
            const dfloat tol, const int MAXIT, const int verbose);

  void PlotFields(memory<dfloat>& Q, std::string fileName);

  void Operator(deviceMemory<double>& o_q, deviceMemory<double>& o_Aq);
  void Operator(deviceMemory<float>& o_q, deviceMemory<float>& o_Aq);

  void BuildOperatorDiagonal(memory<dfloat>& diagA);
  void BuildOperatorDiagonal(deviceMemory<pfloat> &o_invDiagA );
  
  void BuildOperatorDiagonalContinuousTri2D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalContinuousTri3D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalContinuousQuad2D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalContinuousQuad3D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalContinuousTet3D(memory<dfloat>& diagA);
  void BuildOperatorDiagonalContinuousHex3D(memory<dfloat>& diagA);

  void ZeroMean(deviceMemory<double> &o_q);
  void ZeroMean(deviceMemory<float> &o_q);
};

//MassJacobi preconditioner
class MassJacobiPrecon: public operator_t {
private:
	mass_t mass;

  deviceMemory<pfloat> o_invDiagA;

public:
  MassJacobiPrecon() = default;
  MassJacobiPrecon(mass_t& mass);
  void Operator(deviceMemory<pfloat>& o_r, deviceMemory<pfloat>& o_Mr);
  void Update();
};


#endif

