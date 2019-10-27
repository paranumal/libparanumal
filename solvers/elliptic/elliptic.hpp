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

#ifndef ELLIPTIC_HPP
#define ELLIPTIC_HPP 1

#include "core.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "linAlg.hpp"
#include "precon.hpp"
#include "linearSolver.hpp"
#include "parAlmond.hpp"

#define DELLIPTIC LIBP_DIR"/solvers/elliptic/"

class ellipticSettings_t: public settings_t {
public:
  ellipticSettings_t(MPI_Comm& _comm);
  void report();
};

class elliptic_t: public solver_t {
public:
  int Nfields;

  dfloat lambda;
  dfloat tau;

  int disc_ipdg, disc_c0;

  halo_t* traceHalo;

  precon_t* precon;

  dfloat *grad;
  occa::memory o_grad;

  dfloat *weight, *weightG;
  occa::memory o_weight, o_weightG;

  //C0-FEM mask data
  ogs_t *ogsMasked;
  int *mapB;      // boundary flag of face nodes

  dlong Nmasked;
  dlong *maskIds;
  hlong *maskedGlobalIds;
  hlong *maskedGlobalNumbering;
  int   *maskedGlobalOwners;

  occa::memory o_maskIds;
  occa::memory o_mapB;

  int *BCType;
  int *EToB;
  occa::memory o_EToB;

  int allNeumann;
  dfloat allNeumannPenalty;
  dfloat allNeumannScale;

  occa::kernel maskKernel;
  occa::kernel partialAxKernel;
  occa::kernel partialGradientKernel;
  occa::kernel partialIpdgKernel;

  elliptic_t() = delete;
  elliptic_t(mesh_t& _mesh, linAlg_t& _linAlg, dfloat _lambda):
    solver_t(_mesh, _linAlg), lambda(_lambda) {}

  //setup
  static elliptic_t& Setup(mesh_t& mesh, linAlg_t& linAlg, dfloat lambda);

  void BoundarySetup();

  void Run();

  int Solve(linearSolver_t& linearSolver, occa::memory &o_x, occa::memory &o_r,
            const dfloat tol, const int MAXIT, const int verbose);

  void PlotFields(dfloat* Q, char *fileName);

  void Operator(occa::memory& o_q, occa::memory& o_Aq);

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

  void BuildOperatorDiagonal(dfloat *diagA);

  void BuildOperatorDiagonalContinuousTri2D(dfloat *diagA);
  void BuildOperatorDiagonalContinuousTri3D(dfloat *diagA);
  void BuildOperatorDiagonalContinuousQuad2D(dfloat *diagA);
  void BuildOperatorDiagonalContinuousQuad3D(dfloat *diagA);
  void BuildOperatorDiagonalContinuousTet3D(dfloat *diagA);
  void BuildOperatorDiagonalContinuousHex3D(dfloat *diagA);

  void BuildOperatorDiagonalIpdgTri2D(dfloat *diagA);
  void BuildOperatorDiagonalIpdgTri3D(dfloat *diagA);
  void BuildOperatorDiagonalIpdgQuad2D(dfloat *diagA);
  void BuildOperatorDiagonalIpdgQuad3D(dfloat *diagA);
  void BuildOperatorDiagonalIpdgTet3D(dfloat *diagA);
  void BuildOperatorDiagonalIpdgHex3D(dfloat *diagA);

  elliptic_t& SetupNewDegree(mesh_t& meshF);

  elliptic_t& SetupRingPatch(mesh_t& meshPatch);

  void ZeroMean(occa::memory &o_q);
};


#endif

