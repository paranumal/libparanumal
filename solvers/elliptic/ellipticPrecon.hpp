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

#ifndef ELLIPTICPRECON_HPP
#define ELLIPTICPRECON_HPP 1

#include "elliptic.hpp"
#include "parAlmond.hpp"

//Jacobi preconditioner
class JacobiPrecon: public precon_t {
private:
	elliptic_t& elliptic;

  occa::memory o_invDiagA;

public:
  JacobiPrecon(elliptic_t& elliptic);
  void Operator(occa::memory& o_r, occa::memory& o_Mr);
};

//Inverse Mass Matrix preconditioner
class MassMatrixPrecon: public precon_t {
private:
  elliptic_t& elliptic;
  mesh_t& mesh;
  settings_t& settings;

  occa::memory o_rtmp;
  occa::memory o_invMM;

  occa::kernel blockJacobiKernel;
  occa::kernel partialBlockJacobiKernel;

public:
  ~MassMatrixPrecon();
  MassMatrixPrecon(elliptic_t& elliptic);
  void Operator(occa::memory& o_r, occa::memory& o_Mr);
};

//ParAlmond AMG preconditioner
class ParAlmondPrecon: public precon_t {
private:
  elliptic_t& elliptic;
  settings_t& settings;

  parAlmond::parAlmond_t parAlmond;

  dfloat *xG, *rhsG;
  occa::memory o_xG, o_rhsG;

public:
  ~ParAlmondPrecon();
  ParAlmondPrecon(elliptic_t& elliptic);
  void Operator(occa::memory& o_r, occa::memory& o_Mr);
};

// Matrix-free p-Multigrid levels followed by AMG
class MultiGridPrecon: public precon_t {
private:
  elliptic_t& elliptic;
  mesh_t& mesh;
  settings_t& settings;

  parAlmond::parAlmond_t parAlmond;

  bool gather=false;
  occa::memory o_xG, o_rhsG;

public:
  MultiGridPrecon(elliptic_t& elliptic);
  ~MultiGridPrecon() = default;
  void Operator(occa::memory& o_r, occa::memory& o_Mr);
};

// Cast problem into spectrally-equivalent N=1 FEM space and precondition with AMG
class SEMFEMPrecon: public precon_t {
private:
  elliptic_t& elliptic;
  mesh_t& mesh;
  settings_t& settings;

  mesh_t *femMesh;
  elliptic_t* femElliptic;
  parAlmond::parAlmond_t parAlmond;

  occa::memory o_xG, o_rhsG;

  occa::memory o_zFEM, o_rFEM;
  occa::memory o_GzFEM, o_GrFEM;

  ogs_t *FEMogs;

  occa::kernel SEMFEMInterpKernel;
  occa::kernel SEMFEMAnterpKernel;

public:
  ~SEMFEMPrecon();
  SEMFEMPrecon(elliptic_t& elliptic);
  void Operator(occa::memory& o_r, occa::memory& o_Mr);
};

class MGLevel;
// Overlapping additive Schwarz with patch problems consisting of the
//  entire local mesh + 1 ring overlap, solved with a local multigrid
//  precon and coarse problem consisting of the global degree 1
//  problem, solved with parAlmond
class OASPrecon: public precon_t {
private:
  elliptic_t& elliptic;
  mesh_t& mesh;
  settings_t& settings;

  //Patch precon
  mesh_t* meshPatch;
  elliptic_t* ellipticPatch;
  precon_t *preconPatch;
  MGLevel *level;

  ogs_t *ogsMaskedRing; //ogs for 1-ring patch

  //Coarse Precon
  ogs_t *ogsMasked=nullptr;
  parAlmond::parAlmond_t parAlmond;
  occa::memory o_xG, o_rhsG;

  dfloat *rPatch, *zPatch;
  occa::memory o_rPatch, o_zPatch;

  dfloat *rC, *zC;
  occa::memory o_rC, o_zC;

  dfloat *patchWeight;
  occa::memory o_patchWeight;

public:
  ~OASPrecon();
  OASPrecon(elliptic_t& elliptic);
  void Operator(occa::memory& o_r, occa::memory& o_Mr);
};

class MGLevel: public parAlmond::multigridLevel {
public:
  elliptic_t& elliptic;
  mesh_t& mesh;
  linAlg_t& linAlg;

  //prologation
  dfloat *P;
  occa::memory o_P;

  occa::kernel coarsenKernel;
  occa::kernel prolongateKernel;

  bool gatherLevel=false;
  ogs_t *ogsMasked=nullptr;
  occa::memory o_SX, o_GX;

  //masking data
  dlong Nmasked;
  occa::memory o_maskIds;

  //smoothing params
  typedef enum {JACOBI=1,
                CHEBYSHEV=2} SmootherType;
  SmootherType stype;

  dfloat lambda1, lambda0;
  int ChebyshevIterations;

  static size_t smootherResidualBytes;
  static dfloat *smootherResidual;
  static occa::memory o_smootherResidual;
  static occa::memory o_smootherResidual2;
  static occa::memory o_smootherUpdate;

  //jacobi data
  occa::memory o_invDiagA;

  //build a p-multigrid level and connect it to the next one
  MGLevel(elliptic_t& _elliptic, int Nc, int NpCoarse);
  ~MGLevel();

  void Operator(occa::memory &o_X, occa::memory &o_Ax);

  void residual(occa::memory &o_RHS, occa::memory &o_X, occa::memory &o_RES);

  void coarsen(occa::memory &o_X, occa::memory &o_Cx);

  void prolongate(occa::memory &o_X, occa::memory &o_Px);

  //smoother ops
  void smooth(occa::memory &o_RHS, occa::memory &o_X, bool x_is_zero);

  void smoothJacobi    (occa::memory &o_r, occa::memory &o_X, bool xIsZero);
  void smoothChebyshev (occa::memory &o_r, occa::memory &o_X, bool xIsZero);

  void Report();

  void SetupSmoother();
  dfloat maxEigSmoothAx();

  void AllocateStorage();
};


#endif