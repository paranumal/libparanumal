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

#ifndef ELLIPTICPRECON_HPP
#define ELLIPTICPRECON_HPP 1

#include "elliptic.hpp"
#include "parAlmond.hpp"

//Jacobi preconditioner
class JacobiPrecon: public operator_t {
private:
	elliptic_t elliptic;

  deviceMemory<dfloat> o_invDiagA;

public:
  JacobiPrecon() = default;
  JacobiPrecon(elliptic_t& elliptic);
  void Operator(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_Mr);
};

//Inverse Mass Matrix preconditioner
class MassMatrixPrecon: public operator_t {
private:
  elliptic_t elliptic;
  mesh_t mesh;
  settings_t settings;

  deviceMemory<dfloat> o_MrL, o_rtmp;
  deviceMemory<dfloat> o_invMM;

  kernel_t blockJacobiKernel;
  kernel_t partialBlockJacobiKernel;

public:
  MassMatrixPrecon() = default;
  MassMatrixPrecon(elliptic_t& elliptic);
  void Operator(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_Mr);
};

//ParAlmond AMG preconditioner
class ParAlmondPrecon: public operator_t {
private:
  elliptic_t elliptic;
  settings_t settings;

  parAlmond::parAlmond_t parAlmond;

  memory<dfloat> xG, rhsG;
  deviceMemory<dfloat> o_xG, o_rhsG;

public:
  ParAlmondPrecon() = default;
  ParAlmondPrecon(elliptic_t& elliptic);
  void Operator(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_Mr);
};

// Matrix-free p-Multigrid levels followed by AMG
class MultiGridPrecon: public operator_t {
private:
  elliptic_t elliptic;
  mesh_t mesh;
  settings_t settings;

  parAlmond::parAlmond_t parAlmond;

public:
  MultiGridPrecon() = default;
  MultiGridPrecon(elliptic_t& elliptic);
  void Operator(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_Mr);
};

// Cast problem into spectrally-equivalent N=1 FEM space and precondition with AMG
class SEMFEMPrecon: public operator_t {
private:
  elliptic_t elliptic;
  mesh_t mesh;
  settings_t settings;

  mesh_t femMesh;
  elliptic_t femElliptic;
  parAlmond::parAlmond_t parAlmond;

  deviceMemory<dfloat> o_MrL;

  deviceMemory<dfloat> o_zFEM, o_rFEM;
  deviceMemory<dfloat> o_GzFEM, o_GrFEM;

  ogs::ogs_t FEMogs;

  kernel_t SEMFEMInterpKernel;
  kernel_t SEMFEMAnterpKernel;

public:
  SEMFEMPrecon() = default;
  SEMFEMPrecon(elliptic_t& elliptic);
  void Operator(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_Mr);
};


class MGLevel: public parAlmond::multigridLevel {
public:
  elliptic_t elliptic;
  mesh_t mesh;

  //prologation
  memory<dfloat> P;
  deviceMemory<dfloat> o_P;

  kernel_t coarsenKernel, partialCoarsenKernel;
  kernel_t prolongateKernel, partialProlongateKernel;

  //coarse space
  elliptic_t ellipticC;
  mesh_t meshC;

  //smoothing params
  typedef enum {JACOBI=1,
                CHEBYSHEV=2} SmootherType;
  SmootherType stype;

  dfloat lambda1, lambda0;
  int ChebyshevIterations;

  static dlong NsmootherResidual, Nscratch;
  static memory<dfloat> smootherResidual;
  static deviceMemory<dfloat> o_smootherResidual;
  static deviceMemory<dfloat> o_smootherResidual2;
  static deviceMemory<dfloat> o_smootherUpdate;
  static deviceMemory<dfloat> o_transferScratch;

  //jacobi data
  deviceMemory<dfloat> o_invDiagA;

  //build a p-multigrid level and connect it to the next one
  MGLevel() = default;
  MGLevel(elliptic_t& _elliptic,
          dlong _Nrows, dlong _Ncols,
          int Nc, int NpCoarse);

  void Operator(deviceMemory<dfloat> &o_X, deviceMemory<dfloat> &o_Ax);

  void residual(deviceMemory<dfloat> &o_RHS, deviceMemory<dfloat> &o_X, deviceMemory<dfloat> &o_RES);

  void coarsen(deviceMemory<dfloat> &o_X, deviceMemory<dfloat> &o_Cx);

  void prolongate(deviceMemory<dfloat> &o_X, deviceMemory<dfloat> &o_Px);

  //smoother ops
  void smooth(deviceMemory<dfloat> &o_RHS, deviceMemory<dfloat> &o_X, bool x_is_zero);

  void smoothJacobi    (deviceMemory<dfloat> &o_r, deviceMemory<dfloat> &o_X, bool xIsZero);
  void smoothChebyshev (deviceMemory<dfloat> &o_r, deviceMemory<dfloat> &o_X, bool xIsZero);

  void Report();

  void SetupSmoother();
  dfloat maxEigSmoothAx();

  void AllocateStorage();
};

// Overlapping additive Schwarz with patch problems consisting of the
//  entire local mesh + 1 ring overlap, solved with a local multigrid
//  precon and coarse problem consisting of the global degree 1
//  problem, solved with parAlmond
class OASPrecon: public operator_t {
private:
  elliptic_t elliptic;
  mesh_t mesh;
  settings_t settings;

  //Patch precon
  mesh_t meshPatch;
  elliptic_t ellipticPatch;
  precon_t preconPatch;
  MGLevel level;

  ogs::ogs_t ogsMaskedRing; //ogs for 1-ring patch

  //Coarse Precon
  ogs::ogs_t ogsMasked;
  parAlmond::parAlmond_t parAlmond;

  memory<dfloat> rPatch, zPatch;
  memory<dfloat> rPatchL, zPatchL;
  deviceMemory<dfloat> o_rPatch, o_zPatch;
  deviceMemory<dfloat> o_rPatchL, o_zPatchL;

  memory<dfloat> rC, zC;
  deviceMemory<dfloat> o_rC, o_zC;

  memory<dfloat> patchWeight;
  deviceMemory<dfloat> o_patchWeight;

public:
  OASPrecon() = default;
  OASPrecon(elliptic_t& elliptic);
  void Operator(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_Mr);
};


#endif
