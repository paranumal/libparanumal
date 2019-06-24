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
  MassMatrixPrecon(elliptic_t& elliptic);
  void Operator(occa::memory& o_r, occa::memory& o_Mr);
};

//ParAlmond AMG preconditioner
class ParAlmondPrecon: public precon_t {
private:
  elliptic_t& elliptic;
  settings_t& settings;

  parAlmond::solver_t *parAlmondHandle;

  dfloat *xG, *rhsG;
  occa::memory o_xG, o_rhsG;

public:
  ParAlmondPrecon(elliptic_t& elliptic);
  void Operator(occa::memory& o_r, occa::memory& o_Mr);
};

// Matrix-free p-Multigrid levels followed by AMG
class MultiGridPrecon: public precon_t {
private:
  elliptic_t& elliptic;
  mesh_t& mesh;
  settings_t& settings;

  parAlmond::solver_t *parAlmondHandle;

public:
  MultiGridPrecon(elliptic_t& elliptic);
  void Operator(occa::memory& o_r, occa::memory& o_Mr);
};

// Cast problem into spectrally-equivalent N=1 FEM space and precondition with AMG
class SEMFEMPrecon: public precon_t {
private:
  elliptic_t& elliptic;
  mesh_t& mesh;
  settings_t& settings;

  parAlmond::solver_t *parAlmondHandle;

  occa::memory o_xG, o_rhsG;

  occa::memory o_zFEM, o_rFEM;
  occa::memory o_GzFEM, o_GrFEM;

  ogs_t *FEMogs;

  occa::kernel SEMFEMInterpKernel;
  occa::kernel SEMFEMAnterpKernel;

public:
  SEMFEMPrecon(elliptic_t& elliptic);
  void Operator(occa::memory& o_r, occa::memory& o_Mr);
};

#endif