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

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "mesh.hpp"
#include "linAlg.hpp"

class solver_t {
public:
  mesh_t& mesh;

  MPI_Comm& comm;
  occa::device& device;
  settings_t& settings;
  occa::properties& props;
  linAlg_t& linAlg;

  solver_t() = delete;

  solver_t(mesh_t& _mesh, linAlg_t& _linAlg):
    mesh(_mesh),
    comm(_mesh.comm),
    device(_mesh.device),
    settings(_mesh.settings),
    props(_mesh.props),
    linAlg(_linAlg) {};

  virtual ~solver_t(){}

  virtual void Run()=0;
  virtual void Report(dfloat time=0.0, int tstep=0) {};

  virtual void rhsf(occa::memory& o_q, occa::memory& o_rhs, const dfloat time) {};
  virtual void Operator(occa::memory& o_q, occa::memory& o_Aq) {};
};

#endif