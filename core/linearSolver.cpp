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

#include "linearSolver.hpp"

//virtual base time stepper class
linearSolver_t::linearSolver_t(dlong _N, solver_t& _solver):
  N(_N),
  solver(_solver),
  comm(_solver.comm),
  device(_solver.device),
  settings(_solver.settings),
  props(_solver.props) {}

linearSolver_t* linearSolver_t::Setup(dlong N, solver_t& solver) {

  linearSolver_t *linearSolver=NULL;

  if (solver.settings.compareSetting("LINEAR SOLVER","PCG")){
    linearSolver = new pcg(N, solver);
  } else {
    LIBP_ABORT(string("Requested LINEAR SOLVER not found."));
  }

  return linearSolver;
}