/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Anthony Austin

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

namespace libp {

template <typename T>
int linearSolver_t<T>::Solve(operator_t& linearOperator,
			     operator_t& precon,
			     deviceMemory<T>& o_x,
			     deviceMemory<T>& o_rhs,
			     const T tol,
			     const int MAXIT,
			     const int verbose) {
  assertInitialized();
  ig->FormInitialGuess(o_x, o_rhs);
  int iters = ls->Solve(linearOperator, precon, o_x, o_rhs, tol, MAXIT, verbose);
  ig->Update(linearOperator, o_x, o_rhs);

  return iters;
}

template <typename T>
void linearSolver_t<T>::MakeDefaultInitialGuessStrategy() {
  ig = std::make_shared<InitialGuess::Zero<T> >(ls->N, ls->platform, ls->settings, ls->comm);
}

template <typename T>
bool linearSolver_t<T>::isInitialized() {
  return (ls!=nullptr && ig!=nullptr);
}

template <typename T>
void linearSolver_t<T>::assertInitialized() {
  LIBP_ABORT("LinearSolver not initialized",
             ls==nullptr);
  LIBP_ABORT("InitialGuess not initialized",
             ig==nullptr);
}

template class linearSolver_t<double>;
template class linearSolver_t<float>;

} //namespace libp
