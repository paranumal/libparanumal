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

#include "core.hpp"
#include "timeStepper.hpp"

namespace libp {

void timeStepper_t::Run(solver_t& solver,
                        deviceMemory<dfloat>& o_q,
                        dfloat start, dfloat end) {
  assertInitialized();
  ts->Run(solver, o_q, start, end);
}

void timeStepper_t::SetTimeStep(dfloat dt_) {
  assertInitialized();
  ts->SetTimeStep(dt_);
}

dfloat timeStepper_t::GetTimeStep() {
  assertInitialized();
  return ts->dt;
}

dfloat timeStepper_t::GetGamma() {
  assertInitialized();
  return ts->GetGamma();
}

void timeStepper_t::assertInitialized() {
  LIBP_ABORT("timeStepper_t not initialized",
             ts==nullptr);
}

} //namespace libp
