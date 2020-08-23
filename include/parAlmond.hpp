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

#ifndef PARALMOND_HPP
#define PARALMOND_HPP

#include "core.hpp"
#include "ogs.hpp"
#include "settings.hpp"

#include "parAlmond/defines.hpp"
#include "parAlmond/utils.hpp"
#include "parAlmond/kernels.hpp"
#include "parAlmond/vector.hpp"
#include "parAlmond/matrix.hpp"
#include "parAlmond/level.hpp"
#include "parAlmond/agmg.hpp"
#include "parAlmond/coarse.hpp"
#include "parAlmond/solver.hpp"


namespace parAlmond {

solver_t *Init(platform_t& _platform, settings_t& _settings, MPI_Comm _comm);

void AMGSetup(solver_t* M,
             parCOO& A,
             bool nullSpace,
             dfloat nullSpacePenalty);

void Precon(solver_t* M, occa::memory o_x, occa::memory o_rhs);

void Report(solver_t *M);

void Free(solver_t* M);

} //namespace parAlmond

void parAlmondAddSettings(settings_t& settings, const string prefix="");
void parAlmondReportSettings(settings_t& settings);

#endif
