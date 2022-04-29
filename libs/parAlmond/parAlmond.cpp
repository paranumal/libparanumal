/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"
#include "parAlmond/parAlmondKernels.hpp"
#include "parAlmond/parAlmondCoarseSolver.hpp"

namespace libp {

namespace parAlmond {

void parAlmond_t::Setup(platform_t& _platform, settings_t& _settings, comm_t comm) {

  platform = _platform;
  settings = _settings;

  platform.linAlg().InitKernels({"set", "add", "sum", "scale",
                                "axpy", "zaxpy",
                                "amx", "amxpy", "zamxpy",
                                "adx", "adxpy", "zadxpy",
                                "innerProd", "norm2"});

  multigrid = std::make_shared<multigrid_t>(platform, settings, comm);

  //build parAlmond kernels on first construction
  buildParAlmondKernels(platform);
}

void parAlmond_t::Operator(deviceMemory<dfloat>& o_rhs, deviceMemory<dfloat>& o_x) {

  if (multigrid->exact){ //call the linear solver
    int maxIter = 500;
    int verbose = settings.compareSetting("VERBOSE", "TRUE") ? 1 : 0;
    dfloat tol = 1e-8;
    solver_t &A = multigrid->GetLevel<solver_t>(0);
    (void) multigrid->linearSolver.Solve(A, *multigrid, o_x, o_rhs, tol, maxIter, verbose);
  } else { //apply a multigrid cycle
    multigrid->Operator(o_rhs, o_x);
  }
}

void parAlmond_t::Report() {

  if(multigrid->comm.rank()==0) {
    printf("-----------------------------Multigrid Report-----------------------------------------------\n");
    printf("--------------------------------------------------------------------------------------------\n");
    printf("Level |    Type    |    Dimension   |  Per Rank Dim  |   nnz per row   |   Smoother        |\n");
    printf("      |            |                |  (min,max,avg) |  (min,max,avg)  |                   |\n");
    printf("--------------------------------------------------------------------------------------------\n");
  }

  for(int lev=0; lev<multigrid->numLevels-1; lev++) {
    if(multigrid->comm.rank()==0) {printf(" %3d  ", lev);fflush(stdout);}
    multigrid->levels[lev]->Report();
  }

  //base level
  multigrid->coarseSolver->Report(multigrid->numLevels-1);

  if(multigrid->comm.rank()==0)
    printf("--------------------------------------------------------------------------------------------\n");
}

int parAlmond_t::NumLevels() {
  return multigrid->numLevels;
}

dlong parAlmond_t::getNumCols(int k) {
  return multigrid->levels[k]->Ncols;
}

dlong parAlmond_t::getNumRows(int k) {
  return multigrid->levels[k]->Nrows;
}

} //namespace parAlmond

} //namespace libp
