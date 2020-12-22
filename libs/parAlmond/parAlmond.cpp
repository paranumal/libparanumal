/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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
#include "parAlmond/parAlmondMultigrid.hpp"
#include "parAlmond/parAlmondKernels.hpp"

namespace parAlmond {

parAlmond_t::parAlmond_t(platform_t& _platform, settings_t& _settings, MPI_Comm comm):
  platform(_platform), settings(_settings) {

  platform.linAlg.InitKernels({"set", "add", "sum", "scale",
                                "axpy", "zaxpy",
                                "amx", "amxpy", "zamxpy",
                                "adx", "adxpy", "zadxpy",
                                "innerProd", "weightedInnerProd",
                                "norm2", "weightedNorm2"});

  multigrid = new multigrid_t(platform, settings, comm);

  //build parAlmond kernels on first construction
  if (Nrefs==0) buildParAlmondKernels(platform);
  Nrefs++;
}

void parAlmond_t::Operator(occa::memory& o_rhs, occa::memory& o_x) {

  if (multigrid->exact){ //call the linear solver
    int maxIter = 500;
    int verbose = settings.compareSetting("VERBOSE", "TRUE") ? 1 : 0;
    dfloat tol = 1e-8;
    solver_t &A = *(multigrid->levels[0]);
    (void) multigrid->linearSolver->Solve(A, *multigrid, o_x, o_rhs, tol, maxIter, verbose);
  } else { //apply a multigrid cycle
    multigrid->Operator(o_rhs, o_x);
  }
}

//Add level to multigrid heirarchy
void parAlmond_t::AddLevel(multigridLevel* level) {
  multigrid->AddLevel(level);
}

void parAlmond_t::Report() {

  int rank;
  MPI_Comm_rank(multigrid->comm, &rank);

  if(rank==0) {
    printf("-----------------------------Multigrid Report-----------------------------------------------\n");
    printf("--------------------------------------------------------------------------------------------\n");
    printf("Level |    Type    |    Dimension   |  Per Rank Dim  |   nnz per row   |   Smoother        |\n");
    printf("      |            |                |  (min,max,avg) |  (min,max,avg)  |                   |\n");
    printf("--------------------------------------------------------------------------------------------\n");
  }

  for(int lev=0; lev<multigrid->numLevels-1; lev++) {
    if(rank==0) {printf(" %3d  ", lev);fflush(stdout);}
    multigrid->levels[lev]->Report();
  }

  //base level
  multigrid->coarseSolver->Report(multigrid->numLevels-1);

  if(rank==0)
    printf("--------------------------------------------------------------------------------------------\n");
}

dlong parAlmond_t::getNumCols(int k) {
  return multigrid->levels[k]->Ncols;
}

dlong parAlmond_t::getNumRows(int k) {
  return multigrid->levels[k]->Nrows;
}

parAlmond_t::~parAlmond_t() {
  Nrefs--;
  if (Nrefs==0) freeParAlmondKernels();

  delete multigrid;
}

} //namespace parAlmond
