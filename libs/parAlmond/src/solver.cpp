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

#include "solver.hpp"

namespace parAlmond {

solver_t::solver_t(occa::device otherdevice, MPI_Comm othercomm,
                   setupAide otheroptions) {

  device = otherdevice;

  comm = othercomm;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  levels = (multigridLevel **) calloc(MAX_LEVELS,sizeof(multigridLevel *));
  numLevels = 0;

  options = otheroptions;

  if (options.compareArgs("PARALMOND CYCLE", "NONSYM")) {
    ktype = GMRES;
  } else {
    ktype = PCG;
  }

  if(options.compareArgs("PARALMOND CYCLE", "EXACT"))
    exact = true;
  else
    exact = false;

  if(options.compareArgs("PARALMOND CYCLE", "VCYCLE"))
    ctype = VCYCLE;
  else
    ctype = KCYCLE;

  if (options.compareArgs("PARALMOND SMOOTHER", "CHEBYSHEV")) {
    stype = CHEBYSHEV;
    options.getArgs("PARALMOND CHEBYSHEV DEGREE", ChebyshevIterations);
    if (!ChebyshevIterations) ChebyshevIterations=2; //default to 2
  } else { //default to DAMPED_JACOBI
    stype = DAMPED_JACOBI;
  }
}

solver_t::~solver_t() {

  for (int n=0;n<numLevels;n++)
    delete levels[n];

  free(levels);
}

void solver_t::Report() {

  if(rank==0) {
    printf("------------------ParAlmond Report-----------------------------------\n");
    printf("---------------------------------------------------------------------\n");
    printf("level| active ranks |   dimension   |  nnzs         |  nnz/row      |\n");
    printf("     |              | (min,max,avg) | (min,max,avg) | (min,max,avg) |\n");
    printf("---------------------------------------------------------------------\n");
  }

  for(int lev=0; lev<numLevels; lev++) {
    printf(" %3d ", lev);
    levels[lev]->Report();
  }

  if(rank==0)
    printf("---------------------------------------------------------------------\n");
}

}