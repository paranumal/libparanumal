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

#include "elliptic.hpp"

//settings for elliptic solver
ellipticSettings_t::ellipticSettings_t(MPI_Comm& _comm):
  settings_t(_comm) {

  occaAddSettings(*this);
  meshAddSettings(*this);


  newSetting("DATA FILE",
             "data/ellipticSine2D.h",
             "Boundary and Initial conditions header");

  newSetting("LAMBDA",
             "1.0",
             "Coefficient in Screened Poisson Equation");

  newSetting("COEFFICIENT",
             "CONSTANT",
             "Coefficient in Screened Poisson Operator",
             {"CONSTANT", "VARIABLE"});

  newSetting("DISCRETIZATION",
             "CONTINUOUS",
             "Type of Finite Element Discretization",
             {"CONTINUOUS", "IPDG"});

  newSetting("LINEAR SOLVER",
             "PCG",
             "Iterative Linear Solver to use for solve",
             {"PCG", "PCG,FLEXIBLE"});

  newSetting("PRECONDITIONER",
             "NONE",
             "Preconditioning Strategy",
             {"NONE", "JACOBI", "MASSMATRIX", "FULLALMOND", "MULTIGRID", "SEMFEM"});

  /* MULTIGRID options */
  newSetting("MULTIGRID COARSENING",
             "HALFDOFS",
             "p-Multigrid coarsening strategy",
             {"ALLDEGREES", "HALFDEGREES", "HALFDOFS"});

  newSetting("MULTIGRID SMOOTHER",
             "CHEBYSHEV",
             "p-Multigrid smoother",
             {"DAMPEDJACOBI", "CHEBYSHEV"});

  newSetting("MULTIGRID CHEBYSHEV DEGREE",
             "2",
             "Smoothing iterations in Chebyshev smoother");

  newSetting("OUTPUT TO FILE",
             "FALSE",
             "Flag for writing fields to VTU files",
             {"TRUE", "FALSE"});

  newSetting("OUTPUT FILE NAME",
             "elliptic");

  newSetting("VERBOSE",
             "FALSE",
             "Enable verbose output",
             {"TRUE", "FALSE"});

  parAlmondAddSettings(*this);
}

void ellipticSettings_t::report() {

  std::cout << "Settings:\n\n";
  occaReportSettings(*this);
  meshReportSettings(*this);

  reportSetting("DATA FILE");

  reportSetting("LAMBDA");
  reportSetting("COEFFICIENT");
  reportSetting("DISCRETIZATION");
  reportSetting("LINEAR SOLVER");
  reportSetting("PRECONDITIONER");

  if (compareSetting("PRECONDITIONER","MULTIGRID")) {
    reportSetting("MULTIGRID COARSENING");
    reportSetting("MULTIGRID SMOOTHER");
    if (compareSetting("MULTIGRID SMOOTHER","CHEBYSHEV"))
      reportSetting("MULTIGRID CHEBYSHEV DEGREE");
  }

  if (compareSetting("PRECONDITIONER","MULTIGRID")
    ||compareSetting("PRECONDITIONER","FULLALMOND"))
    parAlmondReportSettings(*this);

  reportSetting("OUTPUT TO FILE");
  reportSetting("OUTPUT FILE NAME");
}