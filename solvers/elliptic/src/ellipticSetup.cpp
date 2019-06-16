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

elliptic_t::elliptic_t(mesh_t& _mesh): solver_t(_mesh) {}

elliptic_t& elliptic_t::Setup(mesh_t& mesh){

  elliptic_t* elliptic = new elliptic_t(mesh);

  settings_t& settings = elliptic->settings;

  elliptic->Nfields = 1;

  elliptic->disc_ipdg = settings.compareSetting("DISCRETIZATION","IPDG");
  elliptic->disc_c0   = settings.compareSetting("DISCRETIZATION","CONTINUOUS");

  //setup linear algebra module
  elliptic->linAlg = linAlg_t::Setup(elliptic->device, settings, elliptic->props);
  elliptic->linAlg->InitKernels({"add", "sum",
                                 "axpy", "zaxpy",
                                 "axmy", "zaxmy",
                                 "axdy", "zaxdy",
                                 "innerProd", "weightedInnerProd",
                                 "norm2", "weightedNorm2"});

  //setup linear solver
  elliptic->linearSolver = linearSolver_t::Setup(*elliptic);
  elliptic->linearSolver->Init();

  // Boundary Type translation. Just defaults.
  int BCType[3] = {0,1,2};
  elliptic->BCType = (int*) calloc(3,sizeof(int));
  memcpy(elliptic->BCType,BCType,3*sizeof(int));

  return elliptic;
}
