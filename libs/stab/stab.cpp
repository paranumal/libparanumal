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

#include "stab.hpp"

// using namespace libp;
namespace libp {

void stab_t::Setup(platform_t& _platform, mesh_t &_mesh, 
                  stabSettings_t& _settings, properties_t& kernelInfo){

  platform        = _platform;
  settings        = _settings;
  mesh            = _mesh;
  comm            = _mesh.comm;  // copy base mesh communicator
  props           =  kernelInfo; 

  // Set Types 
  int solverType = 0, detectorType=0, stabType=0; 
  settings.getSetting("STABILIZATION SOLVER TYPE", solverType);
  settings.getSetting("DETECTOR TYPE", detectorType);
  settings.getSetting("STAB TYPE", stabType);

  // Set solver, detector and stabilizer type
  setTypes(Stab::Solver(solverType), 
          Stab::Detector(detectorType), 
          Stab::Type(stabType));

  //Trigger JIT kernel builds
  // ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Max);
  // ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Mul);
  // ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Min);

  platform.linAlg().InitKernels({"innerProd", "max", "sum", "set", "amx"});

  
  // Setup Detector
  detectSetup();

  // Setup Stabilizer
  stabSetup();

  // Copy kernel information as the solver needs stab definitions
  kernelInfo = props;
};

} //namespace libp
