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

#include "fpe.hpp"

//settings for fpe solver
fpeSettings_t::fpeSettings_t(MPI_Comm& _comm):
  settings_t(_comm) {

  newSetting("DATA FILE",
             "data/fpeLinear2D.h",
             "Boundary and Initial conditions header");

  newSetting("VISCOSITY",
             "1.0",
             "Diffusion strength");

  newSetting("ADVECTION TYPE",
             "COLLOCATION",
             "Integration type for flux terms",
             {"COLLOCATION", "CUBATURE"});

  newSetting("TIME INTEGRATOR",
             "DOPRI5",
             "Time integration method",
             {"AB3", "DOPRI5", "LSERK4", "EXTBDF3", "SSBDF3"});

  newSetting("NUMBER OF SUBCYCLES",
             "1",
             "Ratio of full timestep size to subcycling step size");

  newSetting("SUBCYCLING TIME INTEGRATOR",
             "DOPRI5",
             "Time integration method used in subcycling",
             {"AB3", "DOPRI5", "LSERK4"});

  newSetting("START TIME",
             "0",
             "Start time for time integration");

  newSetting("FINAL TIME",
             "10",
             "End time for time integration");

  newSetting("OUTPUT INTERVAL",
             ".1",
             "Time between printing output data");

  newSetting("OUTPUT TO FILE",
             "FALSE",
             "Flag for writing fields to VTU files",
             {"TRUE", "FALSE"});

  newSetting("OUTPUT FILE NAME",
             "fpe");

  ellipticAddSettings(*this, "ELLIPTIC ");
  parAlmondAddSettings(*this, "ELLIPTIC ");
}

void fpeSettings_t::report() {

  int rank;
  MPI_Comm_rank(comm, &rank);

  if (rank==0) {
    std::cout << "Fokker Planck Settings:\n\n";
    reportSetting("DATA FILE");
    reportSetting("TIME INTEGRATOR");

    if (compareSetting("TIME INTEGRATOR","SSBDF3")) {
      reportSetting("NUMBER OF SUBCYCLES");
      reportSetting("SUBCYCLING TIME INTEGRATOR");
    }

    reportSetting("START TIME");
    reportSetting("FINAL TIME");
    reportSetting("OUTPUT INTERVAL");
    reportSetting("OUTPUT TO FILE");
    reportSetting("OUTPUT FILE NAME");

    std::cout << "\nElliptic Solver Settings:\n\n";

    reportSetting("ELLIPTIC DISCRETIZATION");
    reportSetting("ELLIPTIC LINEAR SOLVER");
    reportSetting("ELLIPTIC PRECONDITIONER");

    if (compareSetting("ELLIPTIC PRECONDITIONER","MULTIGRID")) {
      reportSetting("ELLIPTIC MULTIGRID COARSENING");
      reportSetting("ELLIPTIC MULTIGRID SMOOTHER");
      if (compareSetting("ELLIPTIC MULTIGRID SMOOTHER","CHEBYSHEV"))
        reportSetting("ELLIPTIC MULTIGRID CHEBYSHEV DEGREE");
    }

    if (compareSetting("ELLIPTIC PRECONDITIONER","MULTIGRID")
      ||compareSetting("ELLIPTIC PRECONDITIONER","FULLALMOND")) {
      reportSetting("ELLIPTIC PARALMOND CYCLE");
      reportSetting("ELLIPTIC PARALMOND SMOOTHER");
      reportSetting("ELLIPTIC PARALMOND CHEBYSHEV DEGREE");
    }
  }
}

void fpeSettings_t::parseFromFile(occaSettings_t& occaSettings,
                                  meshSettings_t& meshSettings,
                                  const string filename) {
  //read all settings from file
  settings_t s(comm);
  s.readSettingsFromFile(filename);

  for(auto it = s.settings.begin(); it != s.settings.end(); ++it) {
    setting_t* set = it->second;
    const string name = set->getName();
    const string val = set->getVal<string>();
    if (occaSettings.hasSetting(name))
      occaSettings.changeSetting(name, val);
    else if (meshSettings.hasSetting(name))
      meshSettings.changeSetting(name, val);
    else if (hasSetting(name)) //self
      changeSetting(name, val);
    else  {
      stringstream ss;
      ss << "Unknown setting: [" << name << "] requested";
      LIBP_ABORT(ss.str());
    }
  }
}

ellipticSettings_t* fpeSettings_t::extractEllipticSettings() {

  ellipticSettings_t* ellipticSettings = new ellipticSettings_t(comm);

  for(auto it = ellipticSettings->settings.begin(); it != ellipticSettings->settings.end(); ++it) {
    setting_t* set = it->second;
    const string name = set->getName();

    string val;
    getSetting("ELLIPTIC "+name, val);

    set->updateVal(val);
  }

  return ellipticSettings;
}