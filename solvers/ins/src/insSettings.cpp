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

#include "ins.hpp"

//settings for ins solver
insSettings_t::insSettings_t(MPI_Comm& _comm):
  settings_t(_comm) {

  newSetting("DATA FILE",
             "data/insLinear2D.h",
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
             "ins");

  ellipticAddSettings(*this, "VELOCITY ");
  parAlmondAddSettings(*this, "VELOCITY ");

  ellipticAddSettings(*this, "PRESSURE ");
  parAlmondAddSettings(*this, "PRESSURE ");
}

void insSettings_t::report() {

  int rank;
  MPI_Comm_rank(comm, &rank);

  if (rank==0) {
    std::cout << "INS Settings:\n\n";
    reportSetting("DATA FILE");
    reportSetting("ADVECTION TYPE");
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

    std::cout << "\nVelocity Solver Settings:\n\n";

    reportSetting("VELOCITY DISCRETIZATION");
    reportSetting("VELOCITY LINEAR SOLVER");
    reportSetting("VELOCITY PRECONDITIONER");

    if (compareSetting("VELOCITY PRECONDITIONER","MULTIGRID")) {
      reportSetting("VELOCITY MULTIGRID COARSENING");
      reportSetting("VELOCITY MULTIGRID SMOOTHER");
      if (compareSetting("VELOCITY MULTIGRID SMOOTHER","CHEBYSHEV"))
        reportSetting("VELOCITY MULTIGRID CHEBYSHEV DEGREE");
    }

    if (compareSetting("VELOCITY PRECONDITIONER","MULTIGRID")
      ||compareSetting("VELOCITY PRECONDITIONER","FULLALMOND")) {
      reportSetting("VELOCITY PARALMOND CYCLE");
      reportSetting("VELOCITY PARALMOND SMOOTHER");
      reportSetting("VELOCITY PARALMOND CHEBYSHEV DEGREE");
    }

    std::cout << "\nPressure Solver Settings:\n\n";

    reportSetting("PRESSURE DISCRETIZATION");
    reportSetting("PRESSURE LINEAR SOLVER");
    reportSetting("PRESSURE PRECONDITIONER");

    if (compareSetting("PRESSURE PRECONDITIONER","MULTIGRID")) {
      reportSetting("PRESSURE MULTIGRID COARSENING");
      reportSetting("PRESSURE MULTIGRID SMOOTHER");
      if (compareSetting("PRESSURE MULTIGRID SMOOTHER","CHEBYSHEV"))
        reportSetting("PRESSURE MULTIGRID CHEBYSHEV DEGREE");
    }

    if (compareSetting("PRESSURE PRECONDITIONER","MULTIGRID")
      ||compareSetting("PRESSURE PRECONDITIONER","FULLALMOND")) {
      reportSetting("PRESSURE PARALMOND CYCLE");
      reportSetting("PRESSURE PARALMOND SMOOTHER");
      reportSetting("PRESSURE PARALMOND CHEBYSHEV DEGREE");
    }
  }
}

void insSettings_t::parseFromFile(occaSettings_t& occaSettings,
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

ellipticSettings_t* insSettings_t::extractVelocitySettings() {

  ellipticSettings_t* velocitySettings = new ellipticSettings_t(comm);

  for(auto it = velocitySettings->settings.begin(); it != velocitySettings->settings.end(); ++it) {
    setting_t* set = it->second;
    const string name = set->getName();

    string val;
    getSetting("VELOCITY "+name, val);

    set->updateVal(val);
  }

  return velocitySettings;
}

ellipticSettings_t* insSettings_t::extractPressureSettings() {

  ellipticSettings_t* pressureSettings = new ellipticSettings_t(comm);

  for(auto it = pressureSettings->settings.begin(); it != pressureSettings->settings.end(); ++it) {
    setting_t* set = it->second;
    const string name = set->getName();

    string val;
    getSetting("PRESSURE "+name, val);

    set->updateVal(val);
  }

  return pressureSettings;
}