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

#include "wave.hpp"

//settings for wave solver
waveSettings_t::waveSettings_t(comm_t& _comm):
  settings_t(_comm) {

  newSetting("DATA FILE",
             "data/waveHomogeneous2D.h",
             "Boundary and Initial conditions header");

  newSetting("TIME INTEGRATOR",
             "ESDIRK6(5)9L[2]SA", 
             "Time integration method",
             {"TRBDF2-ESDIRK", 
                  "TRX2-ESDIRK", 
                  "ARK3(2)4L[2]SA-ESDIRK", 
                  "Kvaerno(4,2,3)-ESDIRK", 
                  "ESDIRK3(2)4L[2]SA", 
                  "ESDIRK3(2)5L[2]SA", 
                  "ESDIRK3(2I)5L[2]SA", 
                  "Kvaerno(5,3,4)-ESDIRK", 
                  "ESDIRK4(3)6L[2]SA", 
                  "ESDIRK4(3)7L[2]SA", 
                  "ESDIRK4(3I)6L[2]SA", 
                  "ARK4(3)6L[2]SA-ESDIRK", 
                  "ARK4(3)7L[2]SA-ESDIRK", 
                  "ESDIRK5(3)6L[2]SA", 
                  "ESDIRK5(4)7L[2]SA", 
                  "ESDIRK5(4)7L[2]SA2", 
                  "ESDIRK5(4)8L[2]SA", 
                  "ARK5(4)8L[2]SA-ESDIRK", 
                  "ARK5(4)8L[2]SAb-ESDIRK", 
                  "Kvaerno(7,4,5)-ESDIRK", 
                  "ESDIRK6(5)9L[2]SA", 
                  "ESDIRK5(4I)8L[2]SA", 
                  "ESDIRK6(4)7A[2]"});

  newSetting("START TIME",
             "0",
             "Start time for time integration");

  newSetting("FINAL TIME",
             "10",
             "End time for time integration");

  newSetting("TIME STEP",
             "0.1",
             "Time step");

  newSetting("OUTPUT INTERVAL",
             ".1",
             "Time between printing output data");

  newSetting("OUTPUT TO FILE",
             "FALSE",
             "Flag for writing fields to VTU files",
             {"TRUE", "FALSE"});

  newSetting("OUTPUT FILE NAME",
             "wave");

  newSetting("VERBOSE",
             "TRUE",
             "Verbose output setings for wave",
             {"TRUE", "FALSE"});
  
  ellipticAddSettings(*this, "ELLIPTIC ");
  parAlmond::AddSettings(*this, "ELLIPTIC ");
  InitialGuess::AddSettings(*this, "ELLIPTIC ");
  
}

void waveSettings_t::report() {

  if (comm.rank()==0) {
    std::cout << "WAVE Settings:\n\n";
    reportSetting("DATA FILE");
    reportSetting("TIME INTEGRATOR");
    reportSetting("START TIME");
    reportSetting("FINAL TIME");
    reportSetting("TIME STEP");
    reportSetting("OUTPUT INTERVAL");
    reportSetting("OUTPUT TO FILE");
    reportSetting("OUTPUT FILE NAME");

    std::cout << "\n Solver Settings:\n\n";

    reportSetting("ELLIPTIC DISCRETIZATION");
    reportSetting("ELLIPTIC LINEAR SOLVER");
    reportSetting("ELLIPTIC INITIAL GUESS STRATEGY");
    reportSetting("ELLIPTIC INITIAL GUESS HISTORY SPACE DIMENSION");
    reportSetting("ELLIPTIC PRECONDITIONER");
    
    if (compareSetting("ELLIPTIC PRECONDITIONER","MULTIGRID")) {
      reportSetting("ELLIPTIC MULTIGRID COARSENING");
      reportSetting("ELLIPTIC MULTIGRID SMOOTHER");
      if (compareSetting("ELLIPTIC MULTIGRID SMOOTHER","CHEBYSHEV"))
        reportSetting("ELLIPTIC MULTIGRID CHEBYSHEV DEGREE");
    }

    if (compareSetting("ELLIPTIC PRECONDITIONER","MULTIGRID")
      ||compareSetting("ELLIPTIC PRECONDITIONER","PARALMOND")) {
      reportSetting("ELLIPTIC PARALMOND CYCLE");
      reportSetting("ELLIPTIC PARALMOND SMOOTHER");
      reportSetting("ELLIPTIC PARALMOND CHEBYSHEV DEGREE");
    }
  }
}

void waveSettings_t::parseFromFile(platformSettings_t& platformSettings,
                                  meshSettings_t& meshSettings,
                                  const std::string filename) {
  //read all settings from file
  settings_t s(comm);
  s.readSettingsFromFile(filename);

  for(auto it = s.settings.begin(); it != s.settings.end(); ++it) {
    setting_t& set = it->second;
    const std::string name = set.getName();
    const std::string val = set.getVal<std::string>();
    if (platformSettings.hasSetting(name))
      platformSettings.changeSetting(name, val);
    else if (meshSettings.hasSetting(name))
      meshSettings.changeSetting(name, val);
    else if (hasSetting(name)) //self
      changeSetting(name, val);
    else  {
      LIBP_FORCE_ABORT("Unknown setting: [" << name << "] requested");
    }
  }
}


ellipticSettings_t waveSettings_t::extractEllipticSettings() {

  ellipticSettings_t  ellipticSettings(comm);

//  InitialGuess::AddSettings(ellipticSettings);

  for(auto it = ellipticSettings.settings.begin(); it != ellipticSettings.settings.end(); ++it) {
    setting_t& set = it->second;
    const std::string name = set.getName();

    std::string val;
    getSetting("ELLIPTIC "+name, val);

    set.updateVal(val);
  }

  return ellipticSettings;
}

