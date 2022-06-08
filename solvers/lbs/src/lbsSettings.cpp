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

#include "lbs.hpp"

//settings for lbs solver
lbsSettings_t::lbsSettings_t(comm_t& _comm):
  settings_t(_comm) {

  newSetting("DATA FILE",
             "data/lbsTaylorVortex2D.h",
             "Boundary and Initial conditions header");

  newSetting("VELOCITY MODEL",
             "D2Q9",
             "Microscopic Discrete Velocity Model",
             {"D2Q9", "D3Q15", "D3Q19"});

  newSetting("VISCOSITY",
             "1.0",
             "Fluid dynamic viscosity i.e. 1/Re");

  // newSetting("PML PROFILE ORDER",
  //            "4",
  //            "Polynomial degree of PML damping");

  // newSetting("PML SIGMAX MAX",
  //            "100",
  //            "Coefficent of polynomial PML damping in x direction");

  // newSetting("PML SIGMAY MAX",
  //            "100",
  //            "Coefficent of polynomial PML damping in y direction");

  // newSetting("PML SIGMAZ MAX",
  //            "100",
  //            "Coefficent of polynomial PML damping in z direction");

  // newSetting("PML INTEGRATION",
  //            "COLLOCATION",
  //            "Type of integration rule to PML damping profile",
  //            {"COLLOCATION", "CUBATURE"});

  newSetting("TIME INTEGRATOR",
             "LSERK4",
             "Time integration method",
             {"LSERK4"});

  newSetting("CFL NUMBER",
             "1.0",
             "Multiplier for timestep stability bound");

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
             "lbs");
}

void lbsSettings_t::report() {

  if (comm.rank()==0) {
    std::cout << "LBS Settings:\n\n";
    reportSetting("DATA FILE");
    // reportSetting("SPEED OF SOUND");
    reportSetting("VISCOSITY");
    // reportSetting("PML PROFILE ORDER");
    // reportSetting("PML SIGMAX MAX");
    // reportSetting("PML SIGMAY MAX");
    // reportSetting("PML SIGMAZ MAX");
    // reportSetting("PML INTEGRATION");
    reportSetting("TIME INTEGRATOR");
    reportSetting("START TIME");
    reportSetting("FINAL TIME");
    reportSetting("OUTPUT INTERVAL");
    reportSetting("OUTPUT TO FILE");
    reportSetting("OUTPUT FILE NAME");
  }
}

void lbsSettings_t::parseFromFile(platformSettings_t& platformSettings,
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
