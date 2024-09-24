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

#include "cns.hpp"

//settings for cns solver
cnsSettings_t::cnsSettings_t(comm_t& _comm):
  settings_t(_comm) {

  newSetting("DATA FILE", "data/cnsUniform2D.h",
             "Boundary and Initial conditions header");

  newSetting("SOLVER TYPE", "NAVIER-STOKES",
             "Solver to be used", 
            {"EULER", "NAVIER-STOKES"});

 newSetting("NONDIMENSIONAL EQUATIONS", "TRUE",
             "Compute R and Mu from Mach and Reynolds",
             {"FALSE", "TRUE"});

  newSetting("GAMMA", "1.4",
             "Specific heat ratio");

  newSetting("SPECIFIC GAS CONSTANT", "287.058",
             "Specific gas constant (if non-dimensional=FALSE)");

  newSetting("VISCOSITY", "0.0",
             "Bulk viscosity (if non-dimensional=TRUE)");

  newSetting("VISCOSITY TYPE", "CONSTANT",
             "Viscosity treatment", 
             {"CONSTANT", "SUTHERLAND", "POWER-LAW"});

  newSetting("MACH NUMBER", "0.1",
             "Mach number (if non-dimensional=TRUE)");

  newSetting("REYNOLDS NUMBER", "1000.0",
             "Reynolds number (if non-dimensional=TRUE)");

  newSetting("ANGLE OF ATTACK", "0.0",
             "Angle of attack (if non-dimensional=TRUE)");

  newSetting("PRANDTL NUMBER", "0.72",
             "Pranndtl Number");

   newSetting("ISOTHERMAL", "FALSE",
             "Solve isothermal equations", 
             {"FALSE", "TRUE"});

  newSetting("FLOW STATES", "1.0, 1.0, 0.0, 1.0, 1.0",
             "Reference states for problem setting");

   newSetting("GEOMETRIC-TO-PHYSICAL MAP",
             "Maps the geometric IDs to BC IDs");

  newSetting("IC STATE ID", "0",
             "State id defined in flow states to be used for initialization");

  newSetting("BC STATE ID", "0",
             "State id defined in flow states to be used for initialization");

  newSetting("ADVECTION TYPE", "COLLOCATION",
             "Integration type for flux terms",
             {"COLLOCATION", "CUBATURE"});

  newSetting("LDG BETA COEFFICIENT", "0.5",
             "Local DG method alternating upwinding coefficient between 0 and 1");

  newSetting("LDG TAU COEFFICIENT", "1.0",
             "Local DG jump stabilizer between 0  and 1");

  newSetting("TIME INTEGRATOR", "DOPRI5",
             "Time integration method",
             {"AB3", "DOPRI5", "LSERK4"});

  newSetting("CFL NUMBER", "1.0",
             "Multiplier for timestep stability bound");

  newSetting("START TIME", "0.0",
             "Start time for time integration");

  newSetting("FINAL TIME", "1.0",
             "End time for time integration");

  newSetting("OUTPUT INTERVAL", "0.1",
             "Time between printing output data");

  newSetting("OUTPUT TO FILE", "FALSE",
             "Flag for writing fields to VTU files",
             {"TRUE", "FALSE"});

  newSetting("OUTPUT FILE NAME", "output");

  newSetting("REPORT FORCES", "FALSE",
             "Flag for reporting forces",
             {"TRUE", "FALSE"});

  newSetting("REPORT COMPONENT", "FALSE",
             "Flag for reporting forces in component forms i.e. pressure, viscous",
             {"TRUE", "FALSE"});

  newSetting("REPORT GROUPS",
            "Geometric groups for reporting");

  newSetting("REPORT MOMENTS", "FALSE",
             "Flag for reporting moments",
             {"TRUE", "FALSE"});

  newSetting("MOMENT CENTER", "0.0, 0.0, 0.0",
             "Center for the moments");

  newSetting("GEOMETRIC GROUPS",
             "Groups of geometric entities in force/moment report");
 
  newSetting("STABILIZATION SOLVER TYPE", "1",
             "Solver type for stabilization",
             {"1", "2", "3"});

  newSetting("DETECTOR TYPE", "1",
             "Detector type for stabilization: NODETECT: 0 ALL: 1 KLOCKNER:2, PERSSON:3, DUCROS:4",
             {"0", "1", "2", "3", "4"});

  newSetting("STAB TYPE", "0",
             "Stabilization type: NOSTAB:0, FILTER:1, LIMITER:2, ART_DIFF:3, SUBCELL:4",
             {"0", "1", "2", "3", "4"});

  newSetting("ARTDIFF TYPE", "LAPLACE",
             "ARTIFICIAL DIFFUSION TYPE",
             {"LAPLACE", "PHYSICAL"});

  newSetting("ARTDIFF SCALE FACTOR", "3.0",
             "ARTIFICIAL DIFFUSION SCALING FACTOR");

  newSetting("FILTER CUTOFF", "1",
             "Exponential Filter Cutoff Order");

  newSetting("FILTER ORDER", "2",
             "Exponential Filter Order (must be even)");

  newSetting("SUBCELL NUMBER", "5",
             "Subcell number per edge (>=N)");

  newSetting("SUBCELL MINOR GRID", "WARPBLEND",
             "Subcell minor grid");

  newSetting("STAB OUTPUT TO FILE", "FALSE",
             "Detector Output to File ",
             {"TRUE", "FALSE"});
}



// stabSettings_t cnsSettings_t::extractStabSettings() {

//   stabSettings_t stabSettings(comm);

//   for(auto it = stabSettings.settings.begin(); it != stabSettings.settings.end(); ++it) {
//     setting_t& set = it->second;
//     const std::string name = set.getName();
//     std::string val;
//     getSetting(name, val);
//     set.updateVal(val);
//   }
//   return stabSettings;
// }

void cnsSettings_t::report() {

  if (comm.rank()==0) {
    std::cout << "CNS Settings:\n\n";

    reportSetting("DATA FILE");
    reportSetting("GAMMA");
    reportSetting("SPECIFIC GAS CONSTANT");
    reportSetting("VISCOSITY");
    reportSetting("VISCOSITY TYPE");
    reportSetting("SOLVER TYPE");
    reportSetting("MACH NUMBER");
    reportSetting("REYNOLDS NUMBER");
    
    reportSetting("FLOW STATES");
    reportSetting("IC STATE ID");
    reportSetting("BC STATE ID");

    reportSetting("ISOTHERMAL");
    reportSetting("ADVECTION TYPE");
    reportSetting("LDG BETA COEFFICIENT");
    reportSetting("LDG TAU COEFFICIENT");
    reportSetting("TIME INTEGRATOR");
    reportSetting("START TIME");
    reportSetting("FINAL TIME");
    reportSetting("OUTPUT INTERVAL");
    reportSetting("OUTPUT TO FILE");
    reportSetting("OUTPUT FILE NAME");
  }
}

void cnsSettings_t::parseFromFile(platformSettings_t& platformSettings,
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
