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

#include "esdg.hpp"

//settings for esdg solver
esdgSettings_t::esdgSettings_t(comm_t _comm):
  settings_t(_comm) {

  newSetting("DATA FILE",
             "data/esdgTalorVortex2D.h",
             "Boundary and Initial conditioesdg header");

  newSetting("GAMMA",
             "1.4",
             "Specific heat ratio");

  newSetting("MEAN FLOW DENSITY",
             "1",
             "Density of mean flow");

  newSetting("MEAN FLOW XVELOCITY",
             "1",
	     "X component of mean flow velocity");

  newSetting("MEAN FLOW YVELOCITY",
             "0",
	     "Y component of mean flow velocity");

  // mean flow mach number will be = sqrt(ubar*ubar+vbar*vbar)/sqrt(gamma*pbar/rbar)
  newSetting("MEAN FLOW PRESSURE",
             "4",
             "Pressure of mean flow");

  
  newSetting("ADVECTION TYPE",
             "COLLOCATION",
             "Integration type for flux terms",
             {"COLLOCATION", "CUBATURE", "ESDG", "SBP", "CUBDG"});

  newSetting("TIME INTEGRATOR",
             "DOPRI5",
             "Time integration method",
             {"AB3", "DOPRI5", "LSERK4"});

  newSetting("START TIME",
             "0",
             "Start time for time integration");

  newSetting("FINAL TIME",
             "10",
             "End time for time integration");

  
  newSetting("DIFFUSION ONLY TIME",
             "0",
             "Time to run with diffusion only");
  
  newSetting("RELAXATION TIME",
             "0",
             "Time to run with relaxation filter based diffusion");
  
  newSetting("OUTPUT INTERVAL",
             ".1",
             "Time between printing output data");

  newSetting("OUTPUT TO FILE",
             "FALSE",
             "Flag for writing fields to VTU files",
             {"TRUE", "FALSE"});

  newSetting("OUTPUT FILE NAME",
             "esdg");


  newSetting("FLUX DEGREE INCREMENT",
             "0",
             "Degree of flux node set over N");

  newSetting("FILTER TOP MODES",
             "0",
             "Number of top nodes to filter");


  newSetting("LAME LAMBDA",
             "0.0",
             "Lambda Lame parameter");

  newSetting("LAME MU",
             "0.0",
             "Mu Lame parameter");

  newSetting("OUTFLOW MU",
             "0.0",
             "Value of mu Lame parameter at outflow (0 for no outflow ramp)");

  newSetting("OUTFLOW LAYER XMIN",
             "0.0",
             "X-coordinate of start of outflow viscosity ramp");

  newSetting("OUTFLOW LAYER XMAX",
             "0.0",
             "X-coordinate of end of outflow viscosity ramp");


  newSetting("CHECKPOINT SAVE",
             "FALSE",
             "Flag to specify if checkpoint files are saved",
             {"TRUE", "FALSE"});


  newSetting("CHECKPOINT SAVE NAME",
             "default",
             "Root of checkpoint output files");

  newSetting("CHECKPOINT LOAD",
             "FALSE",
             "Flag to specify if checkpoint files are loaded",
             {"TRUE", "FALSE"});
  
  newSetting("CHECKPOINT LOAD NAME",
             "default",
             "Root of checkpoint input files");

  newSetting("GUARD LEVEL",
             "0",
             "Mask for choosing operations to guard");


}

void esdgSettings_t::report() {

  int rank = comm.rank();

  if (rank==0) {
    std::cout << "ESDG Settings:\n\n";

    reportSetting("DATA FILE");

    // flow quantities
    reportSetting("GAMMA");
    reportSetting("MEAN FLOW DENSITY");
    reportSetting("MEAN FLOW XVELOCITY");
    reportSetting("MEAN FLOW YVELOCITY");
    reportSetting("MEAN FLOW PRESSURE");
    reportSetting("LAME MU");
    reportSetting("LAME LAMBDA");

    reportSetting("OUTFLOW MU");
    reportSetting("OUTFLOW LAYER XMIN");
    reportSetting("OUTFLOW LAYER XMAX");

    // discretization types
    reportSetting("ADVECTION TYPE");
    reportSetting("TIME INTEGRATOR");

    // time info
    reportSetting("START TIME");
    reportSetting("FINAL TIME");

    // OUTPUT info
    reportSetting("OUTPUT INTERVAL");
    reportSetting("OUTPUT TO FILE");
    reportSetting("OUTPUT FILE NAME");

    // settings 
    reportSetting("FLUX DEGREE INCREMENT");
    reportSetting("FILTER TOP MODES");

    // checkpoint 
    reportSetting("CHECKPOINT SAVE");
    reportSetting("CHECKPOINT SAVE NAME");

    reportSetting("CHECKPOINT LOAD");
    reportSetting("CHECKPOINT LOAD NAME");

    reportSetting("GUARD LEVEL");

  }
}

void esdgSettings_t::parseFromFile(platformSettings_t& platformSettings,
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
      std::stringstream ss;
      ss << "Unknown setting: [" << name << "] requested";
      LIBP_FORCE_ABORT(ss.str());
    }
  }
}
