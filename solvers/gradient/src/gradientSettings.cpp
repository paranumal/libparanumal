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

#include "gradient.hpp"

//settings for gradient solver
gradientSettings_t::gradientSettings_t(MPI_Comm& _comm):
  settings_t(_comm) {

  newSetting("DATA FILE",
             "data/gradientCos2D.h",
             "Initial conditions header");

  newSetting("OUTPUT TO FILE",
             "FALSE",
             "Flag for writing fields to VTU files",
             {"TRUE", "FALSE"});
}

void gradientSettings_t::report() {

  int rank;
  MPI_Comm_rank(comm, &rank);

  if (rank==0) {
    std::cout << "Gradient Settings:\n\n";
    reportSetting("DATA FILE");
    reportSetting("OUTPUT TO FILE");
  }
}

void gradientSettings_t::parseFromFile(occaSettings_t& occaSettings,
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