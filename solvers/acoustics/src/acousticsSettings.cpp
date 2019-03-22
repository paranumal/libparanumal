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

#include "acoustics.h"

//settings for acoustics solver
acousticsSettings_t::acousticsSettings_t(): settings_t() {

  occaAddSettings(*this);
  meshAddSettings(*this);

  this->newSetting("TIME INTEGRATOR",
                      "DOPRI5",
                      "Time integration method",
                      {"DOPRI5", "LSERK4"});
  this->newSetting("COMPUTE ERROR FLAG",
                      "1");

  this->newSetting("TSTEPS FOR ERROR COMPUTE",
                      "1000",
                      "Number of time steps between error check");

  this->newSetting("TSTEPS FOR SOLUTION OUTPUT",
                      "1000",
                      "Number of time steps between output frames");

  this->newSetting("REPORT FREQUENCY",
                      "100",
                      "Number of time steps between reporting");

  this->newSetting("START TIME",
                      "0");

  this->newSetting("FINAL TIME",
                      "10");

  this->newSetting("OUTPUT INTERVAL",
                      ".1");

  this->newSetting("MAX MRAB LEVELS",
                      "1");

  this->newSetting("FORMAT",
                      "2.0",
                      "Setting file version number",
                      {"2.0"});
}

// void acousticsSettings_t::readSettingsFromFile(const string filename) {
//   this->settings_t::readSettingsFromFile(filename);
// }
