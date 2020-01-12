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

#include "settings.hpp"

void parAlmondAddSettings(settings_t& settings,
                          const string prefix) {

  settings.newSetting(prefix+"PARALMOND CYCLE",
                      "KCYCLE",
                      "Type of Multigrid Cycle",
                      {"VCYCLE", "KCYCLE"});

  settings.newSetting(prefix+"PARALMOND SMOOTHER",
                      "CHEBYSHEV",
                      "Type of Smoother",
                      {"DAMPEDJACOBI", "CHEBYSHEV"});

  settings.newSetting(prefix+"PARALMOND CHEBYSHEV DEGREE",
                      "2",
                      "Number of Chebyshev iteration to run in smoother");

  settings.newSetting(prefix+"PARALMOND PARTITION",
                      "STRONGNODES",
                      "Type of parallel node-distribution in coarse problems",
                      {"STRONGNODES", "DISTRIBUTED", "SATURATE"});

  settings.newSetting(prefix+"PARALMOND AGGREGATION STRATEGY",
                      "DEFAULT",
                      "Type of coarse node aggregation to use",
                      {"DEFAULT", "LPSCN"});

  settings.newSetting(prefix+"PARALMOND LPSCN ORDERING",
                      "NONE",
                      "Type of node ordering to use in LPSCN aggregation",
                      {"MIN", "MAX", "NONE"});

}

void parAlmondReportSettings(settings_t& settings) {

  settings.reportSetting("PARALMOND CYCLE");
  settings.reportSetting("PARALMOND SMOOTHER");

  if (settings.compareSetting("PARALMOND SMOOTHER","CHEBYSHEV"))
    settings.reportSetting("PARALMOND CHEBYSHEV DEGREE");

  settings.reportSetting("PARALMOND PARTITION");
  settings.reportSetting("PARALMOND AGGREGATION STRATEGY");

  if (settings.compareSetting("PARALMOND AGGREGATION STRATEGY", "LPSCN"))
    settings.reportSetting("PARALMOND LPSCN ORDERING");
}