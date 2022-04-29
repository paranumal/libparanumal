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

#include "parAlmond.hpp"

namespace libp {

namespace parAlmond {

void AddSettings(settings_t& settings,
                 const std::string prefix) {

  settings.newSetting(prefix+"PARALMOND CYCLE",
                      "VCYCLE",
                      "Type of Multigrid Cycle",
                      {"VCYCLE", "KCYCLE", "EXACT"});

  settings.newSetting(prefix+"PARALMOND STRENGTH",
                      "SYMMETRIC",
                      "Type of Alegraic Stength-of-Connection Measure",
                      {"RUGESTUBEN", "SYMMETRIC"});

  settings.newSetting(prefix+"PARALMOND AGGREGATION",
                      "SMOOTHED",
                      "Type of Prologation Operator",
                      {"SMOOTHED", "UNSMOOTHED"});

  settings.newSetting(prefix+"PARALMOND SMOOTHER",
                      "CHEBYSHEV",
                      "Type of Smoother",
                      {"DAMPEDJACOBI", "CHEBYSHEV"});

  settings.newSetting(prefix+"PARALMOND CHEBYSHEV DEGREE",
                      "2",
                      "Number of Chebyshev iteration to run in smoother");

}

void ReportSettings(settings_t& settings) {

  settings.reportSetting("PARALMOND CYCLE");
  settings.reportSetting("PARALMOND AGGREGATION");
  settings.reportSetting("PARALMOND SMOOTHER");

  if (settings.compareSetting("PARALMOND SMOOTHER","CHEBYSHEV"))
    settings.reportSetting("PARALMOND CHEBYSHEV DEGREE");
}

} //namespace parAlmond

} //namespace libp
