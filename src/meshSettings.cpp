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

#include "mesh.h"

void meshAddSettings(settings_t& settings) {

  settings.newSetting("MESH FILE",
                      "Box"
                      "Input mesh");
  settings.newSetting("MESH DIMENSION",
                      "2",
                      "Dimension of input mesh",
                      {"2","3"});
  settings.newSetting("ELEMENT TYPE",
                      "3",
                      "Type of mesh elements (number of edges)",
                      {"3","4","6","12"});
  settings.newSetting("POLYNOMIAL DEGREE",
                      "4",
                      "Degree of polynomial finite element space",
                      {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"});
}

