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

#include "mesh.hpp"
#include "parAdogs.hpp"

namespace libp {

meshSettings_t::meshSettings_t(comm_t _comm):
  settings_t(_comm) {

  newSetting("MESH FILE",
             "BOX"
             "Input mesh. Can set to BOX to instead construct a uniform mesh");
  newSetting("MESH DIMENSION",
             "2",
             "Dimension of input mesh",
             {"2","3"});
  newSetting("ELEMENT TYPE",
             "3",
             "Type of mesh elements (number of edges)",
             {"3","4","6","12"});
  newSetting("ELEMENT MAP",
             "ISOPARAMETRIC",
             "Type mapping used to transform each element",
             {"ISOPARAMETRIC","AFFINE"});

  newSetting("BOX DIMX",
             "10",
             "Length of BOX domain in X-dimension");
  newSetting("BOX DIMY",
             "10",
             "Length of BOX domain in Y-dimension");
  newSetting("BOX DIMZ",
             "10",
             "Length of BOX domain in Z-dimension");

  newSetting("BOX NX",
             "10",
             "Number of elements in X-dimension per rank");
  newSetting("BOX NY",
             "10",
             "Number of elements in Y-dimension per rank");
  newSetting("BOX NZ",
             "10",
             "Number of elements in Z-dimension per rank");

  newSetting("BOX GLOBAL NX",
             "0",
             "Global number of elements in X-dimension per rank");
  newSetting("BOX GLOBAL NY",
             "0",
             "Global number of elements in Y-dimension per rank");
  newSetting("BOX GLOBAL NZ",
             "0",
             "Global number of elements in Z-dimension per rank");

  newSetting("BOX BOUNDARY FLAG",
             "1",
             "Type of boundary conditions for BOX domain (-1 for periodic)");

  newSetting("POLYNOMIAL DEGREE",
             "4",
             "Degree of polynomial finite element space",
             {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"});

  paradogs::AddSettings(*this);
}

void meshSettings_t::report() {

  if (comm.rank()==0) {
    std::cout << "Mesh Settings:\n\n";
    if (!compareSetting("MESH FILE","BOX"))
      reportSetting("MESH FILE");

    reportSetting("MESH DIMENSION");
    reportSetting("ELEMENT TYPE");

    if (compareSetting("ELEMENT TYPE","4") ||
        compareSetting("ELEMENT TYPE","12"))
      reportSetting("ELEMENT MAP");

    //report the box settings
    if (compareSetting("MESH FILE","BOX")) {
      reportSetting("BOX DIMX");
      reportSetting("BOX DIMY");
      reportSetting("BOX DIMZ");

      dlong NX, NY, NZ;
      getSetting("BOX GLOBAL NX", NX);
      getSetting("BOX GLOBAL NY", NY);
      getSetting("BOX GLOBAL NZ", NZ);

      if (NX*NY*NZ > 0) {
        reportSetting("BOX GLOBAL NX");
        reportSetting("BOX GLOBAL NY");
        reportSetting("BOX GLOBAL NZ");
      } else {//global sizes not provided, report local size
        reportSetting("BOX NX");
        reportSetting("BOX NY");
        reportSetting("BOX NZ");
      }

      reportSetting("BOX BOUNDARY FLAG");
    }

    reportSetting("POLYNOMIAL DEGREE");

    if (!compareSetting("MESH FILE","BOX")) {
      paradogs::ReportSettings(*this);
    }
  }
}

} //namespace libp
