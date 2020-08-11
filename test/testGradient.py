#!/usr/bin/env python3

#####################################################################################
#
#The MIT License (MIT)
#
#Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#
#####################################################################################

from test import *

gradientData2D = gradientDir + "/data/gradientCos2D.h"
gradientData3D = gradientDir + "/data/gradientCos3D.h"

def gradientSettings(rcformat="2.0", data_file=gradientData2D,
                     mesh="BOX", dim=2, element=4, nx=10, ny=10, nz=10, boundary_flag=1,
                     degree=4, thread_model="Serial", platform_number=0, device_number=0):
  return [setting_t("FORMAT", rcformat),
          setting_t("DATA FILE", data_file),
          setting_t("MESH FILE", mesh),
          setting_t("MESH DIMENSION", dim),
          setting_t("ELEMENT TYPE", element),
          setting_t("BOX NX", nx),
          setting_t("BOX NY", ny),
          setting_t("BOX NZ", nz),
          setting_t("BOX BOUNDARY FLAG", boundary_flag),
          setting_t("POLYNOMIAL DEGREE", degree),
          setting_t("THREAD MODEL", thread_model),
          setting_t("PLATFORM NUMBER", platform_number),
          setting_t("DEVICE NUMBER", device_number),
          setting_t("OUTPUT TO FILE", "FALSE")]

if __name__ == "__main__":
  failCount=0;

  failCount += test(name="testGradientTri",
                    cmd=gradientBin,
                    settings=gradientSettings(element=3,data_file=gradientData2D,dim=2),
                    referenceNorm=4.44288293763069)

  failCount += test(name="testGradientQuad",
                    cmd=gradientBin,
                    settings=gradientSettings(element=4,data_file=gradientData2D,dim=2),
                    referenceNorm=4.44288293763069)

  failCount += test(name="testGradientTet",
                    cmd=gradientBin,
                    settings=gradientSettings(element=6,data_file=gradientData3D,dim=3),
                    referenceNorm=12.1674394705786)

  failCount += test(name="testGradientHex",
                    cmd=gradientBin,
                    settings=gradientSettings(element=12,data_file=gradientData3D,dim=3),
                    referenceNorm=12.1673360264757)

  sys.exit(failCount)