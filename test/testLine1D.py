#!/usr/bin/env python3

#####################################################################################
#
#The MIT License (MIT)
#
#Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus
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

from testTW import *

ellipticData1D = ellipticDir + "/data/ellipticHorns1D.h"
ellipticData2D = ellipticDir + "/data/ellipticSine2D.h"
ellipticData3D = ellipticDir + "/data/ellipticSine3D.h"

def ellipticSettings(rcformat="2.0", data_file=ellipticData1D,
                     mesh="BOX", dim=1, element=2, nx=10, ny=10, nz=10, boundary_flag=1,
                     degree=4, thread_model="CUDA", platform_number=0, device_number=0,
                     Lambda=1.0,
                     discretization="CONTINUOUS",
                     linear_solver="PCG",
                     precon="MULTIGRID",
                     multigrid_smoother="CHEBYSHEV",
                     paralmond_cycle="VCYCLE",
                     paralmond_strength="SYMMETRIC",
                     paralmond_aggregation="SMOOTHED",
                     paralmond_smoother="CHEBYSHEV",
                     output_to_file="FALSE"):
  return [setting_t("FORMAT", rcformat),
          setting_t("DATA FILE", data_file),
          setting_t("MESH FILE", mesh),
          setting_t("MESH DIMENSION", dim),
          setting_t("ELEMENT TYPE", element),
          setting_t("BOX NX", nx),
          setting_t("BOX NY", ny),
          setting_t("BOX NZ", nz),
          setting_t("BOX DIMX", 2),
          setting_t("BOX DIMY", 1),
          setting_t("BOX DIMZ", 1),
          setting_t("BOX BOUNDARY FLAG", boundary_flag),
          setting_t("POLYNOMIAL DEGREE", degree),
          setting_t("THREAD MODEL", thread_model),
          setting_t("PLATFORM NUMBER", platform_number),
          setting_t("DEVICE NUMBER", device_number),
          setting_t("DISCRETIZATION", discretization),
          setting_t("LINEAR SOLVER", linear_solver),
          setting_t("PRECONDITIONER", precon),
          setting_t("MULTIGRID SMOOTHER", multigrid_smoother),
          setting_t("MULTIGRID CHEBYSHEV DEGREE", 2),
          setting_t("PARALMOND CYCLE", paralmond_cycle),
          setting_t("PARALMOND STRENGTH", paralmond_strength),
          setting_t("PARALMOND AGGREGATION", paralmond_aggregation),
          setting_t("PARALMOND SMOOTHER", paralmond_smoother),
          setting_t("PARALMOND CHEBYSHEV DEGREE", 2),
          setting_t("OUTPUT TO FILE", "FALSE"),
          setting_t("VERBOSE", output_to_file)]

def main():
  failCount=0;

  ################
  #C0 tests
  ################

  for L in range(1, 21, 1):
    for P in range(1,10,1):
      nx = 2**L
      degree = P
      failCount += testTW(name="testLine1Da",
                          cmd=ellipticBin,
                          settings=ellipticSettings(element=2,
                                                    degree=degree,
                                                    nx=nx,
                                                    data_file=ellipticData1D,
                                                    dim=1,
                                                    precon="MULTIGRID"))                    


  #clean up
  for file_name in os.listdir(testDir):
    if file_name.endswith('.vtu'):
      os.remove(testDir + "/" + file_name)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)
