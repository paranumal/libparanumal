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

advectionData2D = advectionDir + "/data/advectionLinear2D.h"
advectionData3D = advectionDir + "/data/advectionLinear3D.h"

def advectionSettings(rcformat="2.0", data_file=advectionData2D,
                     mesh="BOX", dim=2, element=4, nx=10, ny=10, nz=10, boundary_flag=-1,
                     degree=4, thread_model="Serial", platform_number=0, device_number=0,
                     time_integrator="DOPRI5", start_time=0.0, final_time=1.0):
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
          setting_t("TIME INTEGRATOR", time_integrator),
          setting_t("START TIME", start_time),
          setting_t("FINAL TIME", final_time),
          setting_t("OUTPUT TO FILE", "FALSE")]

if __name__ == "__main__":
  failCount=0;

  failCount += test(name="testAdvectionTri",
                    cmd=advectionBin,
                    settings=advectionSettings(element=3,data_file=advectionData2D,dim=2),
                    referenceNorm=0.723897834264616)

  failCount += test(name="testAdvectionQuad",
                    cmd=advectionBin,
                    settings=advectionSettings(element=4,data_file=advectionData2D,dim=2),
                    referenceNorm=0.722791610885232)

  failCount += test(name="testAdvectionTet",
                    cmd=advectionBin,
                    settings=advectionSettings(element=6,data_file=advectionData3D,dim=3),
                    referenceNorm=0.835474966415895)

  failCount += test(name="testAdvectionHex",
                    cmd=advectionBin,
                    settings=advectionSettings(element=12,data_file=advectionData3D,dim=3),
                    referenceNorm=0.833820360927384)

  failCount += test(name="testAdvectionTri_MPI", ranks=4,
                    cmd=advectionBin,
                    settings=advectionSettings(element=3,data_file=advectionData2D,dim=2),
                    referenceNorm=0.723627520020827)


  sys.exit(failCount)