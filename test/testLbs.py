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

lbsData2D = lbsDir + "/data/lbsGaussian2D.h"
lbsData3D = lbsDir + "/data/lbsGaussian3D.h"

def lbsSettings(rcformat="2.0", data_file=lbsData2D,
               mesh="BOX", dim=2, element=4, nx=10, ny=10, nz=10, boundary_flag=1,
               degree=4, thread_model=device, platform_number=0, device_number=0,
               viscosity=0.01, velmodel="D2Q9", time_integrator="LSERK4", cfl=0.5, start_time=0.0, final_time=10.0,
               output_to_file="FALSE"):
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
          setting_t("VISCOSITY", viscosity),
          setting_t("VELOCITY MODEL", velmodel),
          # setting_t("SPEED OF SOUND", speed_of_sound),
          # setting_t("PML PROFILE ORDER", pml_order),
          # setting_t("PML SIGMAX MAX", pml_sigx),
          # setting_t("PML SIGMAY MAX", pml_sigy),
          # setting_t("PML SIGMAZ MAX", pml_sigz),
          # setting_t("PML INTEGRATION", pml_type),
          setting_t("TIME INTEGRATOR", time_integrator),
          setting_t("CFL NUMBER", cfl),
          setting_t("START TIME", start_time),
          setting_t("FINAL TIME", final_time),
          setting_t("OUTPUT TO FILE", output_to_file)]

def main():
  failCount=0;

  failCount += test(name="testLbsTri",
                    cmd=lbsBin,
                    settings=lbsSettings(element=3,data_file=lbsData2D,dim=2,velmodel="D2Q9"),
                    referenceNorm=0.444607812999762)

  failCount += test(name="testLbsQuad",
                    cmd=lbsBin,
                    settings=lbsSettings(element=4,data_file=lbsData2D,dim=2,velmodel="D2Q9"),
                    referenceNorm=0.446677648467864)

  failCount += test(name="testLbsTet",
                    cmd=lbsBin,
                    settings=lbsSettings(element=6,data_file=lbsData3D,dim=3, degree=2,velmodel="D3Q15"),
                    referenceNorm=0.816131769809708)

  failCount += test(name="testLbsHex",
                    cmd=lbsBin,
                    settings=lbsSettings(element=12,data_file=lbsData3D,dim=3, degree=2,velmodel="D3Q15"),
                    referenceNorm=0.818015892660822)

  failCount += test(name="testLbsTri_MPI", ranks=4,
                  cmd=lbsBin,
                  settings=lbsSettings(element=3,data_file=lbsData2D,dim=2,velmodel="D2Q9"),
                  referenceNorm=0.44497547998273)

  #clean up
  for file_name in os.listdir(testDir):
    if file_name.endswith('.vtu'):
      os.remove(testDir + "/" + file_name)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)
