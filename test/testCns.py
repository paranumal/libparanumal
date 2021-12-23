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

cnsData2D = cnsDir + "/data/cnsGaussian2D.h"
cnsData3D = cnsDir + "/data/cnsGaussian3D.h"

def cnsSettings(rcformat="2.0", data_file=cnsData2D,
               mesh="BOX", dim=2, element=4, nx=10, ny=10, nz=10, boundary_flag=-1,
               degree=4, thread_model=device, platform_number=0, device_number=0,
               gamma=1.4, viscosity=0.01, isothermal="FALSE",
               advection_type="COLLOCATION",
                time_integrator="DOPRI5", cfl=1.0, start_time=0.0, final_time=1.0,
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
          setting_t("GAMMA", gamma),
          setting_t("VISCOSITY", viscosity),
          setting_t("ISOTHERMAL", isothermal),
          setting_t("ADVECTION TYPE", advection_type),
          setting_t("TIME INTEGRATOR", time_integrator),
          setting_t("CFL NUMBER", cfl),
          setting_t("START TIME", start_time),
          setting_t("FINAL TIME", final_time),
          setting_t("OUTPUT TO FILE", output_to_file)]

def main():
  failCount=0;

  failCount += test(name="testCnsTri",
                    cmd=cnsBin,
                    settings=cnsSettings(element=3,data_file=cnsData2D,dim=2),
                    referenceNorm=27.4583308011963)

  failCount += test(name="testCnsQuad",
                    cmd=cnsBin,
                    settings=cnsSettings(element=4,data_file=cnsData2D,dim=2),
                    referenceNorm=27.4595638641015)

  failCount += test(name="testCnsTet",
                    cmd=cnsBin,
                    settings=cnsSettings(element=6,data_file=cnsData3D,dim=3, degree=2),
                    referenceNorm=85.2379630054725)

  failCount += test(name="testCnsHex",
                    cmd=cnsBin,
                    settings=cnsSettings(element=12,data_file=cnsData3D,dim=3, degree=2),
                    referenceNorm=85.2656652300985)

  failCount += test(name="testCnsTri_cub",
                    cmd=cnsBin,
                    settings=cnsSettings(element=3,data_file=cnsData2D,dim=2,
                                         advection_type="CUBATURE"),
                    referenceNorm=27.4584573000865)

  failCount += test(name="testCnsQuad_cub",
                    cmd=cnsBin,
                    settings=cnsSettings(element=4,data_file=cnsData2D,dim=2,
                                         advection_type="CUBATURE"),
                    referenceNorm=27.4586844667892)

  failCount += test(name="testCnsTet_cub",
                    cmd=cnsBin,
                    settings=cnsSettings(element=6,data_file=cnsData3D,dim=3, degree=2,
                                         advection_type="CUBATURE"),
                    referenceNorm=85.2549800928041)

  failCount += test(name="testCnsHex_cub",
                    cmd=cnsBin,
                    settings=cnsSettings(element=12,data_file=cnsData3D,dim=3, degree=2,
                                         advection_type="CUBATURE"),
                    referenceNorm=85.2652515937174)

  failCount += test(name="testCnsTri_Isothermal",
                    cmd=cnsBin,
                    settings=cnsSettings(element=3,data_file=cnsData2D,dim=2,
                                         isothermal="TRUE"),
                    referenceNorm=10.2200732624953)

  failCount += test(name="testCnsQuad_Isothermal",
                    cmd=cnsBin,
                    settings=cnsSettings(element=4,data_file=cnsData2D,dim=2,
                                         isothermal="TRUE"),
                    referenceNorm=10.2201442207598)

  failCount += test(name="testCnsTet_Isothermal",
                    cmd=cnsBin,
                    settings=cnsSettings(element=6,data_file=cnsData3D,dim=3,
                                         isothermal="TRUE",
                                         nx=6, ny=6, nz=6, degree=3),
                    referenceNorm=31.6748226159438)

  failCount += test(name="testCnsHex_Isothermal",
                    cmd=cnsBin,
                    settings=cnsSettings(element=12,data_file=cnsData3D,dim=3,
                                         isothermal="TRUE",
                                         nx=8, ny=8, nz=8, degree=2),
                    referenceNorm=31.6617290760014)

  failCount += test(name="testCnsTri_Isothermal_cub",
                    cmd=cnsBin,
                    settings=cnsSettings(element=3,data_file=cnsData2D,dim=2,
                                         isothermal="TRUE",
                                         advection_type="CUBATURE"),
                    referenceNorm=10.2197628476831)

  failCount += test(name="testCnsQuad_Isothermal_cub",
                    cmd=cnsBin,
                    settings=cnsSettings(element=4,data_file=cnsData2D,dim=2,
                                         isothermal="TRUE",
                                         advection_type="CUBATURE"),
                    referenceNorm=10.2197250517919)

  failCount += test(name="testCnsTet_Isothermal_cub",
                    cmd=cnsBin,
                    settings=cnsSettings(element=6,data_file=cnsData3D,dim=3,
                                         isothermal="TRUE",
                                         advection_type="CUBATURE",
                                         nx=6, ny=6, nz=6, degree=2),
                    referenceNorm=31.6390441400446)

  failCount += test(name="testCnsHex_Isothermal_cub",
                    cmd=cnsBin,
                    settings=cnsSettings(element=12,data_file=cnsData3D,dim=3,
                                         isothermal="TRUE",
                                         advection_type="CUBATURE",
                                         nx=8, ny=8, nz=8, degree=2),
                    referenceNorm=31.6605632403763)

  failCount += test(name="testCnsTri_MPI", ranks=4,
                    cmd=cnsBin,
                    settings=cnsSettings(element=3,data_file=cnsData2D,dim=2, output_to_file="TRUE"),
                    referenceNorm=27.4600335839337)

  #clean up
  for file_name in os.listdir(testDir):
    if file_name.endswith('.vtu'):
      os.remove(testDir + "/" + file_name)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)
