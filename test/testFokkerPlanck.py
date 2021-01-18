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

fpeData2D = fokkerPlanckDir + "/data/fpeLinear2D.h"
fpeData3D = fokkerPlanckDir + "/data/fpeLinear3D.h"

def fpeSettings(rcformat="2.0", data_file=fpeData2D,
                     mesh="BOX", dim=2, element=4, nx=10, ny=10, nz=10, boundary_flag=-1,
                     degree=4, thread_model=device, platform_number=0, device_number=0,
                     viscosity=0.01,
                     advection_type="COLLOCATION",
                     time_integrator="EXTBDF3", cfl=1.0, start_time=0.0, final_time=1.0,
                     num_subcycles=16, subcycle_integrator="DOPRI5",
                     elliptic_discretization="IPDG",
                     elliptic_linear_solver="PCG",
                     elliptic_precon="MULTIGRID",
                     elliptic_multigrid_smoother="CHEBYSHEV",
                     elliptic_paralmond_cycle="VCYCLE",
                     elliptic_paralmond_smoother="CHEBYSHEV",
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
          setting_t("ADVECTION TYPE", advection_type),
          setting_t("TIME INTEGRATOR", time_integrator),
          setting_t("CFL NUMBER", cfl),
          setting_t("NUMBER OF SUBCYCLES", num_subcycles),
          setting_t("SUBCYCLING TIME INTEGRATOR", subcycle_integrator),
          setting_t("START TIME", start_time),
          setting_t("FINAL TIME", final_time),
          setting_t("ELLIPTIC DISCRETIZATION", elliptic_discretization),
          setting_t("ELLIPTIC LINEAR SOLVER", elliptic_linear_solver),
          setting_t("ELLIPTIC PRECONDITIONER", elliptic_precon),
          setting_t("ELLIPTIC MULTIGRID SMOOTHER", elliptic_multigrid_smoother),
          setting_t("ELLIPTIC PARALMOND CYCLE", elliptic_paralmond_cycle),
          setting_t("ELLIPTIC PARALMOND SMOOTHER", elliptic_paralmond_smoother),
          setting_t("ELLIPTIC VERBOSE", "TRUE"),
          setting_t("OUTPUT TO FILE", output_to_file)]

def main():
  failCount=0;

  #test non-subcycle
  failCount += test(name="testFpeTri",
                    cmd=fpeBin,
                    settings=fpeSettings(element=3,data_file=fpeData2D,dim=2),
                    referenceNorm=0.684376309866456)

  failCount += test(name="testFpeQuad",
                    cmd=fpeBin,
                    settings=fpeSettings(element=4,data_file=fpeData2D,dim=2),
                    referenceNorm=0.684243684323532)

  failCount += test(name="testFpeTet",
                    cmd=fpeBin,
                    settings=fpeSettings(element=6,data_file=fpeData3D,dim=3,
                                         nx=6, ny=6, nz=6, degree=2),
                    referenceNorm=0.492149400541327)

  failCount += test(name="testFpeHex",
                    cmd=fpeBin,
                    settings=fpeSettings(element=12,data_file=fpeData3D,dim=3,
                                         nx=6, ny=6, nz=6, degree=2),
                    referenceNorm=0.455956622483511)

  #cubature tests
  failCount += test(name="testFpeTri_cub",
                    cmd=fpeBin,
                    settings=fpeSettings(element=3,data_file=fpeData2D,dim=2,
                                         advection_type="CUBATURE"),
                    referenceNorm=0.684376309866456)

  failCount += test(name="testFpeQuad_cub",
                    cmd=fpeBin,
                    settings=fpeSettings(element=4,data_file=fpeData2D,dim=2,
                                         advection_type="CUBATURE"),
                    referenceNorm=0.684243684323532)

  failCount += test(name="testFpeTet_cub",
                    cmd=fpeBin,
                    settings=fpeSettings(element=6,data_file=fpeData3D,dim=3,
                                         nx=6, ny=6, nz=6, degree=2,
                                         advection_type="CUBATURE"),
                    referenceNorm=0.492149400541327)

  failCount += test(name="testFpeHex_cub",
                    cmd=fpeBin,
                    settings=fpeSettings(element=12,data_file=fpeData3D,dim=3,
                                         nx=6, ny=6, nz=6, degree=2,
                                         advection_type="CUBATURE"),
                    referenceNorm=0.459185797449144)

  #test subcycle
  failCount += test(name="testFpeTri_ss",
                    cmd=fpeBin,
                    settings=fpeSettings(element=3,data_file=fpeData2D,dim=2,
                                         time_integrator="SSBDF3"),
                    referenceNorm=0.67676248716463)

  failCount += test(name="testFpeQuad_ss",
                    cmd=fpeBin,
                    settings=fpeSettings(element=4,data_file=fpeData2D,dim=2,
                                         time_integrator="SSBDF3"),
                    referenceNorm=0.678629090148834)

  failCount += test(name="testFpeTet_ss",
                    cmd=fpeBin,
                    settings=fpeSettings(element=6,data_file=fpeData3D,dim=3,
                                         nx=6, ny=6, nz=6, degree=2,
                                         time_integrator="SSBDF3"),
                    referenceNorm=0.464060993161937)

  failCount += test(name="testFpeHex_ss",
                    cmd=fpeBin,
                    settings=fpeSettings(element=12,data_file=fpeData3D,dim=3,
                                         nx=6, ny=6, nz=6, degree=2,
                                         time_integrator="SSBDF3"),
                    referenceNorm=0.389099535288541)

  #cubature tests
  failCount += test(name="testFpeTri_ss_cub",
                    cmd=fpeBin,
                    settings=fpeSettings(element=3,data_file=fpeData2D,dim=2,
                                         time_integrator="SSBDF3",
                                         advection_type="CUBATURE"),
                    referenceNorm=0.676762487164629)

  failCount += test(name="testFpeQuad_ss_cub",
                    cmd=fpeBin,
                    settings=fpeSettings(element=4,data_file=fpeData2D,dim=2,
                                         time_integrator="SSBDF3",
                                         advection_type="CUBATURE"),
                    referenceNorm=0.678621302711163)

  failCount += test(name="testFpeTet_ss_cub",
                    cmd=fpeBin,
                    settings=fpeSettings(element=6,data_file=fpeData3D,dim=3,
                                         nx=6, ny=6, nz=6, degree=2,
                                         time_integrator="SSBDF3",
                                         advection_type="CUBATURE"),
                    referenceNorm=0.464060993161937)

  failCount += test(name="testFpeHex_ss_cub",
                    cmd=fpeBin,
                    settings=fpeSettings(element=12,data_file=fpeData3D,dim=3,
                                         nx=6, ny=6, nz=6, degree=2,
                                         time_integrator="SSBDF3",
                                         advection_type="CUBATURE"),
                    referenceNorm=0.390359551549932)

  #test MPI
  failCount += test(name="testFpeTri_MPI", ranks=4,
                    cmd=fpeBin,
                    settings=fpeSettings(element=3,data_file=fpeData2D,dim=2,output_to_file="TRUE"),
                    referenceNorm=0.683519576398795)

  #clean up
  for file_name in os.listdir(testDir):
    if file_name.endswith('.vtu'):
      os.remove(testDir + "/" + file_name)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)
