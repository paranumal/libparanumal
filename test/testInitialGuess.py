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

from test import *

insData2D = insDir + "/data/insVortex2D.h"
insData3D = insDir + "/data/insBeltrami3D.h"

def insSettings(rcformat="2.0", data_file=insData2D,
               mesh="BOX", dim=2, element=4, nx=10, ny=10, nz=10, boundary_flag=2,
               degree=4, thread_model=device, platform_number=0, device_number=0,
               viscosity=0.05,
               advection_type="COLLOCATION",
               time_integrator="EXTBDF3", cfl=1.0, start_time=0.0, final_time=0.1,
               num_subcycles=4, subcycle_integrator="DOPRI5",
               velocity_discretization="CONTINUOUS",
               velocity_linear_solver="PCG",
               velocity_precon="JACOBI",
               velocity_multigrid_smoother="CHEBYSHEV",
               velocity_paralmond_cycle="VCYCLE",
               velocity_paralmond_smoother="CHEBYSHEV",
               velocity_initial_guess_strategy="EXTRAP",
               velocity_initial_guess_history_space_dimension=9,
               velocity_initial_guess_extrap_degree=3,
               velocity_initial_guess_extrap_coeffs_method="MINNORM",
               pressure_discretization="CONTINUOUS",
               pressure_linear_solver="FPCG",
               pressure_precon="MULTIGRID",
               pressure_multigrid_smoother="CHEBYSHEV",
               pressure_paralmond_cycle="KCYCLE",
               pressure_paralmond_smoother="CHEBYSHEV",
               pressure_initial_guess_strategy="EXTRAP",
               pressure_initial_guess_history_space_dimension=9,
               pressure_initial_guess_extrap_degree=3,
               pressure_initial_guess_extrap_coeffs_method="MINNORM",
                output_to_file="FALSE"):
  return [setting_t("FORMAT", rcformat),
          setting_t("DATA FILE", data_file),
          setting_t("MESH FILE", mesh),
          setting_t("MESH DIMENSION", dim),
          setting_t("ELEMENT TYPE", element),
          setting_t("BOX NX", nx),
          setting_t("BOX NY", ny),
          setting_t("BOX NZ", nz),
          setting_t("BOX DIMX", 1),
          setting_t("BOX DIMY", 1),
          setting_t("BOX DIMZ", 1),
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
          setting_t("VELOCITY DISCRETIZATION", velocity_discretization),
          setting_t("VELOCITY LINEAR SOLVER", velocity_linear_solver),
          setting_t("VELOCITY PRECONDITIONER", velocity_precon),
          setting_t("VELOCITY MULTIGRID SMOOTHER", velocity_multigrid_smoother),
          setting_t("VELOCITY PARALMOND CYCLE", velocity_paralmond_cycle),
          setting_t("VELOCITY PARALMOND SMOOTHER", velocity_paralmond_smoother),
          setting_t("VELOCITY VERBOSE", "TRUE"),
          setting_t("VELOCITY INITIAL GUESS STRATEGY", velocity_initial_guess_strategy),
          setting_t("VELOCITY INITIAL GUESS HISTORY SPACE DIMENSION", velocity_initial_guess_history_space_dimension),
          setting_t("VELOCITY INITIAL GUESS EXTRAP DEGREE", velocity_initial_guess_extrap_degree),
          setting_t("VELOCITY INITIAL GUESS EXTRAP COEFFS METHOD",  velocity_initial_guess_extrap_coeffs_method),
          setting_t("PRESSURE DISCRETIZATION", pressure_discretization),
          setting_t("PRESSURE LINEAR SOLVER", pressure_linear_solver),
          setting_t("PRESSURE PRECONDITIONER", pressure_precon),
          setting_t("PRESSURE MULTIGRID SMOOTHER", pressure_multigrid_smoother),
          setting_t("PRESSURE PARALMOND CYCLE", pressure_paralmond_cycle),
          setting_t("PRESSURE PARALMOND SMOOTHER", pressure_paralmond_smoother),
          setting_t("PRESSURE VERBOSE", "TRUE"),
          setting_t("PRESSURE INITIAL GUESS STRATEGY", pressure_initial_guess_strategy),
          setting_t("PRESSURE INITIAL GUESS HISTORY SPACE DIMENSION", pressure_initial_guess_history_space_dimension),
          setting_t("PRESSURE INITIAL GUESS EXTRAP DEGREE", pressure_initial_guess_extrap_degree),
          setting_t("PRESSURE INITIAL GUESS EXTRAP COEFFS METHOD",  pressure_initial_guess_extrap_coeffs_method),
          setting_t("OUTPUT TO FILE", output_to_file)]

def main():
  failCount=0;


  # IG
  #test classic
  failCount += test(name="testInitialGuess_VextrapPqr",
                    cmd=insBin,
                    settings=insSettings(element=3,data_file=insData2D,dim=2,
                                         velocity_initial_guess_strategy="QR",
                                         pressure_initial_guess_strategy="QR"),
                    referenceNorm=0.821033993796193)

  failCount += test(name="testInitialGuess_VextrapPqr",
                    cmd=insBin,
                    settings=insSettings(element=4,data_file=insData2D,dim=2,
                                         velocity_initial_guess_strategy="EXTRAP",
                                         pressure_initial_guess_strategy="QR"),
                    referenceNorm=0.818161264240046)

  failCount += test(name="testInitialGuess_VqrPextrap",
                    cmd=insBin,
                    settings=insSettings(element=6,data_file=insData3D,dim=3,
                                         velocity_initial_guess_strategy="QR",
                                         pressure_initial_guess_strategy="EXTRAP",
                                         nx=6, ny=6, nz=6, degree=2),
                  referenceNorm=1.1985917685975)

  failCount += test(name="testInitialGuess_VspextrapPspextrap",
                    cmd=insBin,
                    settings=insSettings(element=12,data_file=insData3D,dim=3,
                                         velocity_initial_guess_strategy="EXTRAP",
                                         velocity_initial_guess_extrap_coeffs_method="CPQR",
                                         pressure_initial_guess_strategy="EXTRAP",
                                         pressure_initial_guess_extrap_coeffs_method="CPQR",
                                         nx=6, ny=6, nz=6, degree=2),
                    referenceNorm=1.1956470417617)

  failCount += test(name="testInitialGuess_VclassicPclassic",
                    cmd=insBin,
                    settings=insSettings(element=12,data_file=insData3D,dim=3,
                                         nx=6, ny=6, nz=6, degree=2,
                                         advection_type="CUBATURE",
                                         time_integrator="SSBDF3",
                                         velocity_initial_guess_strategy="CLASSIC",
                                         pressure_initial_guess_strategy="CLASSIC"),
                    referenceNorm=1.17790533313334)


  #test wth MPI
  failCount += test(name="testInitialGuess_MPI", ranks=4,
                    cmd=insBin,
                    settings=insSettings(element=3,data_file=insData2D,dim=2,output_to_file="TRUE",
                                         velocity_initial_guess_strategy="QR",
                                         pressure_initial_guess_strategy="QR"),
                    referenceNorm=0.820949431342165)

  #clean up
  for file_name in os.listdir(testDir):
    if file_name.endswith('.vtu'):
      os.remove(testDir + "/" + file_name)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)
