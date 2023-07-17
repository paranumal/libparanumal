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

dataTime2D = "./data/waveQuadraticTime2D.h"
dataTime3D = "./data/waveQuadraticTime3D.h"

dataSpace2D = "./data/waveQuadraticSpace2D.h"
dataSpace3D = "./data/waveQuadraticSpace3D.h"

dataSpaceTime2D = "./data/waveSpaceTime2D.h"
dataSpaceTime3D = "./data/waveSpaceTime3D.h"


def waveSettings(rcformat="2.0",
                 data_file=dataTime2D,
                 mesh="BOX",
                 dim=2,
                 element=4,
                 nx=10,
                 ny=10,
                 nz=10,
                 boundary_flag=-1,
                 degree=4,
                 thread_model=device,
                 platform_number=0,
                 device_number=0,
#                 time_integrator="ESDIRK4(3)6L{2}SA",
                 time_integrator="ESDIRK6(5)9L{2}SA",
                 time_step=0.01,
                 start_time=0.0,
                 final_time=1.1,
                 discretization="IPDG",
                 stopping_criteria="ERRORESTIMATE",
                 convergence_tolerance=1e-11,
                 pmg_chebyshev_degree=2,
                 initial_guess="QR",
                 initial_guess_dim=16,
                 output_to_file="FALSE",
                 output_error_interval=0):
  
  return [setting_t("FORMAT", rcformat),
          setting_t("DATA FILE", data_file),
          setting_t("SOLVER MODE", "TIMEDOMAIN"),
          setting_t("MESH FILE", mesh),
          setting_t("MESH DIMENSION", dim),
          setting_t("ELEMENT TYPE", element),
          setting_t("BOX NX", nx),
          setting_t("BOX NY", ny),
          setting_t("BOX NZ", nz),
          setting_t("BOX DIMX", 2.),
          setting_t("BOX DIMY", 2.),
          setting_t("BOX DIMZ", 2.),
          setting_t("BOX BOUNDARY FLAG", boundary_flag),
          setting_t("POLYNOMIAL DEGREE", degree),
          setting_t("THREAD MODEL", thread_model),
          setting_t("PLATFORM NUMBER", platform_number),
          setting_t("DEVICE NUMBER", device_number),
          setting_t("TIME INTEGRATOR", time_integrator),
          setting_t("TIME STEP", time_step),
          setting_t("START TIME", start_time),
          setting_t("FINAL TIME", final_time),
          setting_t("ELLIPTIC DISCRETIZATION", discretization),
          setting_t("ELLIPTIC LINEAR SOLVER", "PCG"),
          setting_t("ELLIPTIC PRECONDITIONER", "MULTIGRID"),
          setting_t("ELLIPTIC STOPPING CRITERIA", stopping_criteria),
          setting_t("ELLIPTIC ITERATIVE CONVERGENCE TOLERANCE", convergence_tolerance),
          setting_t("ELLIPTIC MULTIGRID COARSENING", "HALFDOFS"),
          setting_t("ELLIPTIC MULTIGRID SMOOTHER", "CHEBYSHEV"),
          setting_t("ELLIPTIC MULTIGRID CHEBYSHEV DEGREE", pmg_chebyshev_degree),
          setting_t("ELLIPTIC PARALMOND STRENGTH", "RUGESTUBEN"),
          setting_t("ELLIPTIC INITIAL GUESS STRATEGY", initial_guess),
          setting_t("ELLIPTIC INITIAL GUESS HISTORY SPACE DIMENSION", initial_guess_dim),
          setting_t("ELLIPTIC PARALMOND CYCLE", "VCYCLE"),
          setting_t("ELLIPTIC PARALMOND AGGREGATION", "SMOOTHED"),
          setting_t("ELLIPTIC PARALMOND SMOOTHER", "CHEBYSHEV"),
          setting_t("ELLIPTIC PARALMOND CHEBYSHEV DEGREE", 2),
          setting_t("OUTPUT TO FILE", output_to_file),
          setting_t("OUTPUT ERROR INTERVAL", output_error_interval)]

def main():
  failCount=0;

  cnt = 1;

  testDim = 3
  testElement = 12
  testDataFile = dataSpaceTime3D
  
  for P in range(1,10):
    for NXP in range(0, 4):
      for dtp in range(0, 6):

        dt = 0.4/(2.**dtp)
        NX = 4*(2**NXP)
        test(name="testSpaceTime"+str(cnt).zfill(5),
             cmd=waveBin,
             settings=waveSettings(element=testElement,
                                   data_file=testDataFile,
                                   thread_model="CUDA",
                                   boundary_flag=1,
                                   dim=testDim,
                                   degree=P,
                                   time_step=dt,
                                   nx=NX,
                                   ny=NX,
                                   output_error_interval=1))
        print("\n")
        cnt = cnt+1;
        

  if 0==1:        
    cnt = 1;
    for dtp in range(0, 4):
      for P in range(4,5):
        for NXP in range(3, 4):
          NX = 4*(2**NXP)
          dt = 0.1/(2.**dtp)
          test(name="testSpace"+str(cnt).zfill(5),
               cmd=waveBin,
               settings=waveSettings(element=3,
                                     data_file=dataSpace2D,
                                     thread_model="CUDA",
                                     boundary_flag=1,
                                     dim=2,
                                     degree=P,
                                     time_step=dt,
                                     nx=NX,
                                     ny=NX,
                                     output_error_interval=1))
          print("\n")
          cnt = cnt+1;


  
          cnt = 1;
          for P in range(1,10):
            for NXP in range(0, 4):
              NX = 4*(2**NXP)
              test(name="testTime"+str(cnt).zfill(5),
                   cmd=waveBin,
                   settings=waveSettings(element=3,
                                         data_file=dataTime2D,
                                         thread_model="CUDA",
                                         boundary_flag=1,
                                         dim=2,
                                         degree=P,
                                         nx=NX,
                                         ny=NX,
                                         output_error_interval=1))
              print("\n")
              cnt = cnt+1;

      
  testDir = "./"
  #clean up
  for file_name in os.listdir(testDir):
    if file_name.endswith('.vtu'):
      os.remove(testDir + "/" + file_name)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)
