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
from testGradient import *

def main():
  failCount=0;

  failCount += test(name="testParAdogsTri_Inertial_MPI", ranks=2,
                    cmd=gradientBin,
                    settings=gradientSettings(element=3,data_file=gradientData2D,dim=2,
                                              mesh=testDir+"/squareTri.msh",
                                              paradogs_partitioning="INERTIAL"),
                    referenceNorm=0.580787485719841)

  failCount += test(name="testParAdogsQuad_Inertial_MPI", ranks=2,
                    cmd=gradientBin,
                    settings=gradientSettings(element=4,data_file=gradientData2D,dim=2,
                                              mesh=testDir+"/squareQuad.msh",
                                              paradogs_partitioning="INERTIAL"),
                    referenceNorm=0.580787485654967)

  failCount += test(name="testParAdogsTet_Inertial_MPI", ranks=2,
                    cmd=gradientBin,
                    settings=gradientSettings(element=6,data_file=gradientData3D,dim=3,
                                              mesh=testDir+"/cubeTet.msh",
                                              paradogs_partitioning="INERTIAL"),
                    referenceNorm=0.942816947760423)

  failCount += test(name="testParAdogsHex_Inertial_MPI", ranks=2,
                    cmd=gradientBin,
                    settings=gradientSettings(element=12,data_file=gradientData3D,dim=3,
                                              mesh=testDir+"/cubeHex.msh",
                                              paradogs_partitioning="INERTIAL"),
                    referenceNorm=0.942816869518335)

  failCount += test(name="testParAdogsTri_Spectral_MPI", ranks=2,
                    cmd=gradientBin,
                    settings=gradientSettings(element=3,data_file=gradientData2D,dim=2,
                                              mesh=testDir+"/squareTri.msh",
                                              paradogs_partitioning="SPECTRAL"),
                    referenceNorm=0.580787485719841)

  failCount += test(name="testParAdogsQuad_Spectral_MPI", ranks=2,
                    cmd=gradientBin,
                    settings=gradientSettings(element=4,data_file=gradientData2D,dim=2,
                                              mesh=testDir+"/squareQuad.msh",
                                              paradogs_partitioning="SPECTRAL"),
                    referenceNorm=0.580787485654967)

  failCount += test(name="testParAdogsTet_Spectral_MPI", ranks=2,
                    cmd=gradientBin,
                    settings=gradientSettings(element=6,data_file=gradientData3D,dim=3,
                                              mesh=testDir+"/cubeTet.msh",
                                              paradogs_partitioning="SPECTRAL"),
                    referenceNorm=0.942816947760423)

  failCount += test(name="testParAdogsHex_Spectral_MPI", ranks=2,
                    cmd=gradientBin,
                    settings=gradientSettings(element=12,data_file=gradientData3D,dim=3,
                                              mesh=testDir+"/cubeHex.msh",
                                              paradogs_partitioning="SPECTRAL"),
                    referenceNorm=0.942816869518335)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)
