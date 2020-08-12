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
from testGradient import *

def main():
  failCount=0;

  testMeshTriP_refnorm =[4.37016024448821, 4.4426446089857, 4.44319809556755, 4.44288849919259, 4.44288258145188, 4.44288293226037, 4.44288293848078, 4.44288293816329, 4.44288293815814, 4.44288293815836, 4.44288293815835, 4.44288293815838, 4.44288293815836, 4.44288293815835, 4.44288293815839]
  testMeshQuadP_refnorm=[4.37016024448821, 4.4426446089857, 4.4428824894834, 4.44288293763069, 4.44288293815795, 4.44288293815837, 4.44288293815837, 4.44288293815837, 4.44288293815837, 4.44288293815837, 4.44288293815837, 4.44288293815837, 4.44288293815837, 4.44288293815837, 4.44288293815837]
  testMeshTetP_refnorm =[11.9681767290924, 12.1666361728703, 12.1707573113819, 12.1674394705786, 12.1673282745486, 12.1673357411262, 12.1673360502229, 12.1673360287753]
  testMeshHexP_refnorm =[11.9681767290924, 12.1666833366009, 12.1673347991738, 12.1673360264757, 12.1673360279197, 12.1673360279208, 12.1673360279208, 12.1673360279208]

  for p in range(1, 11):
    failCount += test(name="testMeshTri_P" + str(p),
                      cmd=gradientBin,
                      settings=gradientSettings(element=3,data_file=gradientData2D,
                                                dim=2,degree=p),
                      referenceNorm=testMeshTriP_refnorm[p-1])

  for p in range(1, 11):
    failCount += test(name="testMeshQuad_P" + str(p),
                      cmd=gradientBin,
                      settings=gradientSettings(element=4,data_file=gradientData2D,
                                                dim=2,degree=p),
                      referenceNorm=testMeshQuadP_refnorm[p-1])

  for p in range(1, 9):
    failCount += test(name="testMeshTet_P" + str(p),
                      cmd=gradientBin,
                      settings=gradientSettings(element=6,data_file=gradientData3D,
                                                dim=3,degree=p),
                      referenceNorm=testMeshTetP_refnorm[p-1])

  for p in range(1, 9):
    failCount += test(name="testMeshHex_P" + str(p),
                      cmd=gradientBin,
                      settings=gradientSettings(element=12,data_file=gradientData3D,
                                                dim=3,degree=p),
                      referenceNorm=testMeshHexP_refnorm[p-1])

  failCount += test(name="testMeshTri_ReadMsh",
                    cmd=gradientBin,
                    settings=gradientSettings(element=3,data_file=gradientData2D,dim=2,
                                              mesh=testDir+"/squareTri.msh"),
                    referenceNorm=0.580787485719841)

  failCount += test(name="testMeshQuad_ReadMsh",
                    cmd=gradientBin,
                    settings=gradientSettings(element=4,data_file=gradientData2D,dim=2,
                                              mesh=testDir+"/squareQuad.msh"),
                    referenceNorm=0.580787485654967)

  failCount += test(name="testMeshTet_ReadMsh",
                    cmd=gradientBin,
                    settings=gradientSettings(element=6,data_file=gradientData3D,dim=3,
                                              mesh=testDir+"/cubeTet.msh"),
                    referenceNorm=0.942816947760423)

  failCount += test(name="testMeshHex_ReadMsh",
                    cmd=gradientBin,
                    settings=gradientSettings(element=12,data_file=gradientData3D,dim=3,
                                              mesh=testDir+"/cubeHex.msh"),
                    referenceNorm=0.942816869518335)

  failCount += test(name="testMeshTri_ReadMsh_MPI", ranks=2,
                    cmd=gradientBin,
                    settings=gradientSettings(element=3,data_file=gradientData2D,dim=2,
                                              mesh=testDir+"/squareTri.msh"),
                    referenceNorm=0.580787485719841)

  failCount += test(name="testMeshQuad_ReadMsh_MPI", ranks=2,
                    cmd=gradientBin,
                    settings=gradientSettings(element=4,data_file=gradientData2D,dim=2,
                                              mesh=testDir+"/squareQuad.msh"),
                    referenceNorm=0.580787485654967)

  failCount += test(name="testMeshTet_ReadMsh_MPI", ranks=2,
                    cmd=gradientBin,
                    settings=gradientSettings(element=6,data_file=gradientData3D,dim=3,
                                              mesh=testDir+"/cubeTet.msh"),
                    referenceNorm=0.942816947760423)

  failCount += test(name="testMeshHex_ReadMsh_MPI", ranks=2,
                    cmd=gradientBin,
                    settings=gradientSettings(element=12,data_file=gradientData3D,dim=3,
                                              mesh=testDir+"/cubeHex.msh"),
                    referenceNorm=0.942816869518335)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)