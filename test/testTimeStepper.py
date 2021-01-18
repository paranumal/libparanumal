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
from testAdvection import *
from testFokkerPlanck import *
from testBns import *

def main():
  failCount=0;

  failCount += test(name="testTimeStepper_ab3",
                    cmd=advectionBin,
                    settings=advectionSettings(element=3,data_file=advectionData2D,
                                               dim=2, time_integrator="AB3", cfl=0.25),
                    referenceNorm=0.723972801309193)

  failCount += test(name="testTimeStepper_dopri5",
                    cmd=advectionBin,
                    settings=advectionSettings(element=3,data_file=advectionData2D,
                                               dim=2, time_integrator="DOPRI5"),
                    referenceNorm=0.723924419144375)

  failCount += test(name="testTimeStepper_lserk4",
                    cmd=advectionBin,
                    settings=advectionSettings(element=3,data_file=advectionData2D,
                                               dim=2, time_integrator="LSERK4"),
                    referenceNorm=0.723924546941676)

  failCount += test(name="testTimeStepper_extbdf3",
                    cmd=fpeBin,
                    settings=fpeSettings(element=3,data_file=fpeData2D,dim=2,
                                         time_integrator="EXTBDF3"),
                    referenceNorm=0.684376309866456)

  failCount += test(name="testTimeStepper_ssbdf3",
                    cmd=fpeBin,
                    settings=fpeSettings(element=3,data_file=fpeData2D,dim=2,
                                         time_integrator="SSBDF3"),
                    referenceNorm=0.67676248716463)

  failCount += test(name="testTimeStepper_ab3_pml",
                    cmd=bnsBin,
                    settings=bnsSettings(element=3,data_file=bnsData2D,dim=2,
                                         time_integrator="AB3", cfl=0.25),
                    referenceNorm=14.2550696512202)

  failCount += test(name="testTimeStepper_saab3_pml",
                    cmd=bnsBin,
                    settings=bnsSettings(element=3,data_file=bnsData2D,dim=2,
                                         time_integrator="SAAB3", cfl=0.25),
                    referenceNorm=14.2550270959095)

  failCount += test(name="testTimeStepper_dopri5_pml",
                    cmd=bnsBin,
                    settings=bnsSettings(element=3,data_file=bnsData2D,dim=2,
                                         time_integrator="DOPRI5"),
                    referenceNorm=14.2550469661139)

  failCount += test(name="testTimeStepper_lserk4_pml",
                    cmd=bnsBin,
                    settings=bnsSettings(element=3,data_file=bnsData2D,dim=2,
                                         time_integrator="LSERK4"),
                    referenceNorm=14.2550469661139)

  failCount += test(name="testTimeStepper_sark4_pml",
                    cmd=bnsBin,
                    settings=bnsSettings(element=3,data_file=bnsData2D,dim=2,
                                         time_integrator="SARK4"),
                    referenceNorm=14.2114867833305)

  failCount += test(name="testTimeStepper_sark5_pml",
                    cmd=bnsBin,
                    settings=bnsSettings(element=3,data_file=bnsData2D,dim=2,
                                         time_integrator="SARK5"),
                    referenceNorm=14.2397306802799)

  failCount += test(name="testTimeStepper_mrab3_pml",
                    cmd=bnsBin,
                    settings=bnsSettings(element=3,data_file=bnsData2D,dim=2,
                                         time_integrator="MRAB3", cfl=0.25),
                    referenceNorm=14.2550696512202)

  failCount += test(name="testTimeStepper_mrsaab3_pml",
                    cmd=bnsBin,
                    settings=bnsSettings(element=3,data_file=bnsData2D,dim=2,
                                         time_integrator="MRSAAB3", cfl=0.25),
                    referenceNorm=14.2550270959095)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)