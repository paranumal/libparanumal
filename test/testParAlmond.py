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
from testElliptic import *

def main():
  failCount=0;

  failCount += test(name="testParAlmond_Vcycle_jacobi",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,
                                              dim=2, precon="PARALMOND",
                                              paralmond_cycle="VCYCLE",
                                              paralmond_smoother="DAMPEDJACOBI"),
                    referenceNorm=0.500000001211135)

  failCount += test(name="testParAlmond_Kcycle_jacobi",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,
                                              dim=2, precon="PARALMOND",
                                              paralmond_cycle="KCYCLE",
                                              paralmond_smoother="DAMPEDJACOBI"),
                    referenceNorm=0.500000001211135)

  failCount += test(name="testParAlmond_Vcycle_cheby",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,
                                              dim=2, precon="PARALMOND",
                                              paralmond_cycle="VCYCLE",
                                              paralmond_smoother="CHEBYSHEV"),
                    referenceNorm=0.500000001211135)

  failCount += test(name="testParAlmond_Kcycle_cheby",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,
                                              dim=2, precon="PARALMOND",
                                              paralmond_cycle="KCYCLE",
                                              paralmond_smoother="CHEBYSHEV"),
                    referenceNorm=0.500000001211135)

  failCount += test(name="testParAlmond_Vcycle_cheby_MPI", ranks=4,
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,
                                              dim=2, precon="PARALMOND",
                                              paralmond_cycle="VCYCLE",
                                              paralmond_smoother="CHEBYSHEV"),
                    referenceNorm=0.500000001211135)

  failCount += test(name="testParAlmond_Kcycle_cheby_MPI", ranks=4,
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,
                                              dim=2, precon="PARALMOND",
                                              paralmond_cycle="KCYCLE",
                                              paralmond_smoother="CHEBYSHEV"),
                    referenceNorm=0.500000001211135)

  # Ruge-Stuben strength
  failCount += test(name="testParAlmond_Vcycle_rugestuben_MPI", ranks=4,
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,
                                              dim=2, precon="PARALMOND",
                                              paralmond_cycle="VCYCLE",
                                              paralmond_strength="RUGESTUBEN",
                                              paralmond_smoother="CHEBYSHEV"),
                    referenceNorm=0.500000000429642)

  # unsmoothed aggregation
  failCount += test(name="testParAlmond_Vcycle_smoothed_MPI", ranks=4,
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,
                                              dim=2, precon="PARALMOND",
                                              paralmond_cycle="KCYCLE",
                                              paralmond_aggregation="SMOOTHED",
                                              paralmond_smoother="CHEBYSHEV"),
                    referenceNorm=0.500000001211135)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)