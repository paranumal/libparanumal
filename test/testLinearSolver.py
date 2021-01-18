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

  failCount += test(name="testLinearSolver_PCG",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="NONE", linear_solver="PCG"),
                    referenceNorm=0.500000001211135)

  failCount += test(name="testLinearSolver_FPCG",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="NONE", linear_solver="FPCG"),
                    referenceNorm=0.500000001211135)

  failCount += test(name="testLinearSolver_NBPCG",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="NONE", linear_solver="NBPCG"),
                    referenceNorm=0.500000001211135)

  failCount += test(name="testLinearSolver_NBFPCG",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="NONE", linear_solver="NBFPCG"),
                    referenceNorm=0.500000001211135)

  failCount += test(name="testLinearSolver_PGMRES",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="NONE", linear_solver="PGMRES"),
                    referenceNorm=0.500000001211135)

  failCount += test(name="testLinearSolver_PMINRES",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="NONE", linear_solver="PMINRES"),
                    referenceNorm=0.500000001211135)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)
