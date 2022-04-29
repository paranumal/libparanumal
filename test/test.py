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

import os
import sys
import subprocess
import re
from pathlib import Path

#can only run the tests from LIBP_DIR/test
pwd = Path(os.getcwd())
testDir = str(pwd)
libPDir = str(pwd.parent)
solverDir = libPDir + "/solvers"

gradientDir      = solverDir + "/gradient"
advectionDir     = solverDir + "/advection"
acousticsDir     = solverDir + "/acoustics"
ellipticDir      = solverDir + "/elliptic"
fokkerPlanckDir  = solverDir + "/fokkerPlanck"
cnsDir           = solverDir + "/cns"
bnsDir           = solverDir + "/bns"
lbsDir           = solverDir + "/lbs"
insDir           = solverDir + "/ins"

gradientBin  = gradientDir      + "/gradientMain"
advectionBin = advectionDir     + "/advectionMain"
acousticsBin = acousticsDir     + "/acousticsMain"
ellipticBin  = ellipticDir      + "/ellipticMain"
fpeBin       = fokkerPlanckDir  + "/fpeMain"
cnsBin       = cnsDir           + "/cnsMain"
bnsBin       = bnsDir           + "/bnsMain"
lbsBin       = lbsDir           + "/lbsMain"
insBin       = insDir           + "/insMain"

inputRC = testDir + "/setup.rc"

TOL = 1.0e-5
alignWidth = 40

numeric_const_pattern = r"[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?"

if len(sys.argv)>1:
  device = sys.argv[1];
  if device!="Serial" and \
     device!="OpenMP" and \
     device!="CUDA"   and \
     device!="HIP"    and \
     device!="OpenCL":
    exit("Invalid mode requested.")
else:
  device="Serial"

class bcolors:
  TEST    = '\033[35m'
  PASS    = '\033[92m'
  WARNING = '\033[93m'
  FAIL    = '\033[91m'
  ENDC    = '\033[0m'

class setting_t:
  def __init__(self, name, value):
    self.name = name
    self.value = value

def writeSetup(filename, settings):
  str_settings=""
  for setting in settings:
    str_settings += "[" + setting.name + "]\n"
    str_settings += str(setting.value) + "\n\n"

  file = open(filename+".rc", "w")
  file.write(str_settings)
  file.close()

def test(name, cmd, settings, referenceNorm, ranks=1):

  #create input file
  writeSetup("setup",settings)

  #print test name
  print(bcolors.TEST + f"{name:.<{alignWidth}}" + bcolors.ENDC, end="", flush=True)

  #run test
  run = subprocess.run(["mpirun", "--oversubscribe", "-np", str(ranks), cmd, inputRC],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  if len(run.stdout.decode().splitlines())==0:
    #this failure is bad, dump the whole output for debug
    print(bcolors.FAIL + "FAIL" + bcolors.ENDC)
    print(bcolors.WARNING + name + " stdout:" + bcolors.ENDC)
    print(run.stdout.decode())
    print(bcolors.WARNING + name + " stderr:" + bcolors.ENDC)
    print(run.stderr.decode())
    #save the setup for reproducibility
    writeSetup(name,settings)
    failed = 1
  else:
    #collect last line of output
    output = run.stdout.decode().splitlines()[-1]

    #check last line's syntax
    failed=0;
    if "Solution norm = " in output:
      norm = float(output.split()[3])
      if abs(norm - referenceNorm) < TOL:
        print(bcolors.PASS + "PASS" + bcolors.ENDC)
      else:
        #failed residual check
        print(bcolors.FAIL + "FAIL" + bcolors.ENDC)
        print(bcolors.WARNING + "Expected Result: " + str(referenceNorm) + bcolors.ENDC)
        print(bcolors.WARNING + "Observed Result: " + str(norm) + bcolors.ENDC)
        #save the setup for reproducibility
        writeSetup(name,settings)
        failed = 1
    else:
      #this failure is worse, so dump the whole output for debug
      print(bcolors.FAIL + "FAIL" + bcolors.ENDC)
      print(bcolors.WARNING + name + " stdout:" + bcolors.ENDC)
      print(run.stdout.decode())
      print(bcolors.WARNING + name + " stderr:" + bcolors.ENDC)
      print(run.stderr.decode())
      #save the setup for reproducibility
      writeSetup(name,settings)
      failed = 1

  # writeSetup(name,settings)
  # print(bcolors.WARNING + name + " stdout:" + bcolors.ENDC)
  # print(run.stdout.decode())
  # print(bcolors.WARNING + name + " stderr:" + bcolors.ENDC)
  # print(run.stderr.decode())

  #clean up
  os.remove(inputRC)

  return failed

if __name__ == "__main__":
  import testMesh
  import testGradient
  import testAdvection
  import testAcoustics
  import testElliptic
  import testFokkerPlanck
  import testCns
  import testBns
  import testLbs
  import testIns
  import testTimeStepper
  import testLinearSolver
  import testParAlmond
  import testParAdogs
  import testInitialGuess

  failCount=0;
  failCount+=testMesh.main()
  failCount+=testParAdogs.main()
  failCount+=testGradient.main()
  failCount+=testAdvection.main()
  failCount+=testAcoustics.main()
  failCount+=testElliptic.main()
  failCount+=testFokkerPlanck.main()
  failCount+=testCns.main()
  failCount+=testBns.main()
  failCount+=testLbs.main()
  failCount+=testIns.main()
  failCount+=testInitialGuess.main()
  failCount+=testTimeStepper.main()
  failCount+=testLinearSolver.main()
  failCount+=testParAlmond.main()

  sys.exit(failCount)
