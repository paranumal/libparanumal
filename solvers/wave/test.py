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
waveBin = "./waveMain"
inputRC = "./setup.rc"

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

def test(name, cmd, settings, ranks=1):

  #create input file
  writeSetup("setup",settings)

  #print test name
  print(bcolors.TEST + f"{name:.<{alignWidth}}" + bcolors.ENDC, end="", flush=True)

  #run test
  run = subprocess.run(["mpirun", "--oversubscribe", "-np", str(ranks), cmd, inputRC],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  writeSetup(name,settings)

  fp = open(name+".out", "w")
  fp.write(run.stdout.decode())
  fp.close()
  
  # writeSetup(name,settings)
  # print(bcolors.WARNING + name + " stdout:" + bcolors.ENDC)
  # print(run.stdout.decode())
  # print(bcolors.WARNING + name + " stderr:" + bcolors.ENDC)
  # print(run.stderr.decode())

  #clean up
  os.remove(inputRC)

  return 

if __name__ == "__main__":
  import testWave

  failCount=0;
  failCount+=testWave.main()
  sys.exit(failCount)
