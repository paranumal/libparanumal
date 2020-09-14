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

ellipticData2D = ellipticDir + "/data/ellipticSine2D.h"
ellipticData3D = ellipticDir + "/data/ellipticSine3D.h"

def ellipticSettings(rcformat="2.0", data_file=ellipticData2D,
                     mesh="BOX", dim=2, element=4, nx=10, ny=10, nz=10, boundary_flag=1,
                     degree=4, thread_model=device, platform_number=0, device_number=0,
                     Lambda=1.0,
                     discretization="CONTINUOUS",
                     linear_solver="PCG",
                     precon="MULTIGRID",
                     multigrid_smoother="CHEBYSHEV",
                     paralmond_cycle="VCYCLE",
                     paralmond_strength="SYMMETRIC",
                     paralmond_aggregation="UNSMOOTHED",
                     paralmond_smoother="CHEBYSHEV",
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
          setting_t("DISCRETIZATION", discretization),
          setting_t("LINEAR SOLVER", linear_solver),
          setting_t("PRECONDITIONER", precon),
          setting_t("MULTIGRID SMOOTHER", multigrid_smoother),
          setting_t("PARALMOND CYCLE", paralmond_cycle),
          setting_t("PARALMOND STRENGTH", paralmond_strength),
          setting_t("PARALMOND AGGREGATION", paralmond_aggregation),
          setting_t("PARALMOND SMOOTHER", paralmond_smoother),
          setting_t("OUTPUT TO FILE", "FALSE"),
          setting_t("VERBOSE", output_to_file)]

def main():
  failCount=0;

  ################
  #C0 tests
  ################
  failCount += test(name="testEllipticTri_C0",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2, precon="NONE"),
                    referenceNorm=0.500000001211135)

  failCount += test(name="testEllipticQuad_C0",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=4,data_file=ellipticData2D,dim=2, precon="NONE"),
                    referenceNorm=0.499999999969716)

  failCount += test(name="testEllipticTet_C0",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3, precon="NONE"),
                    referenceNorm=0.353553400508458)

  failCount += test(name="testEllipticHex_C0",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=12,data_file=ellipticData3D,dim=3, precon="NONE"),
                    referenceNorm=0.353553390458384)

  failCount += test(name="testEllipticQuad3D_C0",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=4,data_file=ellipticData3D,mesh="sphereQuad.msh", dim=3, precon="NONE"),
                    referenceNorm=3.274235742251)


  #C0 precons
  #tri
  failCount += test(name="testEllipticTri_C0_Jacobi",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="JACOBI"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticTri_C0_Massmatrix",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="MASSMATRIX"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticTri_C0_ParAlmond",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="PARALMOND"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticTri_C0_Multigrid",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="MULTIGRID"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticTri_C0_Semfem",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="SEMFEM"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticTri_C0_OAS",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="OAS"),
                    referenceNorm=0.500000001211135)

  #quad
  failCount += test(name="testEllipticQuad_C0_Jacobi",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=4,data_file=ellipticData2D,dim=2,
                                              precon="JACOBI"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticQuad_C0_ParAlmond",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=4,data_file=ellipticData2D,dim=2,
                                              precon="PARALMOND"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticQuad_C0_Multigrid",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=4,data_file=ellipticData2D,dim=2,
                                              precon="MULTIGRID"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticQuad_C0_Semfem",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=4,data_file=ellipticData2D,dim=2,
                                              precon="SEMFEM"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticQuad_C0_OAS",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=4,data_file=ellipticData2D,dim=2,
                                              precon="OAS"),
                    referenceNorm=0.500000001211135)

  #tet
  failCount += test(name="testEllipticTet_C0_Jacobi",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3,
                                              precon="JACOBI"),
                    referenceNorm=0.353553400508458)
  failCount += test(name="testEllipticTet_C0_Massmatrix",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3,
                                              precon="MASSMATRIX"),
                    referenceNorm=0.353553400508458)
  failCount += test(name="testEllipticTet_C0_ParAlmond",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3, degree=2,
                                              precon="PARALMOND"),
                    referenceNorm=0.353474740220582)
  failCount += test(name="testEllipticTet_C0_Multigrid",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3,
                                              precon="MULTIGRID"),
                    referenceNorm=0.353553400508458)
  failCount += test(name="testEllipticTet_C0_Semfem",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3,
                                              precon="SEMFEM"),
                    referenceNorm=0.353553400508458)
  failCount += test(name="testEllipticTet_C0_OAS",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3,
                                              precon="OAS"),
                    referenceNorm=0.353553400508458)

  #hex
  failCount += test(name="testEllipticHex_C0_Jacobi",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=12,data_file=ellipticData3D,dim=3,
                                              precon="JACOBI"),
                    referenceNorm=0.353553400508458)
  failCount += test(name="testEllipticHex_C0_ParAlmond",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=12,data_file=ellipticData3D,dim=3, degree=2,
                                              precon="PARALMOND"),
                    referenceNorm=0.353553400508458)
  failCount += test(name="testEllipticHex_C0_Multigrid",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=12,data_file=ellipticData3D,dim=3,
                                              precon="MULTIGRID"),
                    referenceNorm=0.353553400508458)
  failCount += test(name="testEllipticHex_C0_Semfem",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=12,data_file=ellipticData3D,dim=3,
                                              precon="SEMFEM"),
                    referenceNorm=0.353553400508458)
  failCount += test(name="testEllipticHex_C0_OAS",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=12,data_file=ellipticData3D,dim=3,
                                              precon="OAS"),
                    referenceNorm=0.353553400508458)

  # all Neumann
  failCount += test(name="testEllipticTri_C0_AllNeumann",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              boundary_flag=-1, Lambda=0.0),
                    referenceNorm=0.0962635430608342)

  failCount += test(name="testEllipticQuad_C0_AllNeumann",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=4,data_file=ellipticData2D,dim=2,
                                              boundary_flag=-1, Lambda=0.0),
                    referenceNorm=0.0962635338676865)

  failCount += test(name="testEllipticTet_C0_AllNeumann",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3,
                                              boundary_flag=-1, Lambda=0.0),
                    referenceNorm=0.0595408371412352)

  failCount += test(name="testEllipticHex_C0_AllNeumann",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=12,data_file=ellipticData3D,dim=3,
                                              boundary_flag=-1, Lambda=0.0),
                    referenceNorm=0.059540839002614)

  ################
  #IPDG tests
  ################
  failCount += test(name="testEllipticTri_Ipdg",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="NONE", discretization="IPDG"),
                    referenceNorm=0.500000001211135)

  failCount += test(name="testEllipticQuad_Ipdg",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=4,data_file=ellipticData2D,dim=2,
                                              precon="NONE", discretization="IPDG"),
                    referenceNorm=0.499999999969716)

  failCount += test(name="testEllipticTet_Ipdg",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3,
                                              precon="NONE", discretization="IPDG"),
                    referenceNorm=0.353553400508458)

  failCount += test(name="testEllipticHex_Ipdg",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=12,data_file=ellipticData3D,dim=3,
                                              precon="NONE", discretization="IPDG"),
                    referenceNorm=0.353553390458384)

  #IPDG precons
  #tri
  failCount += test(name="testEllipticTri_Ipdg_Jacobi",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="JACOBI", discretization="IPDG"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticTri_Ipdg_Massmatrix",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="MASSMATRIX", discretization="IPDG"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticTri_Ipdg_ParAlmond",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="PARALMOND", discretization="IPDG"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticTri_Ipdg_Multigrid",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="MULTIGRID", discretization="IPDG"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticTri_Ipdg_OAS",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="OAS", discretization="IPDG"),
                    referenceNorm=0.500000001211135)

  #quad
  failCount += test(name="testEllipticQuad_Ipdg_Jacobi",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=4,data_file=ellipticData2D,dim=2,
                                              precon="JACOBI", discretization="IPDG"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticQuad_Ipdg_ParAlmond",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=4,data_file=ellipticData2D,dim=2,
                                              precon="PARALMOND", discretization="IPDG"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticQuad_Ipdg_Multigrid",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=4,data_file=ellipticData2D,dim=2,
                                              precon="MULTIGRID", discretization="IPDG"),
                    referenceNorm=0.500000001211135)
  failCount += test(name="testEllipticQuad_Ipdg_OAS",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=4,data_file=ellipticData2D,dim=2,
                                              precon="OAS", discretization="IPDG"),
                    referenceNorm=0.500000001211135)

  #tet
  failCount += test(name="testEllipticTet_Ipdg_Jacobi",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3,
                                              precon="JACOBI", discretization="IPDG"),
                    referenceNorm=0.353553400508458)
  failCount += test(name="testEllipticTet_Ipdg_Massmatrix",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3,
                                              precon="MASSMATRIX", discretization="IPDG"),
                    referenceNorm=0.353553400508458)
  failCount += test(name="testEllipticTet_Ipdg_ParAlmond",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3, degree=2,
                                              precon="PARALMOND", discretization="IPDG"),
                    referenceNorm=0.353502126562155)
  failCount += test(name="testEllipticTet_Ipdg_Multigrid",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3,
                                              precon="MULTIGRID", discretization="IPDG"),
                    referenceNorm=0.353553400508458)
  failCount += test(name="testEllipticTet_Ipdg_OAS",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3,
                                              precon="OAS", discretization="IPDG"),
                    referenceNorm=0.353553400508458)

  #hex
  failCount += test(name="testEllipticHex_Ipdg_Jacobi",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=12,data_file=ellipticData3D,dim=3,
                                              precon="JACOBI", discretization="IPDG"),
                    referenceNorm=0.353553400508458)
  failCount += test(name="testEllipticHex_Ipdg_ParAlmond",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=12,data_file=ellipticData3D,dim=3, degree=2,
                                              precon="PARALMOND", paralmond_aggregation="UNSMOOTHED",
                                              discretization="IPDG"),
                    referenceNorm=0.353553400508458)
  failCount += test(name="testEllipticHex_Ipdg_Multigrid",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=12,data_file=ellipticData3D,dim=3,
                                              precon="MULTIGRID", discretization="IPDG"),
                    referenceNorm=0.353553400508458)
  failCount += test(name="testEllipticHex_Ipdg_OAS",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=12,data_file=ellipticData3D,dim=3,
                                              precon="OAS", discretization="IPDG"),
                    referenceNorm=0.353553400508458)

  # all Neumann
  failCount += test(name="testEllipticTri_Ipdg_AllNeumann",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              boundary_flag=-1, Lambda=0.0,
                                              discretization="IPDG"),
                    referenceNorm=0.0962635430608342)

  failCount += test(name="testEllipticQuad_Ipdg_AllNeumann",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=4,data_file=ellipticData2D,dim=2,
                                              boundary_flag=-1, Lambda=0.0,
                                              discretization="IPDG"),
                    referenceNorm=0.0962635338676865)

  failCount += test(name="testEllipticTet_Ipdg_AllNeumann",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=6,data_file=ellipticData3D,dim=3,
                                              boundary_flag=-1, Lambda=0.0,
                                              discretization="IPDG"),
                    referenceNorm=0.0595408371412352)

  failCount += test(name="testEllipticHex_Ipdg_AllNeumann",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=12,data_file=ellipticData3D,dim=3,
                                              boundary_flag=-1, Lambda=0.0,
                                              discretization="IPDG"),
                    referenceNorm=0.059540839002614)

  #large pMG test
  failCount += test(name="testEllipticTri_C0_Multigrid_large",
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              nx=100, ny=100,
                                              precon="MULTIGRID"),
                                              referenceNorm=0.500000001211135)
  #MPI tests
  failCount += test(name="testEllipticTri_C0_Multigrid_MPI", ranks=4,
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="MULTIGRID"),
                    referenceNorm=0.500000001211135)

  failCount += test(name="testEllipticTri_Ipdg_Multigrid_MPI", ranks=4,
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="MULTIGRID", discretization="IPDG", output_to_file="TRUE"),
                    referenceNorm=0.500000001211135)

  failCount += test(name="testEllipticTri_C0_OAS_MPI", ranks=4,
                    cmd=ellipticBin,
                    settings=ellipticSettings(element=3,data_file=ellipticData2D,dim=2,
                                              precon="OAS"),
                    referenceNorm=0.500000001211135)

  #clean up
  for file_name in os.listdir(testDir):
    if file_name.endswith('.vtu'):
      os.remove(testDir + "/" + file_name)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)
