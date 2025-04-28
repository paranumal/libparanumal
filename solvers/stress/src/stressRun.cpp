/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "stress.hpp"
#include "timer.hpp"

void stress_t::Run(){

  //setup linear algebra module
  platform.linAlg().InitKernels({"set"});

  //setup linear solver
  hlong NglobalDofs;
  NglobalDofs = ogsMasked.NgatherGlobal*Nfields;

  linearSolver_t<dfloat> linearSolver;
  linearSolver.Setup<LinearSolver::pcg<dfloat> >(Ndofs, Nhalo, platform, settings, comm);

  properties_t kernelInfo = mesh.props; //copy base occa properties

  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  //add standard boundary functions
  std::string boundaryHeaderFileName;
  if (mesh.dim==2)
    boundaryHeaderFileName = std::string(DSTRESS "/data/stressBoundary2D.h");
  else if (mesh.dim==3)
    boundaryHeaderFileName = std::string(DSTRESS "/data/stressBoundary3D.h");
  kernelInfo["includes"] += boundaryHeaderFileName;

  int Nmax = std::max(mesh.Np, mesh.Nfaces*mesh.Nfp);
  kernelInfo["defines/" "p_Nmax"]= Nmax;

  kernelInfo["defines/" "p_Nfields"]= Nfields;

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();
  
  std::string oklFilePrefix = DSTRESS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  fileName   = oklFilePrefix + "stressRhs" + suffix + oklFileSuffix;
  kernelName = "stressRhs" + suffix;
  kernel_t forcingKernel = platform.buildKernel(fileName, kernelName,
                                                    kernelInfo);

  kernel_t rhsBCKernel, addBCKernel;
  fileName   = oklFilePrefix + "stressRhsBC" + suffix + oklFileSuffix;
  kernelName = "stressRhsBC" + suffix;
  
  rhsBCKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
  
  fileName   = oklFilePrefix + "stressAddBC" + suffix + oklFileSuffix;
  kernelName = "stressAddBC" + suffix;
  
  addBCKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  //create occa buffers
  dlong Nall = Nfields*mesh.Np*(mesh.Nelements+mesh.totalHaloPairs);
  memory<dfloat> rL(Nall);
  memory<dfloat> xL(Nall);
  deviceMemory<dfloat> o_rL = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_xL = platform.malloc<dfloat>(Nall);

  deviceMemory<dfloat> o_r, o_x;
  dlong Ng = ogsMasked.Ngather;
  dlong Nghalo = gHalo.Nhalo;
  dlong Ngall = Nfields*(Ng + Nghalo);
  o_r = platform.malloc<dfloat>(Ngall);
  o_x = platform.malloc<dfloat>(Ngall);

  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  //populate rhs forcing
  forcingKernel(mesh.Nelements,
                mesh.o_wJ,
                mesh.o_MM,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                lambda,
                o_rL);

  //Set x to zero
  platform.linAlg().set(mesh.Nelements*mesh.Np*Nfields, (dfloat)0.0, o_xL);

  rhsBCKernel(mesh.Nelements,
	      mesh.o_wJ,
	      mesh.o_vgeo,
	      mesh.o_sgeo,
	      mesh.o_D,
	      mesh.o_S,
	      mesh.o_MM,
	      mesh.o_vmapM,
	      mesh.o_sM,
	      lambda,
	      o_nut,
	      mesh.o_x,
	      mesh.o_y,
	      mesh.o_z,
	      o_mapB,
	      o_rL);

  // gather rhs to globalDofs if c0
  ogsMasked.Gather(o_r, o_rL, Nfields, ogs::Add, ogs::Trans);
  ogsMasked.Gather(o_x, o_xL, Nfields, ogs::Add, ogs::NoTrans);

  int maxIter = 5000;
  int verbose = settings.compareSetting("VERBOSE", "TRUE") ? 1 : 0;

  timePoint_t start = GlobalPlatformTime(platform);

  //call the solver
  dfloat tol = (sizeof(dfloat)==sizeof(double)) ? 1.0e-8 : 1.0e-5;
  int iter = Solve(linearSolver, o_x, o_r, tol, maxIter, verbose);

  //add the boundary data to the masked nodes
  // scatter x to LocalDofs if c0
  ogsMasked.Scatter(o_xL, o_x, Nfields, ogs::NoTrans);

  //fill masked nodes with BC data
  addBCKernel(mesh.Nelements,
	      mesh.o_x,
	      mesh.o_y,
	      mesh.o_z,
	      o_mapB,
	      o_xL);

  timePoint_t end = GlobalPlatformTime(platform);
  double elapsedTime = ElapsedTime(start, end);

  if ((mesh.rank==0) && verbose){
    printf("%d, " hlongFormat ", %g, %d, %g, %g; global: N, dofs, elapsed, iterations, time per node, nodes*iterations/time %s\n",
           mesh.N,
           NglobalDofs,
           elapsedTime,
           iter,
           elapsedTime/(NglobalDofs),
           NglobalDofs*((dfloat)iter/elapsedTime),
           (char*) settings.getSetting("PRECONDITIONER").c_str());
  }

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

    // copy data back to host
    o_xL.copyTo(xL);

    // output field files
    std::string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d.vtu", name.c_str(), mesh.rank);

    PlotFields(xL, fname);
  }

  // output norm of final solution
  {
    //compute q.M*q
    dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
    deviceMemory<dfloat> o_MxL = platform.reserve<dfloat>(Nentries);
    mesh.MassMatrixApply(o_xL, o_MxL);

    dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_xL, o_MxL, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);
  }
}
