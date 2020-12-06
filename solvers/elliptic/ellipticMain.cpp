/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "elliptic.hpp"
#include "mesh/meshDefines3D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  if(argc!=2)
    LIBP_ABORT(string("Usage: ./ellipticMain setupfile"));

  //create default settings
  platformSettings_t platformSettings(comm);
  meshSettings_t meshSettings(comm);
  ellipticSettings_t ellipticSettings(comm);
  ellipticAddRunSettings(ellipticSettings);

  //load settings from file
  ellipticSettings.parseFromFile(platformSettings, meshSettings,
                                 argv[1]);

  // set up platform
  platform_t platform(platformSettings);

  platformSettings.report();
  meshSettings.report();
  ellipticSettings.report();

  // set up mesh
  mesh_t& mesh = mesh_t::Setup(platform, meshSettings, comm);

  dfloat lambda = 0.0;
  ellipticSettings.getSetting("LAMBDA", lambda);

  // Boundary Type translation. Just defaults.
  int NBCTypes = 3;
  int BCType[NBCTypes] = {0,1,2};

  // set up elliptic solver
  elliptic_t& elliptic = elliptic_t::Setup(platform, mesh, ellipticSettings,
                                           lambda, NBCTypes, BCType);

  // run
  elliptic.Run();

#if 1
  printf("Building Matrix\n");
  
  parAlmond::parCOO A(elliptic.platform, mesh.comm);
  elliptic.BuildOperatorMatrixContinuous(A);

  occa::properties kernelInfo = mesh.props;

  if(sizeof(hlong)==8)
    kernelInfo["defines/hlong"]= "long long int";
  if(sizeof(hlong)==4)
    kernelInfo["defines/hlong"]= "int";

  kernelInfo["defines/" "p_G00ID"]= G00ID;
  kernelInfo["defines/" "p_G01ID"]= G01ID;
  kernelInfo["defines/" "p_G02ID"]= G02ID;
  kernelInfo["defines/" "p_G11ID"]= G11ID;
  kernelInfo["defines/" "p_G12ID"]= G12ID;
  kernelInfo["defines/" "p_G22ID"]= G22ID;
  kernelInfo["defines/" "p_GWJID"]= GWJID;

  char kernelName[BUFSIZ];
  switch(mesh.elementType){
  case TRIANGLES:
    sprintf(kernelName, "ellipticBuildOperatorMatrixContinuousTri2D");
    break;
  case QUADRILATERALS:
    sprintf(kernelName, "ellipticBuildOperatorMatrixContinuousQuad2D");
    break;
  case TETRAHEDRA:
    sprintf(kernelName, "ellipticBuildOperatorMatrixContinuousTet3D");
    break;
  case HEXAHEDRA:
    sprintf(kernelName, "ellipticBuildOperatorMatrixContinuousHex3D");
    break;
  }

  occa::kernel buildMatrixKernel =
    elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticBuildOperatorMatrixContinuous.okl",
				  kernelName,
				  kernelInfo);
  
  occa::memory o_maskedGlobalNumbering =
    platform.malloc(mesh.Np*mesh.Nelements*sizeof(hlong), elliptic.maskedGlobalNumbering);

  occa::memory o_Srr, o_Srs, o_Srt, o_Sss, o_Sst, o_Stt;
  o_Srr = platform.malloc(mesh.Np*mesh.Np*sizeof(dfloat), mesh.Srr);
  o_Srs = platform.malloc(mesh.Np*mesh.Np*sizeof(dfloat), mesh.Srs);
  o_Sss = platform.malloc(mesh.Np*mesh.Np*sizeof(dfloat), mesh.Sss);
  if(mesh.dim==3 && mesh.elementType==TETRAHEDRA){
    o_Srt = platform.malloc(mesh.Np*mesh.Np*sizeof(dfloat), mesh.Srt);
    o_Sst = platform.malloc(mesh.Np*mesh.Np*sizeof(dfloat), mesh.Sst);
    o_Stt = platform.malloc(mesh.Np*mesh.Np*sizeof(dfloat), mesh.Stt);
  }
  
  occa::memory o_AL =
    platform.malloc(mesh.Nelements*mesh.Np*mesh.Np*sizeof(parAlmond::parCOO::nonZero_t));

  switch(mesh.elementType){
  case TRIANGLES:
    buildMatrixKernel(mesh.Nelements, o_maskedGlobalNumbering,
		      o_Srr, o_Srs, o_Sss,
		      mesh.o_MM, mesh.o_ggeo, elliptic.lambda, o_AL);
    break;
  case QUADRILATERALS:
    buildMatrixKernel(mesh.Nelements, o_maskedGlobalNumbering,
		      mesh.o_D,
		      mesh.o_ggeo, elliptic.lambda, o_AL);
    break;
  case TETRAHEDRA:
    buildMatrixKernel(mesh.Nelements, o_maskedGlobalNumbering,
		      o_Srr, o_Srs, o_Srt, o_Sss, o_Sst, o_Stt,
		      mesh.o_MM, mesh.o_ggeo, elliptic.lambda, o_AL);
    break;
  case HEXAHEDRA:
    buildMatrixKernel(mesh.Nelements, o_maskedGlobalNumbering,
		      mesh.o_D,
		      mesh.o_ggeo, elliptic.lambda, o_AL);
    break;
  }


  
#endif
  
  // close down MPI
  MPI_Finalize();
  return LIBP_SUCCESS;
}
