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

#include "occa.hpp"
#include "mesh2D.h"

void occaTestQuad2D(mesh2D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy){

  // build OCCA device
  occa::device device;

  //---[ Device setup with string flags ]-------------------
  device.setup("mode = OpenCL  , platformID = 0, deviceID = 1");

  // OKL: OCCA Kernel Language
  occa::kernel meshGradientQuad2DKernel = 
    device.buildKernelFromSource("src/meshGradientQuad2D.okl",
				 "meshGradientQuad2D");

  // allocate DEVICE arrays
  occa::memory o_q  = 
    device.malloc(mesh->Np*mesh->Nelements*sizeof(float), q);
  occa::memory o_dqdx  = 
    device.malloc(mesh->Np*mesh->Nelements*sizeof(float), dqdx);
  occa::memory o_dqdy  = 
    device.malloc(mesh->Np*mesh->Nelements*sizeof(float), dqdy);

  occa::memory o_vgeo  = 
    device.malloc(mesh->Nvgeo*mesh->Nelements*sizeof(float),
		  mesh->vgeo);

  occa::memory o_Dr  = 
    device.malloc(mesh->Np*mesh->Np*sizeof(float),
		  mesh->Dr);

  occa::memory o_Ds  = 
    device.malloc(mesh->Np*mesh->Np*sizeof(float),
		  mesh->Ds);

  // launch kernel on DEVICE
  meshGradientQuad2DKernel(mesh->Nelements,
			   mesh->Np, 
			   mesh->Nvgeo,
			   o_vgeo,
			   o_Dr,
			   o_Ds,
			   o_q,
			   o_dqdx,
			   o_dqdy);		       

  // copy from DEVICE to HOST array
  o_dqdx.copyTo(dqdx);
  o_dqdy.copyTo(dqdy);

}

