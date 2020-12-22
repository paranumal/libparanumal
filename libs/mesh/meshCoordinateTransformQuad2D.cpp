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

#include "mesh.hpp"
#include "mesh/mesh2D.hpp"

void meshQuad2D::CoordinateTransform(int _cubN, const char *cubatureType){

  /* */
  string mapFileName;
  settings.getSetting("BOX COORDINATE MAP FILE", mapFileName);

  if(mapFileName != "NONE"){
    printf("MAPPING COORDINATES\n");
    
    dfloat epsy = .3;
    occa::properties kernelInfo = props;

    // build kernel
    occa::kernel coordMapKernel = platform.buildKernel(mapFileName, "coordMapKernel", kernelInfo);

    occa::memory o_tmpx, o_tmpy;
    o_tmpx = platform.device.malloc(Np*Nelements*sizeof(dfloat), x);
    o_tmpy = platform.device.malloc(Np*Nelements*sizeof(dfloat), y);
    
    coordMapKernel(Np*Nelements, epsy, o_tmpx, o_tmpy);
    
    o_tmpx.copyTo(x);
    o_tmpy.copyTo(y);
  }
  
  halo->Exchange(x, Np, ogs_dfloat);
  halo->Exchange(y, Np, ogs_dfloat);

  // compute geometric factors
  GeometricFactors();

  // compute surface geofacs
  SurfaceGeometricFactors();

  // compute cubature stuff
  CubatureSetup(_cubN, cubatureType);
  
  // copy to DEVICE
  o_vgeo = platform.malloc((Nelements+totalHaloPairs)*Nvgeo*Np*sizeof(dfloat), vgeo);
  o_sgeo = platform.malloc(Nelements*Nfaces*Nfp*Nsgeo*sizeof(dfloat), sgeo);
  o_ggeo = platform.malloc(Nelements*Np*Nggeo*sizeof(dfloat), ggeo);

}
