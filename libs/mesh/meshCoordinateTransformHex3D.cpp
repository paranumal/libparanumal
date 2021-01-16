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
#include "mesh/mesh3D.hpp"

void meshHex3D::CoordinateTransform(int _cubN, const char *_cubatureType){

  cubatureType = strdup(_cubatureType);
  
  /* */
  string mapFileName;
  settings.getSetting("BOX COORDINATE MAP FILE", mapFileName);

  if(1)
  if(mapFileName != "NONE"){
    
    dfloat epsy = 1., epsz = 1.;
    settings.getSetting("BOX COORDINATE MAP PARAMETER Y", epsy);
    settings.getSetting("BOX COORDINATE MAP PARAMETER Z", epsz);
    
    occa::properties kernelInfo = props;

    // build kernel
    occa::kernel coordMapKernel = platform.buildKernel(mapFileName, "coordMapKernel", kernelInfo);

    occa::memory o_tmpx, o_tmpy, o_tmpz, o_tmpEX, o_tmpEY, o_tmpEZ;
    o_tmpx = platform.device.malloc(Np*Nelements*sizeof(dfloat), x);
    o_tmpy = platform.device.malloc(Np*Nelements*sizeof(dfloat), y);
    o_tmpz = platform.device.malloc(Np*Nelements*sizeof(dfloat), z);
    
    coordMapKernel(Np*Nelements, epsy, epsz, o_tmpx, o_tmpy, o_tmpz);

    o_tmpEX = platform.device.malloc(Nverts*Nelements*sizeof(dfloat), EX);
    o_tmpEY = platform.device.malloc(Nverts*Nelements*sizeof(dfloat), EY);
    o_tmpEZ = platform.device.malloc(Nverts*Nelements*sizeof(dfloat), EZ);

    coordMapKernel(Nverts*Nelements, epsy, epsz, o_tmpEX, o_tmpEY, o_tmpEZ);

    o_tmpx.copyTo(x);
    o_tmpy.copyTo(y);
    o_tmpz.copyTo(z);

    o_tmpEX.copyTo(EX);
    o_tmpEY.copyTo(EY);
    o_tmpEZ.copyTo(EZ);
  }
  
  halo->Exchange(x, Np, ogs_dfloat);
  halo->Exchange(y, Np, ogs_dfloat);
  halo->Exchange(z, Np, ogs_dfloat);

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
