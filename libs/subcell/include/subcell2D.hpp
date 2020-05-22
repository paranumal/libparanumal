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

#ifndef SUBCELL2D_HPP
#define SUBCELL2D_HPP 1
#include "subcellDefines2D.hpp"

class subcellTri2D: public subcell_t {
public:
  // subcellTri2D(mesh_t& _mesh, settings_t& _settings):
  //  subcell_t(_mesh, _settings) {}
  subcellTri2D(mesh_t& _mesh, settings_t& _settings);

   void SetupDetector(); 
   void OccaSetup(); 

   void CreateMinorGrid(); 
   void GeometricFactors(); 
   // void LocalConnect();  
   // void GlobalConnect();  

  // void ParallelReader(const char *fileName);
  // void SetupBox();
  // void SetupPmlBox();
  // void ReferenceNodes(int N);
  // void PhysicalNodes();
  // void GeometricFactors();
  // void SurfaceGeometricFactors();
  // void OccaSetup();

  // void CubatureSetup();
  // void CubatureNodes();
};

class subcellQuad2D: public subcell_t {
public:
  subcellQuad2D(mesh_t& _mesh, settings_t& _settings):
   subcell_t(_mesh, _settings) {printf("created a new subcellQuad2D\n");}
  // void ParallelReader(const char *fileName);
  // void SetupBox();
  // void SetupPmlBox();
  // void ReferenceNodes(int N);
  // void PhysicalNodes();
  // void GeometricFactors();
  // void SurfaceGeometricFactors();
  // void OccaSetup();

  // void CubatureSetup();
  // void CubatureNodes();
};

#endif

