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

// 
Burgers equation in 3D : 
phi_t + 0.5*(phi_x^2 +phi_y^2 + +phi_z^2 + 1) = 0
in [-2, 2] x [-2, 2]x [-2, 2] with periodic boundaries
Solution is smooth at t = 0.5/pi^2 
Discontinous derivatives at t = 1.5/pi^2

*/


// Level-Set function
#define hjsInitialConditions3D(t, x, y, z, q) \
{                                       \
  (*q) = -cos(0.5f*M_PI*(x+y+z));\
}

 
#define hjsComputeHamiltonian3D(t,x,y,z,p1,p2,q1,q2,m1,m2,ham, dhdp, dhdq, dhdm){\
  const dfloat p = 0.5*(p1+p2);\
  const dfloat q = 0.5*(q1+q2);\
  const dfloat m = 0.5*(m1+m2);\
  (*dhdp) = 1.0f; \
  (*dhdq) = 1.0f; \
  (*dhdm) = 1.0f; \
  (*ham)  = 0.5f*(p+q+m+1.f)*(p+q+m+1.f); \
} 

// (*dhdp) = 2.f*fabs(p + q + 1.f); 
// (*dhdq) = 2.f*fabs(p + q + 1.f); 