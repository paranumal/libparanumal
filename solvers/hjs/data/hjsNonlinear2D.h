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
A 2D nonlinear equation: 
phi_t + phi_x*phi_y = 0 ; phi_0 = sin(x) + cons(y)
in [-pi, pi] x [-pi, pi] with periodic boundaries
Solution is smooth at t = 0.8
Discontinous derivatives at t = 1.5
*/

// Level-Set function
#define hjsInitialConditions2D(t, x, y, q) \
{                                       \
  (*q) = sin(x) + cos(y);\
}

#define hjsComputeHamiltonian2D(t,x,y,p1,p2,q1,q2,ham, dhdp, dhdq){\
  const dfloat p = 0.5*(p1+p2);\
  const dfloat q = 0.5*(q1+q2);\
  (*dhdp) = 1.0f; \
  (*dhdq) = 1.0f; \
  (*ham)  = p*q; \
} 
 
// #define hjsComputeHamiltonian2D(t,x,y,p1,p2,q1,q2,ham, dhdp, dhdq){\
//   const dfloat p = 0.5*(p1+p2);\
//   const dfloat q = 0.5*(q1+q2);\
//   (*dhdp) = fabs(q); \
//   (*dhdq) = fabs(p); \
//   (*ham)  = p*q; \
// } 