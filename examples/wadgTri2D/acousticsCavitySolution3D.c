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

#include <math.h>
#include "mesh3D.h"

void acousticsCavitySolution3D(dfloat x, dfloat y, dfloat z, dfloat time,
			       dfloat *u, dfloat *v, dfloat *w, dfloat *p){

  // dudt = -dpdx
  // dvdt = -dpdy
  // dwdt = -dpdz
  // dpdt = -dudx -dvdy -dwdz

  // boundary conditions
  // dpdn = 0
  // u.n = 0

  dfloat pi = M_PI;
  
  *u = sin(pi*time*sqrt(3.)/2.)*cos(pi*x/2.)*sin(pi*y/2.)*sin(pi*z/2.);
  *v = sin(pi*time*sqrt(3.)/2.)*sin(pi*x/2.)*cos(pi*y/2.)*sin(pi*z/2.);
  *w = sin(pi*time*sqrt(3.)/2.)*sin(pi*x/2.)*sin(pi*y/2.)*cos(pi*z/2.);
  *p = -sqrt(3.)*cos(pi*time*sqrt(3.)/2.)*sin(pi*x/2.)*sin(pi*y/2.)*sin(pi*z/2.);

}

