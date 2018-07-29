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


//Ricker pulse
#define ricker(t, f) \
 (1-2*OCCA_PI*OCCA_PI*f*f*t*t)*occaExp(-OCCA_PI*OCCA_PI*f*f*t*t) 

//integrated Ricker pulse
#define intRicker(t, f) \
  t*occaExp(-OCCA_PI*OCCA_PI*f*f*t*t) 

#define acousticsPointSource2D(x, y, t, f, c, u, v, p)  \
{ \
  dfloat r = sqrt(x*x+y*y); \
                              \
  p = ricker((t-r/c),f)/(4*OCCA_PI*c*c*r); \
  u = x*(intRicker((t-r/c),f)/r + ricker((t-r/c),f)/c)/(4*OCCA_PI*c*c*r*r); \
  v = y*(intRicker((t-r/c),f)/r + ricker((t-r/c),f)/c)/(4*OCCA_PI*c*c*r*r); \
}