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

// // Level-Set function
// #define lssInitialConditions2D(t, x, y, q) \
// {                                       \
//   const dfloat phi_hump = 0.25 + 0.25*cos(M_PI * sqrt((x-0.25)*(x-0.25) + (y-0.5)*(y-0.5) ) / 0.15); 			    \
//   const dfloat phi_cone = 1.00 - sqrt((x-0.50)*(x-0.50) + (y-0.25)*(y-0.25) ) / 0.15f; 				\
//   const int test1 = sqrt((x-0.25)*(x-0.25) + (y-0.50)*(y-0.50)) <=0.15 ? 1:0; \
//   const int test2 = sqrt((x-0.50)*(x-0.50) + (y-0.25)*(y-0.25)) <=0.15 ? 1:0; \
//   const int test3 =(sqrt((x-0.50)*(x-0.50) + (y-0.75)*(y-0.75)) <=0.15  && (fabs(x-0.5)>=0.025 || y>=0.85)) ? 1:0; \
//   dfloat phin = 0.0; \
//   if(test1){ phin = phi_hump;}\
//   else if(test2){ phin = phi_cone;}\
//   else if(test3){phin = 1.0;}\
//   else {phin = 0.0;}\
//   *(q) = phin;\
// }


// Level-Set function for 2 humps 2 cones
#define lssInitialConditions2D(t, x, y, q) \
{                                       \
  const dfloat phi_hump1 = 0.25 + 0.25*cos(M_PI * sqrt((x-0.25)*(x-0.25) + (y-0.5)*(y-0.5) ) / 0.15);          \
  const dfloat phi_cone1 = 1.00 - sqrt((x-0.50)*(x-0.50) + (y-0.25)*(y-0.25) ) / 0.15f;        \
  const dfloat phi_hump2 = 0.25 + 0.25*cos(M_PI * sqrt((x-0.75)*(x-0.75) + (y-0.5)*(y-0.5) ) / 0.15);          \
  const dfloat phi_cone2 = 1.00 - sqrt((x-0.50)*(x-0.50) + (y-0.75)*(y-0.75) ) / 0.15f;        \
  const int test1 = sqrt((x-0.25)*(x-0.25) + (y-0.50)*(y-0.50)) <=0.15 ? 1:0; \
  const int test2 = sqrt((x-0.50)*(x-0.50) + (y-0.25)*(y-0.25)) <=0.15 ? 1:0; \
  const int test3 = sqrt((x-0.75)*(x-0.75) + (y-0.50)*(y-0.50)) <=0.15 ? 1:0; \
  const int test4 = sqrt((x-0.50)*(x-0.50) + (y-0.75)*(y-0.75)) <=0.15 ? 1:0; \
  dfloat phin = 0.0; \
  if(test1){ phin = phi_hump1;}\
  else if(test2){ phin = phi_cone1;}\
  else if(test3){phin = phi_hump2;}\
  else if(test4){phin = phi_cone2;}\
  else {phin = 0.0;}\
  *(q) = phin;\
}
// LS Advective field
#define lssAdvectionField2D(t, x, y, q, u, v) \
{                                       \
 *(u) = 2.f*M_PI*(0.5f-y); \
 *(v) = 2.f*M_PI*(x - 0.5f); \
}



// Boundary conditions
/* wall 1, outflow 2 */
#define lssDirichletConditions2D(bc, t, x, y, nx, ny, qM, qB) \
{                                       \
  if(bc==1){                            \
    *(qB) = 0.f;                        \
  } else if(bc==2){                     \
    *(qB) = 0.f;                         \
  }                                     \
}