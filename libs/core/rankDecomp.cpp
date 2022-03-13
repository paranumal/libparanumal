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

#include "core.hpp"

// find a factorization n = nx*ny such that
//  nx>=ny are 'close' to one another
void factor2(const int n, int &nx, int &ny) {
  //start with guessing nx ~= n^1/2
  nx = round(sqrt(n));
  ny = 1;

  for (;nx<n;nx++) {
    if (n % nx ==0) { //if nx divides n
      ny = n / nx; //divide out nx

      //swap if needed
      if (ny>nx) std::swap(nx,ny);

      return;
    }
  }

  //if we made it this far, n is prime
  nx = n;
}

// find a factorization n = nx*ny*nz such that
//  nx>=ny>=nz are all 'close' to one another
void factor3(const int n, int &nx, int &ny, int &nz) {
  //start with guessing nx ~= n^1/3
  nx = round(std::cbrt(n));
  ny = nz = 1;

  for (;nx<n;nx++) {
    if (n % nx ==0) { //if nx divides n
      const int f = n / nx; //divide out nx

      ny = round(sqrt(f)); //guess ny ~= sqrt(f)
      for (;ny<f;ny++) {
        if (f % ny == 0) { //if ny divides f
          nz = f/ny; //divide out ny

          //sort
          if (ny>nx) std::swap(nx,ny);
          if (nz>ny) std::swap(ny,nz);
          if (ny>nx) std::swap(nx,ny);

          return;
        }
      }

      //if we're here, f is prime
      ny = f;
      nz = 1;

      //swap if needed
      if (ny>nx) std::swap(nx,ny);

      return;
    }
  }

  //if we made it this far, n is prime
  nx = n;
}