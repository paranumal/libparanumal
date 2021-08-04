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

#include "lbs.hpp"

void lbs_t::latticeSetup(){

 settings.getSetting("VISCOSITY", nu);

 if (settings.compareSetting("VELOCITY MODEL", "D2Q9")) {
  // speed of sound and number of dicrete velocity fields
  c        = 1.0/sqrt(3.0);
  Nfields  = 9; 
  Nmacro = mesh.dim + 1;
    
  LBM   = (dfloat *) calloc(Nfields*Nmacro,sizeof(dfloat)); 
  LMAP  = (int *) calloc(Nfields,sizeof(int)); 

  LBM[0 + 0*Nfields] = 4.0/9.0; 
  LBM[0 + 1*Nfields] = 0.0; 
  LBM[0 + 2*Nfields] = 0.0; 

  LMAP[0]= 0;  

  // Construct lattice model
  for(int n=1; n<Nfields; n++){
    int even = (n%2)==0  ? 1 :0; 
    // weights 
    LBM[n + 0*Nfields] = even==0 ? 1.0/9.0 : 1.0/36.0; 
    // x-velocity
    LBM[n + 1*Nfields] = even==0 ? cos(M_PI*(n-1)/4.0) : sqrt(2.0)*cos(M_PI*(n-1)/4.0); 
    // y-velocity
    LBM[n + 2*Nfields] = even==0 ? sin(M_PI*(n-1)/4.0) : sqrt(2.0)*sin(M_PI*(n-1)/4.0); 

     LMAP[n] = n==4  ? (Nfields-1) : (n+4)%(Nfields-1); 

    // printf("%d \n",  LMAP[n]);
  } 

  }else if(settings.compareSetting("VELOCITY MODEL", "D3Q15")){




  }

  alpha    = nu/(c*c); // Relaxation parameter 
  tauInv   = 1.0/alpha; // need to check that....
  RT       = c*c; // remove this later, no need!!!!!
}
