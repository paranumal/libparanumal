/*

  The MIT License (MIT)

  Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

  if (settings.compareSetting("VELOCITY MODEL", "D2Q9") && mesh.dim==2) {
    // speed of sound and number of dicrete velocity fields
    c        = 1.0/sqrt(3.0);
    Nfields  = 9; 
    Nmacro = mesh.dim + 1;
    
    LBM.malloc(Nfields*Nmacro);
    LMAP.malloc(Nfields);

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
    } 

  }else if(settings.compareSetting("VELOCITY MODEL", "D3Q15") && mesh.dim==3){

    c        = 1.0/sqrt(3.0);
    Nfields  = 15; 
    Nmacro   = mesh.dim + 1;

    LBM.malloc(Nfields*Nmacro);
    LMAP.malloc(Nfields);
    // Weights
    LBM[0  + 0*Nfields] = 2.0/ 9.0; LMAP[0 ] = 0; 

    LBM[1  + 0*Nfields] = 1.0/ 9.0; LMAP[1 ] = 2; 
    LBM[2  + 0*Nfields] = 1.0/ 9.0; LMAP[2 ] = 1; 
    LBM[3  + 0*Nfields] = 1.0/ 9.0; LMAP[3 ] = 4; 
    LBM[4  + 0*Nfields] = 1.0/ 9.0; LMAP[4 ] = 3; 
    LBM[5  + 0*Nfields] = 1.0/ 9.0; LMAP[5 ] = 6; 
    LBM[6  + 0*Nfields] = 1.0/ 9.0; LMAP[6 ] = 5; 

    LBM[7  + 0*Nfields] = 1.0/72.0; LMAP[7 ] = 8; 
    LBM[8  + 0*Nfields] = 1.0/72.0; LMAP[8 ] = 7; 
    LBM[9  + 0*Nfields] = 1.0/72.0; LMAP[9 ] = 10; 
    LBM[10 + 0*Nfields] = 1.0/72.0; LMAP[10] = 9; 
    LBM[11 + 0*Nfields] = 1.0/72.0; LMAP[11] = 12; 
    LBM[12 + 0*Nfields] = 1.0/72.0; LMAP[12] = 11; 
    LBM[13 + 0*Nfields] = 1.0/72.0; LMAP[13] = 14; 
    LBM[14 + 0*Nfields] = 1.0/72.0; LMAP[14] = 13; 

    // Velocities
    LBM[0  + 1*Nfields] =  0.0; LBM[0  + 2*Nfields] =  0.0; LBM[0  + 3*Nfields] =  0.0;

    LBM[1  + 1*Nfields] =  1.0; LBM[1  + 2*Nfields] =  0.0; LBM[1  + 3*Nfields] =  0.0;
    LBM[2  + 1*Nfields] = -1.0; LBM[2  + 2*Nfields] =  0.0; LBM[2  + 3*Nfields] =  0.0;
    LBM[3  + 1*Nfields] =  0.0; LBM[3  + 2*Nfields] =  1.0; LBM[3  + 3*Nfields] =  0.0;
    LBM[4  + 1*Nfields] =  0.0; LBM[4  + 2*Nfields] = -1.0; LBM[4  + 3*Nfields] =  0.0;
    LBM[5  + 1*Nfields] =  0.0; LBM[5  + 2*Nfields] =  0.0; LBM[5  + 3*Nfields] =  1.0;
    LBM[6  + 1*Nfields] =  0.0; LBM[6  + 2*Nfields] =  0.0; LBM[6  + 3*Nfields] = -1.0;

    LBM[7   + 1*Nfields] =  1.0; LBM[7   + 2*Nfields] =  1.0; LBM[7   + 3*Nfields] =  1.0;
    LBM[8   + 1*Nfields] = -1.0; LBM[8   + 2*Nfields] = -1.0; LBM[8   + 3*Nfields] = -1.0;
    LBM[9   + 1*Nfields] = -1.0; LBM[9   + 2*Nfields] =  1.0; LBM[9   + 3*Nfields] =  1.0;
    LBM[10  + 1*Nfields] =  1.0; LBM[10  + 2*Nfields] = -1.0; LBM[10  + 3*Nfields] = -1.0;
    LBM[11  + 1*Nfields] =  1.0; LBM[11  + 2*Nfields] = -1.0; LBM[11  + 3*Nfields] =  1.0;
    LBM[12  + 1*Nfields] = -1.0; LBM[12  + 2*Nfields] =  1.0; LBM[12  + 3*Nfields] = -1.0;
    LBM[13  + 1*Nfields] =  1.0; LBM[13  + 2*Nfields] =  1.0; LBM[13  + 3*Nfields] = -1.0;
    LBM[14  + 1*Nfields] = -1.0; LBM[14  + 2*Nfields] = -1.0; LBM[14  + 3*Nfields] =  1.0;

  }else if(settings.compareSetting("VELOCITY MODEL", "D3Q19") && mesh.dim==3){
    c       = 1.0/sqrt(3.0);
    Nfields  = 19; 
    Nmacro   = mesh.dim + 1;

    LBM.malloc(Nfields*Nmacro);
    LMAP.malloc(Nfields);
    // Weights
    LBM[0  + 0*Nfields] = 1.0/ 3.0; LMAP[0 ] = 0; 

    LBM[1  + 0*Nfields] = 1.0/ 18.0; LMAP[1 ] = 2; 
    LBM[2  + 0*Nfields] = 1.0/ 18.0; LMAP[2 ] = 1; 
    LBM[3  + 0*Nfields] = 1.0/ 18.0; LMAP[3 ] = 4; 
    LBM[4  + 0*Nfields] = 1.0/ 18.0; LMAP[4 ] = 3; 
    LBM[5  + 0*Nfields] = 1.0/ 18.0; LMAP[5 ] = 6; 
    LBM[6  + 0*Nfields] = 1.0/ 18.0; LMAP[6 ] = 5; 

    LBM[7  + 0*Nfields] = 1.0/36.0; LMAP[7 ] = 8; 
    LBM[8  + 0*Nfields] = 1.0/36.0; LMAP[8 ] = 7; 
    LBM[9  + 0*Nfields] = 1.0/36.0; LMAP[9 ] = 10; 
    LBM[10 + 0*Nfields] = 1.0/36.0; LMAP[10] = 9; 
    LBM[11 + 0*Nfields] = 1.0/36.0; LMAP[11] = 12; 
    LBM[12 + 0*Nfields] = 1.0/36.0; LMAP[12] = 11; 
    LBM[13 + 0*Nfields] = 1.0/36.0; LMAP[13] = 14; 
    LBM[14 + 0*Nfields] = 1.0/36.0; LMAP[14] = 13; 
    LBM[15 + 0*Nfields] = 1.0/36.0; LMAP[15] = 16; 
    LBM[16 + 0*Nfields] = 1.0/36.0; LMAP[16] = 15; 
    LBM[17 + 0*Nfields] = 1.0/36.0; LMAP[17] = 18; 
    LBM[18 + 0*Nfields] = 1.0/36.0; LMAP[18] = 17; 

    // Velocities
    LBM[0  + 1*Nfields] =  0.0; LBM[0  + 2*Nfields] =  0.0; LBM[0  + 3*Nfields] =  0.0;

    LBM[1  + 1*Nfields] =  1.0; LBM[1  + 2*Nfields] =  0.0; LBM[1  + 3*Nfields] =  0.0;
    LBM[2  + 1*Nfields] = -1.0; LBM[2  + 2*Nfields] =  0.0; LBM[2  + 3*Nfields] =  0.0;
    LBM[3  + 1*Nfields] =  0.0; LBM[3  + 2*Nfields] =  1.0; LBM[3  + 3*Nfields] =  0.0;
    LBM[4  + 1*Nfields] =  0.0; LBM[4  + 2*Nfields] = -1.0; LBM[4  + 3*Nfields] =  0.0;
    LBM[5  + 1*Nfields] =  0.0; LBM[5  + 2*Nfields] =  0.0; LBM[5  + 3*Nfields] =  1.0;
    LBM[6  + 1*Nfields] =  0.0; LBM[6  + 2*Nfields] =  0.0; LBM[6  + 3*Nfields] = -1.0;
    // x-y 
    LBM[7   + 1*Nfields] =  1.0; LBM[7   + 2*Nfields] =  1.0; LBM[7   + 3*Nfields] =  0.0;
    LBM[8   + 1*Nfields] = -1.0; LBM[8   + 2*Nfields] = -1.0; LBM[8   + 3*Nfields] =  0.0;
    LBM[9   + 1*Nfields] = -1.0; LBM[9   + 2*Nfields] =  1.0; LBM[9   + 3*Nfields] =  0.0;
    LBM[10  + 1*Nfields] =  1.0; LBM[10  + 2*Nfields] = -1.0; LBM[10  + 3*Nfields] =  0.0;
    // x-z
    LBM[11  + 1*Nfields] =  1.0; LBM[11  + 2*Nfields] =  0.0; LBM[11  + 3*Nfields] =  1.0;
    LBM[12  + 1*Nfields] = -1.0; LBM[12  + 2*Nfields] =  0.0; LBM[12  + 3*Nfields] = -1.0;
    LBM[13  + 1*Nfields] = -1.0; LBM[13  + 2*Nfields] =  0.0; LBM[13  + 3*Nfields] =  1.0;
    LBM[14  + 1*Nfields] =  1.0; LBM[14  + 2*Nfields] =  0.0; LBM[14  + 3*Nfields] = -1.0;
    // y-z
    LBM[15  + 1*Nfields] =  0.0; LBM[15  + 2*Nfields] =  1.0; LBM[15  + 3*Nfields] =  1.0;
    LBM[16  + 1*Nfields] =  0.0; LBM[16  + 2*Nfields] = -1.0; LBM[16  + 3*Nfields] = -1.0;
    LBM[17  + 1*Nfields] =  0.0; LBM[17  + 2*Nfields] = -1.0; LBM[17  + 3*Nfields] =  1.0;
    LBM[18  + 1*Nfields] =  0.0; LBM[18  + 2*Nfields] =  1.0; LBM[18  + 3*Nfields] = -1.0;

  }else {
    LIBP_FORCE_ABORT("Requested VELOCIY MODEL not found.");
  }

  alpha    = nu/(c*c); // Relaxation parameter 
  tauInv   = 1.0/alpha; // need to check that....
  RT       = c*c; // remove this later, no need!!!!!

  // Lattice-Boltzmann Model
  o_LBM = platform.malloc<dfloat>(LBM);
  o_LMAP = platform.malloc<dfloat>(LMAP);
}
