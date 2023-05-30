
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

#include "wave.hpp"
#include "timer.hpp"
#include "ellipticPrecon.hpp"

void wave_t::waveHoltz(deviceMemory<dfloat> &o_rDL,
                       deviceMemory<dfloat> &o_rPL,
                       deviceMemory<dfloat> &o_rFL){
  
  // start with zeroed P, L, and filtered pressure
  platform.linAlg().set(Nall, (dfloat)0, o_DL);
  platform.linAlg().set(Nall, (dfloat)0, o_PL);

  // arbitrarily try 100 outer iterations
  for(int it=1;it<=100;++it){

    std::cout << "STARTING OUTER ITERATION: " << it << std::endl;
    
    platform.linAlg().set(Nall, (dfloat)0, o_DL);
    platform.linAlg().set(Nall, (dfloat)0, o_FPL);

    // solve to T = 2pi/omega
    Solve(o_DL, o_PL, o_FL);

    // PL <= FPL
    o_FPL.copyTo(o_PL);
    
    if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {
      static int slice=0;
      ++slice;
      
      // copy data back to host
      o_PL.copyTo(PL);
      o_DL.copyTo(DL);
      
      // output field files
      std::string name;
      settings.getSetting("OUTPUT FILE NAME", name);
      char fname[BUFSIZ];
      sprintf(fname, "ITER_DP_%04d_%04d.vtu",  mesh.rank, slice);
      PlotFields(DL, PL, fname);
    }
  }
}
