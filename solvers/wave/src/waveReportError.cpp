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

void wave_t::ReportError(dfloat t, dfloat elapsedTime, deviceMemory<dfloat> &o_DLin, deviceMemory<dfloat> &o_PLin){

  memory<dfloat> DLin(Nall);
  memory<dfloat> PLin(Nall);
  
  deviceMemory<dfloat> o_exactDL = platform.malloc<dfloat>(Nall);
  deviceMemory<dfloat> o_exactPL = platform.malloc<dfloat>(Nall);

  memory<dfloat> exactDL(Nall);
  memory<dfloat> exactPL(Nall);

  waveInitialConditionsKernel(Nall, t, mesh.o_x, mesh.o_y, mesh.o_z, o_exactDL, o_exactPL);
  
  o_exactDL.copyTo(exactDL);
  o_exactPL.copyTo(exactPL);

  o_DLin.copyTo(DLin);
  o_PLin.copyTo(PLin);
  
  dfloat errDLMax = 0;
  dfloat errPLMax = 0;
  for(int n=0;n<Nall;++n){
    dfloat errDLn = fabs(DLin[n] - exactDL[n]);
    dfloat errPLn = fabs(PLin[n] - exactPL[n]);
    
    errDLMax = std::max(errDLMax, errDLn);
    errPLMax = std::max(errPLMax, errPLn);
  }

  std::cout << std::setprecision (6);
  std::cout << std::scientific;
  std::cout << "errDLMax = " << errDLMax << " errPLMax = " << errPLMax << std::endl;
  
  dfloat errPL2 = 0, errDL2 = 0;
  dfloat normExactPL2 = 0, normExactDL2 = 0;
  // change for other elements
  if(mesh.elementType==Mesh::TRIANGLES ||
     mesh.elementType==Mesh::TETRAHEDRA){
    for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.Np;++n){
        dlong idn = e*mesh.Np + n;
        dfloat errPL2n = 0, normExactPL2n = 0;
        dfloat errDL2n = 0, normExactDL2n = 0;
        for(int m=0;m<mesh.Np;++m){
          dlong idm = e*mesh.Np + m;
          dfloat MMnm = mesh.MM[n*mesh.Np+m];
          errPL2n += MMnm*(PLin[idm]-exactPL[idm]);
          errDL2n += MMnm*(DLin[idm]-exactDL[idm]);
          normExactPL2n += MMnm*exactPL[idm];
          normExactDL2n += MMnm*exactDL[idm];
        }
        errPL2 += WJ[e]*(PLin[idn]-exactPL[idn])*errPL2n;
        errDL2 += WJ[e]*(DLin[idn]-exactDL[idn])*errDL2n;
        normExactPL2 += WJ[e]*(exactPL[idn])*normExactPL2n;
        normExactDL2 += WJ[e]*(exactDL[idn])*normExactDL2n;
      }
    }
  }
  
  if(mesh.elementType==Mesh::QUADRILATERALS ||
     mesh.elementType==Mesh::HEXAHEDRA){
    for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.Np;++n){
        dlong idn = e*mesh.Np + n;
        dfloat WJn = WJ[idn];
        errPL2 += WJn*pow(PLin[idn]-exactPL[idn],2);
        errDL2 += WJn*pow(DLin[idn]-exactDL[idn],2);
        
        normExactPL2 += WJn*pow(exactPL[idn],2);
        normExactDL2 += WJn*pow(exactDL[idn],2);
      }
    }
  }

  dfloat eps = 1e-16;
  dfloat relErrDL2 = sqrt(errDL2/(eps+normExactDL2));
  dfloat relErrPL2 = sqrt(errPL2/(eps+normExactPL2));
  
  std::cout << t << "," << mesh.N << ", " << mesh.Nelements << ", " << dt << ", " << relErrDL2 << ", " <<  relErrPL2 << ", " << elapsedTime <<
     "; %% t, N, Nelements, dt, relErrDL2, relErrPL2, elapsedTime" << std::endl;

}
