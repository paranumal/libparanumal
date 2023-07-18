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

#include "acoustics.hpp"
#include "timer.hpp"

void acoustics_t::ReportError(dfloat t, dfloat dt, dfloat elapsedTime, int iterations, deviceMemory<dfloat> &o_qin){

  dlong Nall = mesh.Np*mesh.Nelements;
  
  memory<dfloat> qin(Nall*Nfields);
  memory<dfloat> exactq(Nall*Nfields);
  
  deviceMemory<dfloat> o_exactq = platform.malloc<dfloat>(Nall*Nfields);

  initialConditionKernel(mesh.Nelements, t, mesh.o_x, mesh.o_y, mesh.o_z, o_exactq);
  
  o_exactq.copyTo(exactq);
  o_qin.copyTo(qin);
  
  dfloat errRMax = 0;
  dfloat errUMax = 0, errVMax = 0, errWMax = 0;

  for(dlong e=0;e<mesh.Nelements;++e){
    for(int n=0;n<mesh.Np;++n){
      dlong id = e*mesh.Np*Nfields + n;
      dfloat errRn = fabs(qin[id+0*mesh.Np] - exactq[id+0*mesh.Np]);
      dfloat errUn = fabs(qin[id+1*mesh.Np] - exactq[id+1*mesh.Np]);
      dfloat errVn = fabs(qin[id+2*mesh.Np] - exactq[id+2*mesh.Np]);
      dfloat errWn = fabs(qin[id+3*mesh.Np] - exactq[id+3*mesh.Np]);
      
      errRMax = std::max(errRMax, errRn);
      errUMax = std::max(errUMax, errUn);
      errVMax = std::max(errVMax, errVn);
      errWMax = std::max(errWMax, errWn);
    }
  }

  std::cout << std::setprecision (6);
  std::cout << std::scientific;
  std::cout << "errRMax = " << errRMax <<
     " errUMax = " << errUMax <<
     " errVMax = " << errVMax <<
     " errWMax = " << errWMax <<
     std::endl;
  
  dfloat errRL2 = 0, errUL2 = 0;
  dfloat normExactUL2 = 0, normExactRL2 = 0;
  // change for other elements
  if(mesh.elementType==Mesh::TRIANGLES ||
     mesh.elementType==Mesh::TETRAHEDRA){
    for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.Np;++n){
        dlong idn = e*mesh.Np + n;
        dfloat errRL2n = 0, normExactRL2n = 0;
        dfloat errUL2n = 0, normExactUL2n = 0;
        dfloat errVL2n = 0, normExactVL2n = 0;
        dfloat errWL2n = 0, normExactWL2n = 0;
        for(int m=0;m<mesh.Np;++m){
          dlong idm = e*mesh.Np*Nfields + m;
          dfloat MMnm = mesh.MM[n*mesh.Np+m];
          errRL2n += MMnm*(qin[idm+0*mesh.Np]-exactq[idm+0*mesh.Np]);
          errUL2n += MMnm*(qin[idm+1*mesh.Np]-exactq[idm+1*mesh.Np]);
          errVL2n += MMnm*(qin[idm+2*mesh.Np]-exactq[idm+2*mesh.Np]);
          errWL2n += MMnm*(qin[idm+3*mesh.Np]-exactq[idm+3*mesh.Np]);
                              
          normExactRL2n += MMnm*exactq[idm+0*mesh.Np];
          normExactUL2n += MMnm*exactq[idm+1*mesh.Np];
          normExactVL2n += MMnm*exactq[idm+2*mesh.Np];
          normExactWL2n += MMnm*exactq[idm+3*mesh.Np];
        }
        errRL2 += mesh.wJ[e]*(qin[idn+0*mesh.Np]-exactq[idn+0*mesh.Np])*errRL2n;
        errUL2 += mesh.wJ[e]*(qin[idn+1*mesh.Np]-exactq[idn+1*mesh.Np])*errUL2n;
        errUL2 += mesh.wJ[e]*(qin[idn+2*mesh.Np]-exactq[idn+2*mesh.Np])*errVL2n;
        errUL2 += mesh.wJ[e]*(qin[idn+3*mesh.Np]-exactq[idn+3*mesh.Np])*errWL2n;
        normExactRL2 += mesh.wJ[e]*(exactq[idn+0*mesh.Np])*normExactRL2n;
        normExactUL2 += mesh.wJ[e]*(exactq[idn+1*mesh.Np])*normExactUL2n;
        normExactUL2 += mesh.wJ[e]*(exactq[idn+2*mesh.Np])*normExactVL2n;
        normExactUL2 += mesh.wJ[e]*(exactq[idn+3*mesh.Np])*normExactWL2n;
      }
    }
  }
  
  if(mesh.elementType==Mesh::QUADRILATERALS ||
     mesh.elementType==Mesh::HEXAHEDRA){
    for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.Np;++n){
        dlong gidn = e*mesh.Np + n;
        dfloat wJn = mesh.wJ[gidn];
        dlong idn = e*mesh.Np*Nfields + n;
        errRL2 += wJn*pow(qin[idn+0*mesh.Np]-exactq[idn+0*mesh.Np],2);
        errUL2 += wJn*pow(qin[idn+1*mesh.Np]-exactq[idn+1*mesh.Np],2);
        errUL2 += wJn*pow(qin[idn+2*mesh.Np]-exactq[idn+2*mesh.Np],2);
        errUL2 += wJn*pow(qin[idn+3*mesh.Np]-exactq[idn+3*mesh.Np],2);
        
        normExactRL2 += wJn*pow(exactq[idn+0*mesh.Np],2);
        normExactUL2 += wJn*pow(exactq[idn+1*mesh.Np],2);
        normExactUL2 += wJn*pow(exactq[idn+2*mesh.Np],2);
        normExactUL2 += wJn*pow(exactq[idn+3*mesh.Np],2);
      }
    }
  }

  dfloat eps = 1e-16;
  dfloat relErrRL2 = sqrt(errRL2/(eps+normExactRL2));
  dfloat relErrUL2 = sqrt(errUL2/(eps+normExactUL2));
  
  std::cout << t << "," << mesh.N << ", " << mesh.Nelements << ", " << dt << ", " << iterations << ", " << relErrRL2 << ", " <<  relErrUL2 << ", " << elapsedTime <<
     "; %% t, N, Nelements, dt, relErrRL2, relErrUL2, elapsedTime" << std::endl;

}
