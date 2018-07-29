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

#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshGradient2D(mesh2D *mesh,
		    dfloat *q,
		    dfloat *dqdx,
		    dfloat *dqdy
		    ){
  
  int cnt = 0;
  for(int e=0;e<mesh->Nelements;++e){
    
    dfloat drdx = mesh->vgeo[e*mesh->Nvgeo + RXID];
    dfloat drdy = mesh->vgeo[e*mesh->Nvgeo + RYID];
    dfloat dsdx = mesh->vgeo[e*mesh->Nvgeo + SXID];
    dfloat dsdy = mesh->vgeo[e*mesh->Nvgeo + SYID];
    
    for(int n=0;n<mesh->Np;++n){
      dfloat dqdr = 0, dqds = 0;
      
      for(int m=0;m<mesh->Np;++m){
	
	dqdr += mesh->Dr[n*mesh->Np + m]*q[m + e*mesh->Np];
	dqds += mesh->Ds[n*mesh->Np + m]*q[m + e*mesh->Np];
      }
      
      dqdx[n+e*mesh->Np] = drdx*dqdr + dsdx*dqds;
      dqdy[n+e*mesh->Np] = drdy*dqdr + dsdy*dqds;
    }
  }
}
