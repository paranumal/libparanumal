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

#include "boltzmann2D.h"

void boltzmannComputeVorticity2D(mesh2D *mesh, dfloat *q, int outfld, int Nfields){
  
  // compute vorticity
  dfloat *u = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *v = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(int e=0;e<mesh->Nelements;++e){
    
    for(int n=0;n<mesh->Np;++n){
      int base = mesh->Nfields*(n + e*mesh->Np);
      dfloat rho = mesh->q[base];
      u[n] = mesh->q[1 + base]*mesh->sqrtRT/rho;
      v[n] = mesh->q[2 + base]*mesh->sqrtRT/rho;
    }
    for(int n=0;n<mesh->Np;++n){
      dfloat dudr = 0, duds = 0, dvdr = 0, dvds = 0;
      for(int m=0;m<mesh->Np;++m){
	dudr += mesh->Dr[n*mesh->Np+m]*u[m];
	duds += mesh->Ds[n*mesh->Np+m]*u[m];
	dvdr += mesh->Dr[n*mesh->Np+m]*v[m];
	dvds += mesh->Ds[n*mesh->Np+m]*v[m];
      }
      dfloat dudy =
	mesh->vgeo[mesh->Nvgeo*e+RYID]*dudr +
	mesh->vgeo[mesh->Nvgeo*e+SYID]*duds;
      dfloat dvdx =
	mesh->vgeo[mesh->Nvgeo*e+RXID]*dvdr +
	mesh->vgeo[mesh->Nvgeo*e+SXID]*dvds;
      q[outfld + Nfields*(n+e*mesh->Np)] = dudy-dvdx;
    }
  }
  free(u);
  free(v);
}
