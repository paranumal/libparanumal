#include <stdio.h>
#include <math.h>
#include "mesh2D.h"

void acousticsComputeVorticity2D(mesh2D *mesh, dfloat *q, int outfld, int Nfields){
  
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
	mesh->vgeo[mesh->Nvgeo*e+RYID]*dvdr +
	mesh->vgeo[mesh->Nvgeo*e+SYID]*dvds;
      q[outfld + Nfields*(n+e*mesh->Np)] = dudy-dvdx;
    }
  }
  free(u);
  free(v);
}
