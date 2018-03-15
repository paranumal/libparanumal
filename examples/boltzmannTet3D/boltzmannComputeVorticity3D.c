#include "boltzmann3D.h"

void boltzmannComputeVorticity3D(mesh3D *mesh, dfloat *q, int Nfields){
  
  // compute vorticity
  dfloat *u = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *v = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *w = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  //
  for(int e=0;e<mesh->Nelements;++e){
    
    for(int n=0;n<mesh->Np;++n){
      int base = mesh->Nfields*(n + e*mesh->Np);
      dfloat rho = mesh->q[base + 0];
      u[n] = mesh->q[1 + base]*mesh->sqrtRT/rho;
      v[n] = mesh->q[2 + base]*mesh->sqrtRT/rho;
      w[n] = mesh->q[3 + base]*mesh->sqrtRT/rho;
    }

    for(int n=0;n<mesh->Np;++n){
      dfloat dudr = 0, duds = 0, dudt = 0; 
      dfloat dvdr = 0, dvds = 0, dvdt =0 ;
      dfloat dwdr = 0, dwds = 0, dwdt =0 ;
      for(int m=0;m<mesh->Np;++m){

      	dudr += mesh->Dr[n*mesh->Np+m]*u[m];
      	duds += mesh->Ds[n*mesh->Np+m]*u[m];
        dudt += mesh->Dt[n*mesh->Np+m]*u[m];
        //
      	dvdr += mesh->Dr[n*mesh->Np+m]*v[m];
        dvds += mesh->Ds[n*mesh->Np+m]*v[m];
        dvdt += mesh->Dt[n*mesh->Np+m]*v[m];
        //
        dwdr += mesh->Dr[n*mesh->Np+m]*w[m];
        dwds += mesh->Ds[n*mesh->Np+m]*w[m];
        dwdt += mesh->Dt[n*mesh->Np+m]*w[m];
      }

      dfloat dwdy = mesh->vgeo[mesh->Nvgeo*e+RYID]*dwdr + 
                    mesh->vgeo[mesh->Nvgeo*e+SYID]*dwds +
                    mesh->vgeo[mesh->Nvgeo*e+TYID]*dwdt ;

      dfloat dvdz = mesh->vgeo[mesh->Nvgeo*e+RZID]*dvdr + 
                    mesh->vgeo[mesh->Nvgeo*e+SZID]*dvds +
                    mesh->vgeo[mesh->Nvgeo*e+TZID]*dvdt ;

      q[4 + Nfields*(n+e*mesh->Np)] = dwdy-dvdz;

      dfloat dwdx = mesh->vgeo[mesh->Nvgeo*e+RXID]*dwdr + 
                    mesh->vgeo[mesh->Nvgeo*e+SXID]*dwds +
                    mesh->vgeo[mesh->Nvgeo*e+TXID]*dwdt ;

      dfloat dudz = mesh->vgeo[mesh->Nvgeo*e+RZID]*dudr + 
                    mesh->vgeo[mesh->Nvgeo*e+SZID]*duds +
                    mesh->vgeo[mesh->Nvgeo*e+TZID]*dudt;

      q[5 + Nfields*(n+e*mesh->Np)] = -dwdx + dudz;

      dfloat dvdx = mesh->vgeo[mesh->Nvgeo*e+RXID]*dvdr + 
                    mesh->vgeo[mesh->Nvgeo*e+SXID]*dvds +
                    mesh->vgeo[mesh->Nvgeo*e+TXID]*dvdt;

      dfloat dudy = mesh->vgeo[mesh->Nvgeo*e+RYID]*dudr + 
                    mesh->vgeo[mesh->Nvgeo*e+SYID]*duds +
                    mesh->vgeo[mesh->Nvgeo*e+TYID]*dudt;
                    
      q[6 + Nfields*(n+e*mesh->Np)] = dvdx-dudy;
    }
  }
  free(u);
  free(v);
  free(w);

}
