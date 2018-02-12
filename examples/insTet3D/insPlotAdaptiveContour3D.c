#include "ins3D.h"

// interpolate data to plot nodes and save to file (one per process
void insPlotAdaptiveContour3D(ins_t *ins, char *fileName){

  mesh3D *mesh = ins->mesh;

  int Nlevels = 10;
  dfloat levels[10] = {1.,2.,3.,4.,5.,6.,7.,8.,9.,10.};
  dfloat tol = 1E-3;

  dfloat *Vort = (dfloat*) calloc(mesh->Np*mesh->Nelements,sizeof(dfloat));

  // calculate vorticity magnitude
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      dfloat dUdr = 0, dUds = 0, dUdt = 0 ;
      dfloat dVdr = 0, dVds = 0, dVdt = 0 ;
      dfloat dWdr = 0, dWds = 0, dWdt = 0 ; 
      for(iint m=0;m<mesh->Np;++m){
        iint id = m+e*mesh->Np;
        id += ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);

        dUdr += mesh->Dr[n*mesh->Np+m]*ins->U[id];
        dUds += mesh->Ds[n*mesh->Np+m]*ins->U[id];
        dUdt += mesh->Dt[n*mesh->Np+m]*ins->U[id];

        dVdr += mesh->Dr[n*mesh->Np+m]*ins->V[id];
        dVds += mesh->Ds[n*mesh->Np+m]*ins->V[id];
        dVdt += mesh->Dt[n*mesh->Np+m]*ins->V[id];

        dWdr += mesh->Dr[n*mesh->Np+m]*ins->W[id];
        dWds += mesh->Ds[n*mesh->Np+m]*ins->W[id];
        dWdt += mesh->Dt[n*mesh->Np+m]*ins->W[id];
      }

      dfloat rx = mesh->vgeo[e*mesh->Nvgeo+RXID];
      dfloat ry = mesh->vgeo[e*mesh->Nvgeo+RYID];
      dfloat rz = mesh->vgeo[e*mesh->Nvgeo+RZID];    
      
      dfloat sx = mesh->vgeo[e*mesh->Nvgeo+SXID];
      dfloat sy = mesh->vgeo[e*mesh->Nvgeo+SYID];
      dfloat sz = mesh->vgeo[e*mesh->Nvgeo+SZID];    
     
      dfloat tx = mesh->vgeo[e*mesh->Nvgeo+TXID];
      dfloat ty = mesh->vgeo[e*mesh->Nvgeo+TYID];
      dfloat tz = mesh->vgeo[e*mesh->Nvgeo+TZID];    

      dfloat dUdx = rx*dUdr + sx*dUds + tx*dUdt;
      dfloat dUdy = ry*dUdr + sy*dUds + ty*dUdt;
      dfloat dUdz = rz*dUdr + sz*dUds + tz*dUdt;
    
      dfloat dVdx = rx*dVdr + sx*dVds + tx*dVdt;
      dfloat dVdy = ry*dVdr + sy*dVds + ty*dVdt;
      dfloat dVdz = rz*dVdr + sz*dVds + tz*dVdt;
      
      dfloat dWdx = rx*dWdr + sx*dWds + tx*dWdt;
      dfloat dWdy = ry*dWdr + sy*dWds + ty*dWdt;
      dfloat dWdz = rz*dWdr + sz*dWds + tz*dWdt;
      
      // Compute vorticity Vector
      dfloat Vx = dWdy-dVdz;
      dfloat Vy = dUdz-dWdx;
      dfloat Vz = dVdx-dUdy;

      Vort[e*mesh->Np+n] = sqrt(Vx*Vx+Vy*Vy+Vz*Vz);
    }
  }
  
  meshPlotAdaptiveContour3D(mesh, fileName, Vort, Nlevels, levels, tol);
}
