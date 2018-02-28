#include "acoustics3D.h"

void acousticsUpdate3D(mesh3D *mesh, dfloat rka, dfloat rkb){
  
  // Low storage Runge Kutta time step update
  for(int n=0;n<mesh->Nelements*mesh->Np*mesh->Nfields;++n){

    mesh->resq[n] = rka*mesh->resq[n] + mesh->dt*mesh->rhsq[n];
    
    mesh->q[n] += rkb*mesh->resq[n];
  }
}

