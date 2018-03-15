#include "acoustics2D.h"

void acousticsPmlUpdate2D(mesh2D *mesh, dfloat rka, dfloat rkb){
  
  // Low storage Runge Kutta time step update
  for(iint n=0;n<mesh->pmlNelements*mesh->Np*mesh->pmlNfields;++n){

    mesh->pmlresq[n] = rka*mesh->pmlresq[n] + mesh->dt*mesh->pmlrhsq[n];
    
    mesh->pmlq[n] += rkb*mesh->pmlresq[n];
  }
}

