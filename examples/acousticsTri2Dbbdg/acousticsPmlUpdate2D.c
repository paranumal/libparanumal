#include "acoustics2D.h"

void acousticsMRABpmlUpdate2D(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, iint lev, dfloat dt){

  for(iint m=0;m<mesh->MRABpmlNelements[lev];m++){
    iint pmlId = mesh->MRABpmlIds[lev][m];
    for(iint n=0;n<mesh->Np;++n){
      iint id = mesh->pmlNfields*(pmlId*mesh->Np + n);

      iint rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->pmlNfields;
      iint rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->pmlNfields;
      iint rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->pmlNfields;
      for (iint fld=0;fld<mesh->pmlNfields;fld++)
        mesh->pmlq[id+fld] += dt*(a1*mesh->pmlrhsq[rhsId1+fld] + a2*mesh->pmlrhsq[rhsId2+fld] + a3*mesh->pmlrhsq[rhsId3+fld]);
    }
  }
}