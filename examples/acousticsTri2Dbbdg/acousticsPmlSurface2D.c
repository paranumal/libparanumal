#include <math.h>
#include <stdlib.h>
#include "mesh2D.h"

void boundaryConditions2D(iint bc, dfloat t, dfloat x, dfloat y,
                          dfloat uM, dfloat vM, dfloat pM,
                          dfloat *uP, dfloat *vP, dfloat *pP);
  
void acousticsPmlSurface2Dbbdg(mesh2D *mesh, iint lev, dfloat t){

  // temporary storage for flux terms
  dfloat *fluxu = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxv = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxpx = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxpy = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));

  dfloat *fluxu_copy = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxv_copy = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxpx_copy = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxpx_copy = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));

  // for all elements
  for(iint et=0;et<mesh->MRABpmlNelements[lev];++et){
    iint e = mesh->MRABpmlElementIds[lev][et];
    iint pmlId = mesh->MRABpmlElementIds[lev][et];
    // for all face nodes of all elements
    for(iint n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      // find face that owns this node
      iint face = n/mesh->Nfp;

      // load surface geofactors for this face
      iint sid = mesh->Nsgeo*(e*mesh->Nfaces+face);
      dfloat nx = mesh->sgeo[sid+0];
      dfloat ny = mesh->sgeo[sid+1];
      dfloat sJ = mesh->sgeo[sid+2];
      dfloat invJ = mesh->sgeo[sid+3];

      iint id = n + e*mesh->Nfaces*mesh->Nfp;
      iint idM = id*mesh->Nfields;
      iint idP = mesh->mapP[id]*mesh->Nfields;

      // load negative trace node values of q
      dfloat uM = mesh->fQM[idM+0];
      dfloat vM = mesh->fQM[idM+1];
      dfloat pM = mesh->fQM[idM+2];

      // load positive trace node values of q
      dfloat uP = mesh->fQP[idP+0]; 
      dfloat vP = mesh->fQP[idP+1];
      dfloat pP = mesh->fQP[idP+2];

      // find boundary type
      iint boundaryType = mesh->EToB[e*mesh->Nfaces+face];
      if(boundaryType>0) {
        iint idM = mesh->vmapM[id];
        boundaryConditions2D(boundaryType, t, mesh->x[idM], mesh->y[idM], 
                              uM, vM, pM, &uP, &vP, &pP);
      }

      // compute (q^* - q^-)
      dfloat duS = 0.5f*(uP-uM) + mesh->Lambda2*(-nx)*(pP-pM);
      dfloat dvS = 0.5f*(vP-vM) + mesh->Lambda2*(-ny)*(pP-pM);
      dfloat dpxS = 0.5f*(pP-pM) + mesh->Lambda2*(-nx*(uP-uM));
      dfloat dpyS = 0.5f*(pP-pM) + mesh->Lambda2*(-ny*(vP-vM));

      // evaluate "flux" terms: (sJ/J)*(A*nx+B*ny)*(q^* - q^-)
      fluxu[n] = invJ*sJ*(-nx*dpS);
      fluxv[n] = invJ*sJ*(-ny*dpS);
      fluxpx[n] = invJ*sJ*(-nx*duS);
      fluxpy[n] = invJ*sJ*(-ny*dvS);
    }

    // apply L0 to fluxes. use fact that L0 = tridiagonal in 2D
    for(iint n=0;n<mesh->Nfp*mesh->Nfaces;++n){
    
      iint id = n % mesh->Nfp;  // warning: redundant reads
      dfloat L0val = mesh->L0vals[3*id+1]; 

      dfloat utmpflux = L0val * fluxu[n];
      dfloat vtmpflux = L0val * fluxv[n];
      dfloat pxtmpflux = L0val * fluxpx[n];
      dfloat pytmpflux = L0val * fluxpy[n];

      if (id > 0){     
        utmpflux += mesh->L0vals[3*id]*fluxu[n-1]; // add previous term
        vtmpflux += mesh->L0vals[3*id]*fluxv[n-1]; // add previous term
        pxtmpflux += mesh->L0vals[3*id]*fluxpx[n-1]; // add previous term
        pytmpflux += mesh->L0vals[3*id]*fluxpy[n-1]; // add previous term
      }
      if (id < mesh->Nfp-1){
        utmpflux += mesh->L0vals[3*id+2]*fluxu[n+1];// add next term
        vtmpflux += mesh->L0vals[3*id+2]*fluxv[n+1];// add next term
        pxtmpflux += mesh->L0vals[3*id+2]*fluxpx[n+1];// add next term
        pytmpflux += mesh->L0vals[3*id+2]*fluxpy[n+1];// add next term
      }
      fluxu_copy[n] = utmpflux;
      fluxv_copy[n] = vtmpflux;
      fluxpx_copy[n] = pxtmpflux;
      fluxpy_copy[n] = pytmpflux;
    }

    // apply lift reduction and accumulate RHS
    for(iint n=0;n<mesh->Np;++n){
      iint rhsId = 3*mesh->Nfields*(mesh->Np*e + n) + mesh->Nfields*mesh->MRABshiftIndex[lev];
      iint pmlrhsId = 3*mesh->pmlNfields*(mesh->Np*pmlId + n) + mesh->pmlNfields*mesh->MRABshiftIndex[lev];

      // load RHS
      dfloat rhsqnu = mesh->rhsq[rhsId+0];
      dfloat rhsqnv = mesh->rhsq[rhsId+1];
      dfloat rhsqnpx = mesh->pmlrhsq[pmlrhsId+0];
      dfloat rhsqnpy = mesh->pmlrhsq[pmlrhsId+1];

      // rhs += LIFT*((sJ/J)*(A*nx+B*ny)*(q^* - q^-))
      for (int m = 0; m < mesh->max_EL_nnz; ++m){
        iint idm = m + n*mesh->max_EL_nnz;
        dfloat ELval = mesh->ELvals[idm];
        iint ELid = mesh->ELids[idm];
        rhsqnu += ELval * fluxu_copy[ELid];
        rhsqnv += ELval * fluxv_copy[ELid];
        rhsqnpx += ELval * fluxpx_copy[ELid];
        rhsqnpy += ELval * fluxpy_copy[ELid];
      }
      
      // store incremented rhs
      mesh->rhsq[rhsId+0] = rhsqnu;
      mesh->rhsq[rhsId+1] = rhsqnv;
      mesh->pmlrhsq[pmlrhsId+0] = rhsqnpx;  
      mesh->pmlrhsq[pmlrhsId+1] = rhsqnpy;  
    }
  }

  // free temporary flux storage
  free(fluxu); free(fluxv); free(fluxpx); free(fluxpy);
  free(fluxu_copy); free(fluxv_copy); free(fluxpx_copy); free(fluxpy_copy);
}
