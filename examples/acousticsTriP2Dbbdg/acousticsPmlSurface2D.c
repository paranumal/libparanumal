#include "acoustics2D.h"

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
  dfloat *fluxpy_copy = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));

  // for all elements
  for(iint et=0;et<mesh->MRABpmlNelements[lev];++et){
    iint e = mesh->MRABpmlElementIds[lev][et];
    iint pmlId = mesh->MRABpmlIds[lev][et];
    iint N = mesh->N[e];
    // for all face nodes of all elements
    for (iint f=0;f<mesh->Nfaces;f++) {
      for(iint n=0;n<mesh->Nfp[N];++n){
        // load surface geofactors for this face
        iint sid = mesh->Nsgeo*(e*mesh->Nfaces+f);
        dfloat nx   = mesh->sgeo[sid+NXID];
        dfloat ny   = mesh->sgeo[sid+NYID];
        dfloat sJ   = mesh->sgeo[sid+SJID];
        dfloat invJ = mesh->sgeo[sid+IJID];

        iint id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax + n;
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
        iint boundaryType = mesh->EToB[e*mesh->Nfaces+f];
        if(boundaryType>0) {
          iint idM = mesh->vmapM[id];
          boundaryConditions2D(boundaryType, t, mesh->x[idM], mesh->y[idM],
                                uM, vM, pM, &uP, &vP, &pP);
        }

        // compute (q^* - q^-)
        dfloat duS = 0.5f*(uP-uM) + mesh->Lambda2*(-nx)*(pP-pM);
        dfloat dvS = 0.5f*(vP-vM) + mesh->Lambda2*(-ny)*(pP-pM);
        dfloat dpS = 0.5f*(pP-pM) + mesh->Lambda2*(-nx*(uP-uM)-ny*(vP-vM));

        // evaluate "flux" terms: (sJ/J)*(A*nx+B*ny)*(q^* - q^-)
        fluxu[f*mesh->Nfp[N] +n] = invJ*sJ*(-nx*dpS);
        fluxv[f*mesh->Nfp[N] +n] = invJ*sJ*(-ny*dpS);
        fluxpx[f*mesh->Nfp[N] +n] = invJ*sJ*(-nx*duS);
        fluxpy[f*mesh->Nfp[N] +n] = invJ*sJ*(-ny*dvS);
      }
    }

    // apply L0 to fluxes. use fact that L0 = tridiagonal in 2D
    for(iint n=0;n<mesh->Nfp[N]*mesh->Nfaces;++n){

      iint id = n % mesh->Nfp[N];  // warning: redundant reads
      dfloat L0val = mesh->L0vals[N][3*id+1];

      dfloat utmpflux = L0val * fluxu[n];
      dfloat vtmpflux = L0val * fluxv[n];
      dfloat pxtmpflux = L0val * fluxpx[n];
      dfloat pytmpflux = L0val * fluxpy[n];

      if (id > 0){
        utmpflux += mesh->L0vals[N][3*id]*fluxu[n-1]; // add previous term
        vtmpflux += mesh->L0vals[N][3*id]*fluxv[n-1]; // add previous term
        pxtmpflux += mesh->L0vals[N][3*id]*fluxpx[n-1]; // add previous term
        pytmpflux += mesh->L0vals[N][3*id]*fluxpy[n-1]; // add previous term
      }
      if (id < mesh->Nfp[N]-1){
        utmpflux += mesh->L0vals[N][3*id+2]*fluxu[n+1];// add next term
        vtmpflux += mesh->L0vals[N][3*id+2]*fluxv[n+1];// add next term
        pxtmpflux += mesh->L0vals[N][3*id+2]*fluxpx[n+1];// add next term
        pytmpflux += mesh->L0vals[N][3*id+2]*fluxpy[n+1];// add next term
      }
      fluxu_copy[n] = utmpflux;
      fluxv_copy[n] = vtmpflux;
      fluxpx_copy[n] = pxtmpflux;
      fluxpy_copy[n] = pytmpflux;
    }

    // apply lift reduction and accumulate RHS
    for(iint n=0;n<mesh->Np[N];++n){
      iint rhsId = 3*mesh->Nfields*(mesh->NpMax*e + n) + mesh->Nfields*mesh->MRABshiftIndex[lev];
      iint pmlrhsId = 3*mesh->pmlNfields*(mesh->NpMax*pmlId + n) + mesh->pmlNfields*mesh->MRABshiftIndex[lev];

      // load RHS
      dfloat rhsqnu = mesh->rhsq[rhsId+0];
      dfloat rhsqnv = mesh->rhsq[rhsId+1];
      dfloat rhsqnpx = mesh->pmlrhsq[pmlrhsId+0];
      dfloat rhsqnpy = mesh->pmlrhsq[pmlrhsId+1];

      // rhs += LIFT*((sJ/J)*(A*nx+B*ny)*(q^* - q^-))
      for (int m = 0; m < mesh->max_EL_nnz[N]; ++m){
        iint idm = m + n*mesh->max_EL_nnz[N];
        dfloat ELval = mesh->ELvals[N][idm];
        iint ELid = mesh->ELids[N][idm];
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
