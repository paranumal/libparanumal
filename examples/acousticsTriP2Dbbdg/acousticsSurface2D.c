#include "acoustics2D.h"

void boundaryConditions2D(iint bc, dfloat t, dfloat x, dfloat y,
      dfloat uM, dfloat vM, dfloat pM,
            dfloat *uP, dfloat *vP, dfloat *pP){
  if(1){//bc==1){
    *uP = -uM;  
    *vP = -vM;  
    *pP = pM; 
  }   
  if(0){ // (bc==2){  
    dfloat dx = 1.f/sqrt(2.f);
    dfloat dy = 1.f/sqrt(2.f);
    dfloat omega = 10.f*M_PI;
    dfloat wave = cos(omega*(t-(x*dx+y*dy))); 
    dfloat uI = dx*wave;
    dfloat vI = dy*wave;
    dfloat pI = wave; 
    *uP = -uM -2.f*uI;  
    *vP = -vM -2.f*vI;
    *pP = pM;   
  }
}

// function to compute surface contributions 
// for nodal DG acoustics right hand side
void acousticsSurface2Dbbdg(mesh2D *mesh, iint lev, dfloat t){

  // temporary storage for flux terms
  dfloat *fluxu = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxv = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxp = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));

  dfloat *fluxu_copy = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxv_copy = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxp_copy = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));

  // for all elements
  for(iint et=0;et<mesh->MRABNelements[lev];++et){
    iint e = mesh->MRABelementIds[lev][et];
    iint N = mesh->N[e];
    
    // for all face nodes of all elements
    for (iint f=0;f<mesh->Nfaces;f++) {
      for(iint n=0;n<mesh->Nfp[N];++n){
        // load surface geofactors for this face
        iint sid = mesh->Nsgeo*(e*mesh->Nfaces+f);
        dfloat nx   = mesh->sgeo[sid+0];
        dfloat ny   = mesh->sgeo[sid+1];
        dfloat sJ   = mesh->sgeo[sid+2];
        dfloat invJ = mesh->sgeo[sid+3];

        // indices of negative and positive traces of face node
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
          boundaryConditions2D(boundaryType, t, mesh->x[idM], mesh->y[idM], uM, vM, pM, &uP, &vP, &pP);
        }

        // compute (q^* - q^-)
        dfloat duS = 0.5f*(uP-uM) + mesh->Lambda2*(-nx)*(pP-pM);
        dfloat dvS = 0.5f*(vP-vM) + mesh->Lambda2*(-ny)*(pP-pM);
        dfloat dpS = 0.5f*(pP-pM) + mesh->Lambda2*(-nx*(uP-uM)-ny*(vP-vM));

        // evaluate "flux" terms: (sJ/J)*(A*nx+B*ny)*(q^* - q^-)
        fluxu[f*mesh->Nfp[N] +n] = invJ*sJ*(-nx*dpS);
        fluxv[f*mesh->Nfp[N] +n] = invJ*sJ*(-ny*dpS);
        fluxp[f*mesh->Nfp[N] +n] = invJ*sJ*(-nx*duS-ny*dvS);
      }
    }

    // apply L0 to fluxes. use fact that L0 = tridiagonal in 2D
    for(iint n=0;n<mesh->Nfp[N]*mesh->Nfaces;++n){
    
      iint id = n % mesh->Nfp[N];  // warning: redundant reads
      dfloat L0val = mesh->L0vals[N][3*id+1]; 

      dfloat utmpflux = L0val * fluxu[n];
      dfloat vtmpflux = L0val * fluxv[n];
      dfloat ptmpflux = L0val * fluxp[n];

      if (id > 0){     
        utmpflux += mesh->L0vals[N][3*id]*fluxu[n-1]; // add previous term
        vtmpflux += mesh->L0vals[N][3*id]*fluxv[n-1]; // add previous term
        ptmpflux += mesh->L0vals[N][3*id]*fluxp[n-1]; // add previous term
      }
      if (id < mesh->Nfp[N]-1){
        utmpflux += mesh->L0vals[N][3*id+2]*fluxu[n+1];// add next term
        vtmpflux += mesh->L0vals[N][3*id+2]*fluxv[n+1];// add next term
        ptmpflux += mesh->L0vals[N][3*id+2]*fluxp[n+1];// add next term
      }
      fluxu_copy[n] = utmpflux;
      fluxv_copy[n] = vtmpflux;
      fluxp_copy[n] = ptmpflux;
    }

    // apply lift reduction and accumulate RHS
    for(iint n=0;n<mesh->Np[N];++n){
      iint id = 3*mesh->Nfields*(mesh->NpMax*e + n) + mesh->Nfields*mesh->MRABshiftIndex[lev];
      
      // load RHS
      dfloat rhsqnu = mesh->rhsq[id+0];
      dfloat rhsqnv = mesh->rhsq[id+1];
      dfloat rhsqnp = mesh->rhsq[id+2];

      // rhs += LIFT*((sJ/J)*(A*nx+B*ny)*(q^* - q^-))
      for (int m = 0; m < mesh->max_EL_nnz[N]; ++m){
        iint id = m + n*mesh->max_EL_nnz[N];
        dfloat ELval = mesh->ELvals[N][id];
        iint ELid = mesh->ELids[N][id];
        rhsqnu += ELval * fluxu_copy[ELid];
        rhsqnv += ELval * fluxv_copy[ELid];
        rhsqnp += ELval * fluxp_copy[ELid];
      }
      
      // store incremented rhs
      mesh->rhsq[id+0] = rhsqnu;
      mesh->rhsq[id+1] = rhsqnv;
      mesh->rhsq[id+2] = rhsqnp;  
    }
  }

  // free temporary flux storage
  free(fluxu); free(fluxv); free(fluxp);
  free(fluxu_copy); free(fluxv_copy); free(fluxp_copy);
}
