#include "acoustics3D.h"

void boundaryConditions3D(int bc, dfloat time, dfloat x, dfloat y, dfloat z,
			  dfloat uM, dfloat vM, dfloat wM, dfloat pM,
			  dfloat *uP, dfloat *vP, dfloat *wP, dfloat *pP){
  if(1){//bc==1){
    *uP = -uM;	
    *vP = -vM;
    *wP = -wM;
    *pP = pM;	
  }		
  if(0){//bc==2){	
    dfloat dx = 1.f/sqrt(2.f);
    dfloat dy = 1.f/sqrt(2.f);
    dfloat dz = 0;
    dfloat omega = 10.f*M_PI;
    dfloat wave = cos(omega*(time-(x*dx+y*dy+z*dz)));	
    dfloat uI = dx*wave;
    dfloat vI = dy*wave;
    dfloat wI = dz*wave;
    dfloat pI = wave;	
    *uP = -uM -2.f*uI;	
    *vP = -vM -2.f*vI;
    *wP = -wM -2.f*wI;	
    *pP = pM;		
  }
}
    
void acousticsSurface3Dbbdg(mesh3D *mesh, int lev, dfloat time){

  // temporary storage for flux terms
  dfloat *fluxu = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxv = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxw = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxp = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));

  dfloat *fluxu_copy = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxv_copy = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxw_copy = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxp_copy = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));

  // for all elements
  for(int et=0;et<mesh->MRABNelements[lev];++et){
    int e = mesh->MRABelementIds[lev][et];
    int N = mesh->N[e];

    // for all face nodes of all elements
    for (int f=0;f<mesh->Nfaces;f++){
      for(int n=0;n<mesh->Nfp[N];++n){
        // load surface geofactors for this face
        int  sid = mesh->Nsgeo*(e*mesh->Nfaces+f);
        dfloat nx = mesh->sgeo[sid+NXID];
        dfloat ny = mesh->sgeo[sid+NYID];
        dfloat nz = mesh->sgeo[sid+NZID];
        dfloat sJ = mesh->sgeo[sid+SJID];
        dfloat invJ = mesh->sgeo[sid+IJID];

        // indices of negative and positive traces of face node
        int id  = e*mesh->NfpMax*mesh->Nfaces + f*mesh->NfpMax +n;
        int idM = id*mesh->Nfields;
        int idP = mesh->mapP[id]*mesh->Nfields;
        
        // load negative trace node values of q
        dfloat uM = mesh->fQM[idM+0];
        dfloat vM = mesh->fQM[idM+1];
        dfloat wM = mesh->fQM[idM+2];
        dfloat pM = mesh->fQM[idM+3];

        // load positive trace node values of q
        dfloat uP = mesh->fQP[idP+0]; 
        dfloat vP = mesh->fQP[idP+1]; 
        dfloat wP = mesh->fQP[idP+2]; 
        dfloat pP = mesh->fQP[idP+3]; 

        // find boundary type
        int boundaryType = mesh->EToB[e*mesh->Nfaces+f];
        if(boundaryType>0) {
          int idM = mesh->vmapM[id];
          boundaryConditions3D(boundaryType, time,
                   mesh->x[idM], mesh->y[idM], mesh->z[idM],
                   uM, vM, wM, pM,
                   &uP, &vP,&wP, &pP);
        }
        
        // compute (q^* - q^-)
        dfloat duS = 0.5*(uP-uM) + mesh->Lambda2*(-nx*(pP-pM));
        dfloat dvS = 0.5*(vP-vM) + mesh->Lambda2*(-ny*(pP-pM));
        dfloat dwS = 0.5*(wP-wM) + mesh->Lambda2*(-nz*(pP-pM));
        dfloat dpS = 0.5*(pP-pM) + mesh->Lambda2*(-nx*(uP-uM)-ny*(vP-vM)-nz*(wP-wM));
        
        // evaluate "flux" terms: (sJ/J)*(A*nx+B*ny)*(q^* - q^-)
        fluxu[f*mesh->Nfp[N]+n] = invJ*sJ*(-nx*dpS);
        fluxv[f*mesh->Nfp[N]+n] = invJ*sJ*(-ny*dpS);
        fluxw[f*mesh->Nfp[N]+n] = invJ*sJ*(-nz*dpS);
        fluxp[f*mesh->Nfp[N]+n] = invJ*sJ*(-nx*duS-ny*dvS-nz*dwS);
      }
    }

    // apply L0 to fluxes.
    for(int n=0;n<mesh->Nfp[N]*mesh->Nfaces;++n){
    
      int id = n%mesh->Nfp[N];  
      int f  = n/mesh->Nfp[N];

      dfloat utmpflux = 0.0;
      dfloat vtmpflux = 0.0;
      dfloat wtmpflux = 0.0;
      dfloat ptmpflux = 0.0;

      // sparse application of L0
      for (int m = 0; m < 7; ++m){
        int   L0id  = mesh->L0ids [N][7*id+m];
        dfloat L0val = mesh->L0vals[N][7*id+m];
        
        utmpflux += L0val * fluxu[L0id+f*mesh->Nfp[N]];
        vtmpflux += L0val * fluxv[L0id+f*mesh->Nfp[N]];
        wtmpflux += L0val * fluxw[L0id+f*mesh->Nfp[N]];
        ptmpflux += L0val * fluxp[L0id+f*mesh->Nfp[N]];
      }

      fluxu_copy[n] = utmpflux;
      fluxv_copy[n] = vtmpflux;
      fluxw_copy[n] = wtmpflux;
      fluxp_copy[n] = ptmpflux;
    }

    // apply lift reduction and accumulate RHS
    for(int n=0;n<mesh->Np[N];++n){
      int id = 3*mesh->Nfields*(mesh->NpMax*e + n) + mesh->Nfields*mesh->MRABshiftIndex[lev];
      
      // load RHS
      dfloat rhsqnu = mesh->rhsq[id+0];
      dfloat rhsqnv = mesh->rhsq[id+1];
      dfloat rhsqnw = mesh->rhsq[id+2];
      dfloat rhsqnp = mesh->rhsq[id+3];

      // sparse application of EL
      for (int m = 0; m < mesh->max_EL_nnz[N]; ++m){
        int id = m + n*mesh->max_EL_nnz[N];
        dfloat ELval = mesh->ELvals[N][id];
        int   ELid  = mesh->ELids [N][id];

        rhsqnu += ELval * fluxu_copy[ELid];
        rhsqnv += ELval * fluxv_copy[ELid];
        rhsqnw += ELval * fluxw_copy[ELid];
        rhsqnp += ELval * fluxp_copy[ELid];
      }

      // store incremented rhs
      mesh->rhsq[id+0] = rhsqnu;
      mesh->rhsq[id+1] = rhsqnv;
      mesh->rhsq[id+2] = rhsqnw;
      mesh->rhsq[id+3] = rhsqnp;

    }
  }

  // free temporary flux storage
  free(fluxu); free(fluxv); free(fluxw); free(fluxp);
  free(fluxu_copy); free(fluxv_copy); free(fluxw_copy); free(fluxp_copy);
}
