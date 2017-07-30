#include "acoustics3D.h"

void boundaryConditions3D(iint bc, dfloat time, dfloat x, dfloat y, dfloat z,
			  dfloat uM, dfloat vM, dfloat wM, dfloat pM,
			  dfloat *uP, dfloat *vP, dfloat *wP, dfloat *pP){
  if(bc==1){
    *uP = -uM;
    *vP = -vM;
    *wP = -wM;
    *pP = pM;
  }
  if(bc==2){
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


void acousticsSurface3Dbbdg(mesh3D *mesh, iint lev, dfloat time){

  // temporary storage for flux terms
  dfloat *fluxu = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxv = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxw = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxp = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));

  dfloat *fluxu_copy = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxv_copy = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxw_copy = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxp_copy = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));

  // for all elements
  for(iint et=0;et<mesh->MRABNelements[lev];++et){
    iint e = mesh->MRABelementIds[lev][et];
    // for all face nodes of all elements
    for(iint n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      // find face that owns this node
      iint face = n/mesh->Nfp;

      // load surface geofactors for this face
      iint  sid = mesh->Nsgeo*(e*mesh->Nfaces+face);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat nz = mesh->sgeo[sid+NZID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat invJ = mesh->sgeo[sid+IJID];

      // indices of negative and positive traces of face node
      iint id  = e*mesh->Nfp*mesh->Nfaces + n;
      iint idM = id*mesh->Nfields;
      iint idP = mesh->mapP[id]*mesh->Nfields;

      if(idP<0) idP = idM;

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
      iint boundaryType = mesh->EToB[e*mesh->Nfaces+face];
      if(boundaryType>0)
        boundaryConditions3D(boundaryType, time,
                 mesh->x[idM], mesh->y[idM], mesh->z[idM],
                 uM, vM, wM, pM,
                 &uP, &vP,&wP, &pP);

      // compute (q^* - q^-)
      dfloat duS = 0.5*(uP-uM) + mesh->Lambda2*(-nx*(pP-pM));
      dfloat dvS = 0.5*(vP-vM) + mesh->Lambda2*(-ny*(pP-pM));
      dfloat dwS = 0.5*(wP-wM) + mesh->Lambda2*(-nz*(pP-pM));
      dfloat dpS = 0.5*(pP-pM) + mesh->Lambda2*(-nx*(uP-uM)-ny*(vP-vM)-nz*(wP-wM));

      // evaluate "flux" terms: (sJ/J)*(A*nx+B*ny)*(q^* - q^-)
      fluxu[n] = invJ*sJ*(-nx*dpS);
      fluxv[n] = invJ*sJ*(-ny*dpS);
      fluxw[n] = invJ*sJ*(-nz*dpS);
      fluxp[n] = invJ*sJ*(-nx*duS-ny*dvS-nz*dwS);
    }

    // apply L0 to fluxes.
    for(iint n=0;n<mesh->Nfp*mesh->Nfaces;++n){

      iint id = n%mesh->Nfp;
      iint f  = n/mesh->Nfp;

      dfloat utmpflux = 0.0;
      dfloat vtmpflux = 0.0;
      dfloat wtmpflux = 0.0;
      dfloat ptmpflux = 0.0;

      // sparse application of L0
      for (int m = 0; m < 7; ++m){
        iint   L0id  = mesh->L0ids [7*id+m];
        dfloat L0val = mesh->L0vals[7*id+m];

        utmpflux += L0val * fluxu[L0id+f*mesh->Nfp];
        vtmpflux += L0val * fluxv[L0id+f*mesh->Nfp];
        wtmpflux += L0val * fluxw[L0id+f*mesh->Nfp];
        ptmpflux += L0val * fluxp[L0id+f*mesh->Nfp];
      }

      fluxu_copy[n] = utmpflux;
      fluxv_copy[n] = vtmpflux;
      fluxw_copy[n] = wtmpflux;
      fluxp_copy[n] = ptmpflux;
    }

    // apply lift reduction and accumulate RHS
    for(iint n=0;n<mesh->Np;++n){
      iint id = 3*mesh->Nfields*(mesh->Np*e + n) + mesh->Nfields*mesh->MRABshiftIndex[lev];

      // load RHS
      dfloat rhsqnu = mesh->rhsq[id+0];
      dfloat rhsqnv = mesh->rhsq[id+1];
      dfloat rhsqnw = mesh->rhsq[id+2];
      dfloat rhsqnp = mesh->rhsq[id+3];

      // sparse application of EL
      for (int m = 0; m < mesh->max_EL_nnz; ++m){
        iint id = m + n*mesh->max_EL_nnz;
        dfloat ELval = mesh->ELvals[id];
        iint   ELid  = mesh->ELids [id];

        rhsqnu += ELval * fluxu_copy[ELid];
        rhsqnv += ELval * fluxv_copy[ELid];
        rhsqnw += ELval * fluxw_copy[ELid];
        rhsqnp += ELval * fluxp_copy[ELid];
      }

      // store incremented rhs
      mesh->rhsq[id]   = rhsqnu;
      mesh->rhsq[id+1] = rhsqnv;
      mesh->rhsq[id+2] = rhsqnw;
      mesh->rhsq[id+3] = rhsqnp;

    }
  }

  // free temporary flux storage
  free(fluxu); free(fluxv); free(fluxw); free(fluxp);
  free(fluxu_copy); free(fluxv_copy); free(fluxw_copy); free(fluxp_copy);
}
