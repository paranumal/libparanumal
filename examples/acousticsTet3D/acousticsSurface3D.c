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



// function to compute surface contributions 
// for nodal DG acoustics right hand side
void acousticsSurface3D(mesh3D *mesh, dfloat time){

  // temporary storage for flux terms
  dfloat *fluxu = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxv = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxw = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxp = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));

  // for all elements
  for(iint e=0;e<mesh->Nelements;++e){
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
      iint idM = mesh->Nfields*mesh->vmapM[id];
      iint idP = mesh->Nfields*mesh->vmapP[id];

      if(idP<0) idP = idM;
      
      // load negative trace node values of q
      dfloat uM = mesh->q[idM+0];
      dfloat vM = mesh->q[idM+1];
      dfloat wM = mesh->q[idM+2];
      dfloat pM = mesh->q[idM+3];

      // load positive trace node values of q
      dfloat uP = mesh->q[idP+0]; 
      dfloat vP = mesh->q[idP+1];
      dfloat wP = mesh->q[idP+2];
      dfloat pP = mesh->q[idP+3];

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

    // for each node in the element 
    for(iint n=0;n<mesh->Np;++n){
      iint id = mesh->Nfields*(mesh->Np*e + n);

      // load rhs data from volume fluxes
      dfloat rhsu = mesh->rhsq[id];
      dfloat rhsv = mesh->rhsq[id+1];
      dfloat rhsw = mesh->rhsq[id+2];
      dfloat rhsp = mesh->rhsq[id+3];
      
      // rhs += LIFT*((sJ/J)*(A*nx+B*ny)*(q^* - q^-))
      for(int m=0;m<(mesh->Nfp*mesh->Nfaces);++m){
	dfloat L = mesh->LIFT[n*mesh->Nfp*mesh->Nfaces+m];
	rhsu += L*fluxu[m];
	rhsv += L*fluxv[m];
	rhsw += L*fluxw[m];
	rhsp += L*fluxp[m];
      }

      // store incremented rhs
      mesh->rhsq[id]   = rhsu;
      mesh->rhsq[id+1] = rhsv;
      mesh->rhsq[id+2] = rhsw;
      mesh->rhsq[id+3] = rhsp;

    }
  }

  // free temporary flux storage
  free(fluxu); free(fluxv); free(fluxw); free(fluxp);
}
    
