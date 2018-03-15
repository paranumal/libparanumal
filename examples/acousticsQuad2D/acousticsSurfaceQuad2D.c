#include <stdlib.h>
#include "acousticsQuad2D.h"

// function to compute surface contributions 
// for nodal DG acoustics right hand side
void acousticsSurfaceQuad2D(mesh2D *mesh){

  dfloat LIFT = mesh->N*(mesh->N+1.)/2.;
  
  // for all elements
  for(int e=0;e<mesh->Nelements;++e){
    // for all face nodes of all elements
    for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      // load surface geofactors for this face
      int sid = mesh->Nsgeo*(e*mesh->Nfaces*mesh->Nfp+n);
      dfloat nx = mesh->sgeo[sid+0];
      dfloat ny = mesh->sgeo[sid+1];
      dfloat sJ = mesh->sgeo[sid+2];
      dfloat invJ = mesh->sgeo[sid+3];

      // indices of negative and positive traces of face node
      int id  = e*mesh->Nfp*mesh->Nfaces + n;
      int idM = mesh->Nfields*mesh->vmapM[id];
      int idP = mesh->Nfields*mesh->vmapP[id];
      
      // load negative trace node values of q
      dfloat uM = mesh->q[idM+0];
      dfloat vM = mesh->q[idM+1];
      dfloat pM = mesh->q[idM+2];
      
      // load positive trace node values of q
      dfloat uP = mesh->q[idP+0]; // fix BCs later
      dfloat vP = mesh->q[idP+1];
      dfloat pP = mesh->q[idP+2];

      // boundary conditions 
      if(idM==idP || idP<0){
	// assert Neumann for pressure and no penetration for velocity
	uP = -uM;
	vP = -vM;
	pP =  pM;
      }
      
      // compute (q^* - q^-)
      dfloat duS = 0.5*(uP-uM) + mesh->Lambda2*(-nx)*(pP-pM);
      dfloat dvS = 0.5*(vP-vM) + mesh->Lambda2*(-ny)*(pP-pM);
      dfloat dpS = 0.5*(pP-pM) + mesh->Lambda2*(-nx*(uP-uM)-ny*(vP-vM));
      
      // evaluate "flux" terms: (sJ/J)*(A*nx+B*ny)*(q^* - q^-)
      
      mesh->rhsq[idM+0] += LIFT*invJ*sJ*(-nx*dpS);
      mesh->rhsq[idM+1] += LIFT*invJ*sJ*(-ny*dpS);
      mesh->rhsq[idM+2] += LIFT*invJ*sJ*(-nx*duS-ny*dvS);
    }
  }
}

