#include <math.h>
#include <stdlib.h>
#include "mesh2D.h"

void boundaryConditions2D(int bc, dfloat t, dfloat x, dfloat y,
			  dfloat q1M, dfloat q2M, dfloat q3M,
			  dfloat q4M, dfloat q5M, dfloat q6M,
			  dfloat *q1P, dfloat *q2P, dfloat *q3P,
			  dfloat *q4P, dfloat *q5P, dfloat *q6P){
  if(bc==1){
    // wall
    *q1P = q1M;
    *q2P = -q2M;
    *q3P = -q3M;
    *q4P = q4M;
    *q5P = q5M;
    *q6P = q6M;
  }		
  if(bc==2){	
    // uniform flow
  }
}


// function to compute surface contributions 
// for nodal DG boltzmann right hand side
void meshBoltzmannSurface2D(mesh2D *mesh, dfloat t){

  // temporary storage for flux terms
  dfloat *fluxq1 = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxq2 = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxq3 = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxq4 = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxq5 = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));
  dfloat *fluxq6 = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces,sizeof(dfloat));

  // for all elements
  for(int e=0;e<mesh->Nelements;++e){
    // for all face nodes of all elements
    for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      // find face that owns this node
      int face = n/mesh->Nfp;

      // load surface geofactors for this face
      int sid = mesh->Nsgeo*(e*mesh->Nfaces+face);
      dfloat nx = mesh->sgeo[sid+0];
      dfloat ny = mesh->sgeo[sid+1];
      dfloat sJ = mesh->sgeo[sid+2];
      dfloat invJ = mesh->sgeo[sid+3];

      // indices of negative and positive traces of face node
      int id  = e*mesh->Nfp*mesh->Nfaces + n;
      int idM = mesh->vmapM[id];
      int idP = mesh->vmapP[id];
      if(idP<0) idP=idM;
      int qidM = mesh->Nfields*idM;
      int qidP = mesh->Nfields*idP;
      
      // load negative trace node values of q
      dfloat q1M = mesh->q[qidM+0];
      dfloat q2M = mesh->q[qidM+1];
      dfloat q3M = mesh->q[qidM+2];
      dfloat q4M = mesh->q[qidM+3];
      dfloat q5M = mesh->q[qidM+4];
      dfloat q6M = mesh->q[qidM+5];

      // load positive trace node values of q
      dfloat q1P = mesh->q[qidP+0];
      dfloat q2P = mesh->q[qidP+1];
      dfloat q3P = mesh->q[qidP+2];
      dfloat q4P = mesh->q[qidP+3];
      dfloat q5P = mesh->q[qidP+4];
      dfloat q6P = mesh->q[qidP+5];

      // find boundary type
      int boundaryType = mesh->EToB[e*mesh->Nfaces+face];
      if(boundaryType>0)
	boundaryConditions2D(boundaryType, t, mesh->x[idM], mesh->y[idM],
			     q1M, q2M, q3M, q4M, q5M, q6M,
			     &q1P, &q2P, &q3P, &q4P, &q5P, &q6P);

      // compute (q^* - q^-)
      dfloat dq1S, dq2S, dq3S, dq4S, dq5S, dq6S;
      dfloat sqrt2 = sqrt(2.), sqrtRT = mesh->sqrtRT;
      
      dq1S = 0.5f*(q1P-q1M)
	+ sqrtRT*mesh->Lambda2*(-nx*(q2P-q2M)-ny*(q3P-q3M));   
      dq2S = 0.5f*(q2P-q2M)
	+ sqrtRT*mesh->Lambda2*(-nx*(q1P-q1M+sqrt2*(q5P-q5M))-ny*(q4P-q4M));		    
      dq3S = 0.5f*(q3P-q3M)
	+ sqrtRT*mesh->Lambda2*(-nx*(q4P-q4M)-ny*(q1P-q1M+sqrt2*(q6P-q6M)));

      dq4S = 0.5f*(q4P-q4M)
	+ sqrtRT*mesh->Lambda2*(-nx*(q3P-q3M)-ny*(q2P-q2M));
      
      dq5S = 0.5f*(q5P-q5M)
	+ sqrtRT*mesh->Lambda2*(-nx*sqrt2*(q2P-q2M));
      
      dq6S = 0.5f*(q6P-q6M)
	+ sqrtRT*mesh->Lambda2*(-ny*sqrt2*(q3P-q3M));          
      
      // evaluate "flux" terms: (sJ/J)*(A*nx+B*ny)*(q^* - q^-)
      fluxq1[n] = invJ*sJ*sqrtRT*(-nx*dq2S-ny*dq3S);
      fluxq2[n] = invJ*sJ*sqrtRT*(-nx*(dq1S+sqrt2*dq5S)-ny*dq4S);
      fluxq3[n] = invJ*sJ*sqrtRT*(-nx*dq4S-ny*(dq1S+sqrt2*dq6S));
      fluxq4[n] = invJ*sJ*sqrtRT*(-nx*dq3S-ny*dq2S);
      fluxq5[n] = invJ*sJ*sqrtRT*(-nx*sqrt2*dq2S);
      fluxq6[n] = invJ*sJ*sqrtRT*(-ny*sqrt2*dq3S);
    }
    
    // for each node in the element 
    for(int n=0;n<mesh->Np;++n){
      int id = mesh->Nfields*(mesh->Np*e + n);

      // load rhs data from volume fluxes
      dfloat rhsq1 = mesh->rhsq[id];
      dfloat rhsq2 = mesh->rhsq[id+1];
      dfloat rhsq3 = mesh->rhsq[id+2];
      dfloat rhsq4 = mesh->rhsq[id+3];
      dfloat rhsq5 = mesh->rhsq[id+4];
      dfloat rhsq6 = mesh->rhsq[id+5];
      
      // rhs += LIFT*((sJ/J)*(A*nx+B*ny)*(q^* - q^-))
      for(int m=0;m<mesh->Nfp*mesh->Nfaces;++m){
	dfloat L = mesh->LIFT[n*mesh->Nfp*mesh->Nfaces+m];
	rhsq1 += L*fluxq1[m];
	rhsq2 += L*fluxq2[m];
	rhsq3 += L*fluxq3[m];
	rhsq4 += L*fluxq4[m];
	rhsq5 += L*fluxq5[m];
	rhsq6 += L*fluxq6[m];
      }

      // store incremented rhs
      mesh->rhsq[id]   = rhsq1;
      mesh->rhsq[id+1] = rhsq2;
      mesh->rhsq[id+2] = rhsq3;
      mesh->rhsq[id+3] = rhsq4;
      mesh->rhsq[id+4] = rhsq5;
      mesh->rhsq[id+5] = rhsq6;
    }
  }

  // free temporary flux storage
  free(fluxq1); free(fluxq2); free(fluxq3);
  free(fluxq4); free(fluxq5); free(fluxq6);
}
    
